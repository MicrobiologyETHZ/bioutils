#!/usr/bin/env python3
"""
Unified Literature Monitoring Tool

This script combines PubMed and bioRxiv scraping with source-aware scoring
for automated literature monitoring in microbiome research.
"""

import requests
import pandas as pd
import numpy as np
from datetime import datetime, timedelta
import re
import json
import time
from concurrent.futures import ThreadPoolExecutor
import argparse
import os
import sys
import xml.etree.ElementTree as ET
from urllib.parse import quote_plus
from abc import ABC, abstractmethod

try:
    import yaml
except ImportError:
    print("ERROR: PyYAML is required for configuration files.")
    print("Please install it with: pip install pyyaml")
    sys.exit(1)


class LiteratureSource(ABC):
    """Abstract base class for literature sources"""
    
    def __init__(self, config):
        self.config = config
        self.today = datetime.now()
        self.start_date = self.today - timedelta(days=config["days_back"])
        
    @abstractmethod
    def search_papers(self, query):
        """Search for papers matching the query"""
        pass
        
    @abstractmethod
    def score_papers(self, papers_df):
        """Apply source-specific scoring to papers"""
        pass


class PubMedSource(LiteratureSource):
    """PubMed literature source using NCBI E-utilities"""
    
    def __init__(self, config):
        super().__init__(config)
        self.base_esearch = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        self.base_efetch = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        self.date_from = self.start_date.strftime("%Y/%m/%d")
        self.date_to = self.today.strftime("%Y/%m/%d")
    
    def search_papers(self, query):
        """Search PubMed for papers using NCBI E-utilities"""
        print(f"Searching PubMed for '{query}'...")
        
        search_query = f"{query} AND {self.date_from}:{self.date_to}[Date - Publication]"
        
        params = {
            "db": "pubmed",
            "term": search_query,
            "retmax": self.config.get("result_limit", 100),
            "usehistory": "y",
            "retmode": "json"
        }
        
        if self.config.get("ncbi_api_key"):
            params["api_key"] = self.config["ncbi_api_key"]
        
        max_retries = 3
        for attempt in range(max_retries):
            try:
                response = requests.get(self.base_esearch, params=params, timeout=30)
                response.raise_for_status()
                
                search_data = response.json()
                
                if 'esearchresult' in search_data and 'idlist' in search_data['esearchresult']:
                    id_list = search_data['esearchresult']['idlist']
                    
                    if not id_list:
                        print(f"No results found for '{query}'")
                        return pd.DataFrame()
                    
                    print(f"Found {len(id_list)} papers for '{query}'")
                    
                    web_env = search_data['esearchresult']['webenv']
                    query_key = search_data['esearchresult']['querykey']
                    
                    return self._fetch_paper_details(id_list, web_env, query_key, query)
                    
                else:
                    print(f"No results or unexpected response format for '{query}'")
                    return pd.DataFrame()
                
            except requests.exceptions.RequestException as e:
                if attempt < max_retries - 1:
                    wait_time = 2 ** attempt
                    print(f"API request failed, retrying in {wait_time} seconds...")
                    time.sleep(wait_time)
                else:
                    print(f"Error searching PubMed for '{query}': {e}")
                    return pd.DataFrame()
        
        return pd.DataFrame()
    
    def _fetch_paper_details(self, id_list, web_env, query_key, original_query):
        """Fetch detailed paper information using EFetch"""
        params = {
            "db": "pubmed",
            "WebEnv": web_env,
            "query_key": query_key,
            "retmode": "xml",
            "rettype": "abstract"
        }
        
        if self.config.get("ncbi_api_key"):
            params["api_key"] = self.config["ncbi_api_key"]
        
        batch_size = 50
        all_results = []
        
        for start in range(0, len(id_list), batch_size):
            end = min(start + batch_size, len(id_list))
            batch_ids = id_list[start:end]
            
            params["id"] = ",".join(batch_ids)
            
            max_retries = 3
            for attempt in range(max_retries):
                try:
                    response = requests.get(self.base_efetch, params=params, timeout=30)
                    response.raise_for_status()
                    
                    papers = self._parse_pubmed_xml(response.text, original_query)
                    all_results.extend(papers)
                    
                    time.sleep(0.5)
                    break
                    
                except requests.exceptions.RequestException as e:
                    if attempt < max_retries - 1:
                        wait_time = 2 ** attempt
                        print(f"API request failed, retrying in {wait_time} seconds...")
                        time.sleep(wait_time)
                    else:
                        print(f"Error fetching paper details: {e}")
        
        if all_results:
            df = pd.DataFrame(all_results)
            df['source'] = 'PubMed'
            return df
        else:
            return pd.DataFrame()
    
    def _parse_pubmed_xml(self, xml_text, original_query):
        """Parse PubMed XML response to extract paper details"""
        papers = []
        
        try:
            root = ET.fromstring(xml_text)
            
            for article in root.findall(".//PubmedArticle"):
                try:
                    pmid = article.find(".//PMID").text if article.find(".//PMID") is not None else ""
                    
                    title_element = article.find(".//ArticleTitle")
                    title = title_element.text if title_element is not None else "No title"
                    
                    abstract_elements = article.findall(".//AbstractText")
                    abstract = " ".join([elem.text for elem in abstract_elements if elem.text]) if abstract_elements else ""
                    
                    author_list = article.findall(".//Author")
                    authors = []
                    
                    for author in author_list:
                        last_name = author.find("LastName")
                        fore_name = author.find("ForeName")
                        initials = author.find("Initials")
                        
                        if last_name is not None:
                            if fore_name is not None:
                                authors.append(f"{last_name.text} {fore_name.text}")
                            elif initials is not None:
                                authors.append(f"{last_name.text} {initials.text}")
                            else:
                                authors.append(last_name.text)
                    
                    author_string = ", ".join(authors)
                    
                    pub_date = article.find(".//PubDate")
                    year = pub_date.find("Year").text if pub_date.find("Year") is not None else ""
                    month = pub_date.find("Month").text if pub_date.find("Month") is not None else ""
                    day = pub_date.find("Day").text if pub_date.find("Day") is not None else ""
                    
                    date = f"{year}-{month}-{day}" if day else (f"{year}-{month}" if month else year)
                    
                    journal_element = article.find(".//Journal/Title")
                    journal = journal_element.text if journal_element is not None else ""
                    
                    doi_element = article.find(".//ArticleId[@IdType='doi']")
                    doi = doi_element.text if doi_element is not None else ""
                    
                    url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else ""
                    
                    paper = {
                        'title': title,
                        'authors': author_string,
                        'date': date,
                        'pmid': pmid,
                        'doi': doi,
                        'abstract': abstract,
                        'journal': journal,
                        'url': url,
                        'search_term': original_query
                    }
                    
                    papers.append(paper)
                    
                except Exception as e:
                    print(f"Error parsing article: {e}")
                    continue
            
        except Exception as e:
            print(f"Error parsing XML: {e}")
        
        return papers
    
    def score_papers(self, papers_df):
        """Apply PubMed-specific scoring with configurable weights"""
        if papers_df.empty:
            return papers_df
        
        # Get scoring weights from config
        scoring = self.config.get("scoring", {}).get("pubmed", {})
        
        # Initialize base score
        papers_df['score'] = 0
        
        # Ensure string columns
        string_columns = ['title', 'authors', 'abstract', 'journal']
        for col in string_columns:
            if col in papers_df.columns:
                papers_df[col] = papers_df[col].fillna('').astype(str)
            else:
                papers_df[col] = ''
        
        # System relevance
        system_pattern = '|'.join(self.config["relevant_systems"])
        papers_df['system_match_title'] = papers_df['title'].str.lower().str.contains(
            system_pattern, regex=True, case=False, na=False).astype(int) * scoring.get("system_relevance_title", 3)
        papers_df['system_match_abstract'] = papers_df['abstract'].str.lower().str.contains(
            system_pattern, regex=True, case=False, na=False).astype(int) * scoring.get("system_relevance_abstract", 2)
        papers_df['system_match'] = papers_df[['system_match_title', 'system_match_abstract']].max(axis=1)
        papers_df['score'] += papers_df['system_match']
        
        # Computational methods
        comp_pattern = '|'.join(self.config["computational_keywords"])
        papers_df['computational_match_title'] = papers_df['title'].str.lower().str.contains(
            comp_pattern, regex=True, case=False, na=False).astype(int) * scoring.get("computational_title", 2)
        papers_df['computational_match_abstract'] = papers_df['abstract'].str.lower().str.contains(
            comp_pattern, regex=True, case=False, na=False).astype(int) * scoring.get("computational_abstract", 1)
        papers_df['computational_match'] = papers_df[['computational_match_title', 'computational_match_abstract']].max(axis=1)
        papers_df['score'] += papers_df['computational_match']
        
        # PubMed-specific: Journal prestige
        top_journals = self.config.get("top_tier_journals", ["Nature", "Science", "Cell"])
        good_journals = self.config.get("important_journals", [])
        
        papers_df['top_journal'] = papers_df['journal'].str.contains(
            '|'.join(top_journals), regex=True, case=False, na=False).astype(int) * scoring.get("top_tier_journal", 3)
        papers_df['good_journal'] = papers_df['journal'].str.contains(
            '|'.join(good_journals), regex=True, case=False, na=False).astype(int) * scoring.get("important_journal", 2)
        papers_df['journal_score'] = papers_df[['top_journal', 'good_journal']].max(axis=1)
        papers_df['score'] += papers_df['journal_score']
        
        # Important authors
        author_pattern = '|'.join(self.config["important_authors"])
        papers_df['author_match'] = papers_df['authors'].str.contains(
            author_pattern, regex=True, case=False, na=False).astype(int) * scoring.get("important_author", 1)
        papers_df['score'] += papers_df['author_match']
        
        # Dataset availability
        if not papers_df['abstract'].str.len().eq(0).all():
            papers_df['dataset_available'] = papers_df['abstract'].str.lower().str.contains(
                'dataset|data available|github|available at|repository', regex=True, na=False).astype(int) * scoring.get("dataset_available", 1)
            papers_df['score'] += papers_df['dataset_available']
        
        # Review penalty/bonus
        papers_df['is_review'] = papers_df['title'].str.lower().str.contains(
            'review|overview|survey', regex=True, na=False).astype(int)
        papers_df['comprehensive_review'] = papers_df['title'].str.lower().str.contains(
            'comprehensive|systematic', regex=True, na=False).astype(int)
        
        # Apply review penalty
        papers_df['review_penalty'] = papers_df['is_review'].astype(int) * scoring.get("review_penalty", -1)
        papers_df['score'] += papers_df['review_penalty']
        
        # Apply comprehensive review bonus (cancels penalty + adds bonus)
        papers_df['comprehensive_bonus'] = papers_df['comprehensive_review'].astype(int) * (abs(scoring.get("review_penalty", -1)) + scoring.get("comprehensive_review_bonus", 1))
        papers_df['score'] += papers_df['comprehensive_bonus']
        
        return papers_df.sort_values('score', ascending=False)


class BioRxivSource(LiteratureSource):
    """bioRxiv literature source using web scraping"""
    
    def __init__(self, config):
        super().__init__(config)
    
    def search_papers(self, query):
        """Search bioRxiv using their official API"""
        date_str = self.start_date.strftime("%Y-%m-%d")
        today_str = self.today.strftime("%Y-%m-%d")
        
        print(f"Searching bioRxiv API for papers from {date_str} to {today_str}")
        
        # bioRxiv API endpoint - gets ALL papers in date range
        api_url = f"https://api.biorxiv.org/details/biorxiv/{date_str}/{today_str}"
        
        max_retries = 3
        for attempt in range(max_retries):
            try:
                print(f"Making API request: {api_url}")
                response = requests.get(api_url, timeout=60)  # Longer timeout for API
                
                if response.status_code == 200:
                    try:
                        data = response.json()
                        
                        if not data.get('messages') or data['messages'][0].get('status') != 'ok':
                            print(f"API returned error status for date range {date_str} to {today_str}")
                            return pd.DataFrame()
                        
                        if not data.get('collection'):
                            print(f"No papers found in date range {date_str} to {today_str}")
                            return pd.DataFrame()
                        
                        print(f"Found {len(data['collection'])} total papers, filtering by query '{query}'")
                        
                        # Filter papers by query terms
                        results = []
                        query_terms = query.lower().replace('+', ' ').replace('AND', '').split()
                        query_terms = [term.strip() for term in query_terms if term.strip()]
                        
                        for paper in data['collection']:
                            try:
                                title = paper.get('title', '').lower()
                                abstract = paper.get('abstract', '').lower()
                                authors = paper.get('authors', '').lower()
                                
                                # Check if any query term appears in title, abstract, or authors
                                text_to_search = f"{title} {abstract} {authors}"
                                
                                if any(term in text_to_search for term in query_terms):
                                    paper_data = {
                                        'title': paper.get('title', ''),
                                        'authors': paper.get('authors', ''),
                                        'date': paper.get('date', ''),
                                        'doi': paper.get('doi', ''),
                                        'abstract': paper.get('abstract', ''),
                                        'journal': 'bioRxiv',
                                        'url': f"https://www.biorxiv.org/content/10.1101/{paper.get('doi', '')}" if paper.get('doi') else '',
                                        'search_term': query
                                    }
                                    results.append(paper_data)
                                    
                            except Exception as e:
                                print(f"Error processing paper: {e}")
                                continue
                        
                        print(f"Found {len(results)} papers matching query '{query}'")
                        
                        if results:
                            df = pd.DataFrame(results)
                            df['source'] = 'bioRxiv'
                            return df
                        else:
                            return pd.DataFrame()
                            
                    except json.JSONDecodeError as e:
                        print(f"Error parsing JSON response: {e}")
                        return pd.DataFrame()
                        
                else:
                    print(f"API Error: {response.status_code} - {response.text[:200]}")
                    if attempt < max_retries - 1:
                        wait_time = 5 * (attempt + 1)
                        print(f"Retrying in {wait_time} seconds...")
                        time.sleep(wait_time)
                    else:
                        return pd.DataFrame()
                        
            except requests.exceptions.RequestException as e:
                print(f"Request error: {e}")
                if attempt < max_retries - 1:
                    wait_time = 5 * (attempt + 1)
                    print(f"Retrying in {wait_time} seconds...")
                    time.sleep(wait_time)
                else:
                    return pd.DataFrame()
        
        return pd.DataFrame()
    
    def score_papers(self, papers_df):
        """Apply bioRxiv-specific scoring with configurable weights"""
        if papers_df.empty:
            return papers_df
        
        # Get scoring weights from config
        scoring = self.config.get("scoring", {}).get("biorxiv", {})
        
        # Initialize base score
        papers_df['score'] = 0
        
        # Ensure string columns
        string_columns = ['title', 'authors', 'abstract']
        for col in string_columns:
            if col in papers_df.columns:
                papers_df[col] = papers_df[col].fillna('').astype(str)
            else:
                papers_df[col] = ''
        
        # System relevance
        system_pattern = '|'.join(self.config["relevant_systems"])
        papers_df['system_match_title'] = papers_df['title'].str.lower().str.contains(
            system_pattern, regex=True, case=False, na=False).astype(int) * scoring.get("system_relevance_title", 3)
        papers_df['system_match_abstract'] = papers_df['abstract'].str.lower().str.contains(
            system_pattern, regex=True, case=False, na=False).astype(int) * scoring.get("system_relevance_abstract", 2)
        papers_df['system_match'] = papers_df[['system_match_title', 'system_match_abstract']].max(axis=1)
        papers_df['score'] += papers_df['system_match']
        
        # bioRxiv-specific: Higher weight on computational methods
        # Preprints often showcase novel computational approaches
        comp_pattern = '|'.join(self.config["computational_keywords"])
        papers_df['computational_match_title'] = papers_df['title'].str.lower().str.contains(
            comp_pattern, regex=True, case=False, na=False).astype(int) * scoring.get("computational_title", 3)
        papers_df['computational_match_abstract'] = papers_df['abstract'].str.lower().str.contains(
            comp_pattern, regex=True, case=False, na=False).astype(int) * scoring.get("computational_abstract", 2)
        papers_df['computational_match'] = papers_df[['computational_match_title', 'computational_match_abstract']].max(axis=1)
        papers_df['score'] += papers_df['computational_match']
        
        # bioRxiv-specific: Preprint novelty bonus
        # Recent preprints get boost for being "first to market"
        papers_df['preprint_bonus'] = scoring.get("preprint_novelty_bonus", 1)
        papers_df['score'] += papers_df['preprint_bonus']
        
        # Important authors (higher weight for bioRxiv since no journal signal)
        author_pattern = '|'.join(self.config["important_authors"])
        papers_df['author_match'] = papers_df['authors'].str.contains(
            author_pattern, regex=True, case=False, na=False).astype(int) * scoring.get("important_author", 2)
        papers_df['score'] += papers_df['author_match']
        
        # Dataset availability (more important for preprints)
        if not papers_df['abstract'].str.len().eq(0).all():
            papers_df['dataset_available'] = papers_df['abstract'].str.lower().str.contains(
                'dataset|data available|github|available at|repository', regex=True, na=False).astype(int) * scoring.get("dataset_available", 2)
            papers_df['score'] += papers_df['dataset_available']
        
        # Review penalty/bonus
        papers_df['is_review'] = papers_df['title'].str.lower().str.contains(
            'review|overview|survey', regex=True, na=False).astype(int)
        papers_df['comprehensive_review'] = papers_df['title'].str.lower().str.contains(
            'comprehensive|systematic', regex=True, na=False).astype(int)
        
        # Apply review penalty
        papers_df['review_penalty'] = papers_df['is_review'].astype(int) * scoring.get("review_penalty", -1)
        papers_df['score'] += papers_df['review_penalty']
        
        # Apply comprehensive review bonus (cancels penalty + adds bonus)
        papers_df['comprehensive_bonus'] = papers_df['comprehensive_review'].astype(int) * (abs(scoring.get("review_penalty", -1)) + scoring.get("comprehensive_review_bonus", 1))
        papers_df['score'] += papers_df['comprehensive_bonus']
        
        return papers_df.sort_values('score', ascending=False)


class LiteratureMonitor:
    """Unified literature monitoring system"""
    
    def __init__(self, config_file=None):
        """Initialize with default config or from a file"""
        # Default configuration with scoring weights
        self.config = {
            "search_terms": [
                "microbiome", "microbial community",
                "metatranscriptomics",
                "bacterial growth", "bacterial metabolism", 
            ],
            "combined_searches": [
                "bacterial AND gene AND expression",
                "microbiome AND transcriptomics",
                "microbiota AND growth",
                "bacterial AND transcriptomics AND computational"
            ],
            "important_authors": [
                "Knight R", "Segata N", "Huttenhower C",
                "Gilbert J", "Turnbaugh P", "Fischbach M"
            ],
            "important_journals": [
                "Nature Microbiology", "Cell Host Microbe", "Microbiome",
                "ISME Journal", "mSystems", "PLoS Comput Biol"
            ],
            "top_tier_journals": [
                "Nature", "Science", "Cell"
            ],
            "relevant_systems": [
                "soil microbiome", 
                "marine microbiome", "plant microbiome", "oral microbiome",
                "ocean microbiome", "coral microbiome", "gut microbiome"
            ],
            "computational_keywords": [
                "machine learning", "deep learning", "network analysis",
                "bayesian", "simulation", "modeling", "computational",
                "metatranscriptomics"
            ],
            "scoring": {
                "pubmed": {
                    "system_relevance_title": 3,
                    "system_relevance_abstract": 2,
                    "computational_title": 2,
                    "computational_abstract": 1,
                    "top_tier_journal": 3,
                    "important_journal": 2,
                    "important_author": 1,
                    "dataset_available": 1,
                    "review_penalty": -1,
                    "comprehensive_review_bonus": 1
                },
                "biorxiv": {
                    "system_relevance_title": 3,
                    "system_relevance_abstract": 2,
                    "computational_title": 3,
                    "computational_abstract": 2,
                    "preprint_novelty_bonus": 1,
                    "important_author": 2,
                    "dataset_available": 2,
                    "review_penalty": -1,
                    "comprehensive_review_bonus": 1
                }
            },
            "defaults": {
                "days_back": 3,
                "min_score": 3,
                "result_limit": 100,
                "sources": ["pubmed", "biorxiv"]
            },
            "api_keys": {
                "ncbi_api_key": ""
            }
        }
        
        # Load custom configuration if provided
        if config_file:
            with open(config_file, 'r') as f:
                if config_file.endswith('.yaml') or config_file.endswith('.yml'):
                    custom_config = yaml.safe_load(f)
                else:
                    # Support JSON for backward compatibility
                    custom_config = json.load(f)
                self._deep_update(self.config, custom_config)
        
        # Apply defaults for backward compatibility
        self.config.update({
            "days_back": self.config["defaults"]["days_back"],
            "min_score": self.config["defaults"]["min_score"],
            "result_limit": self.config["defaults"]["result_limit"],
            "sources": self.config["defaults"]["sources"],
            "ncbi_api_key": self.config["api_keys"]["ncbi_api_key"]
        })
        
        # Initialize sources
        self.sources = {}
        if "pubmed" in self.config["sources"]:
            self.sources["pubmed"] = PubMedSource(self.config)
        if "biorxiv" in self.config["sources"]:
            self.sources["biorxiv"] = BioRxivSource(self.config)
        
        # Results storage
        self.all_papers = pd.DataFrame()
        self.source_papers = {}  # Separate storage by source
    
    def _deep_update(self, base_dict, update_dict):
        """Recursively update nested dictionaries"""
        for key, value in update_dict.items():
            if key in base_dict and isinstance(base_dict[key], dict) and isinstance(value, dict):
                self._deep_update(base_dict[key], value)
            else:
                base_dict[key] = value
    
    def run_searches(self):
        """Run all configured searches across all sources"""
        print(f"Starting literature search across {len(self.sources)} sources...")
        
        all_source_results = []
        
        for source_name, source in self.sources.items():
            print(f"\n=== Searching {source_name.upper()} ===")
            
            # Run searches for this source
            source_dfs = []
            
            # Use ThreadPoolExecutor for parallel searches within each source
            with ThreadPoolExecutor(max_workers=3) as executor:
                # Submit individual term searches
                term_futures = {executor.submit(source.search_papers, term): term 
                               for term in self.config["search_terms"]}
                
                # Submit combined searches
                combined_futures = {executor.submit(source.search_papers, query): query 
                                   for query in self.config["combined_searches"]}
                
                # Collect results from individual terms
                for future in term_futures:
                    term = term_futures[future]
                    try:
                        df = future.result()
                        if not df.empty:
                            print(f"Found {len(df)} papers for term '{term}'")
                            source_dfs.append(df)
                    except Exception as e:
                        print(f"Error processing results for '{term}': {e}")
                
                # Collect results from combined searches
                for future in combined_futures:
                    query = combined_futures[future]
                    try:
                        df = future.result()
                        if not df.empty:
                            print(f"Found {len(df)} papers for query '{query}'")
                            source_dfs.append(df)
                    except Exception as e:
                        print(f"Error processing results for '{query}': {e}")
            
            # Combine results for this source
            if source_dfs:
                source_combined = pd.concat(source_dfs)
                
                # Deduplicate within source
                source_combined = self._deduplicate_papers(source_combined, source_name)
                
                # Apply source-specific scoring
                source_combined = source.score_papers(source_combined)
                
                # Store source-specific results
                self.source_papers[source_name] = source_combined
                all_source_results.append(source_combined)
                
                print(f"Total unique papers from {source_name}: {len(source_combined)}")
            
        # Combine all sources
        if all_source_results:
            self.all_papers = pd.concat(all_source_results, ignore_index=True)
            
            # Cross-source deduplication (more sophisticated)
            self.all_papers = self._cross_source_deduplicate(self.all_papers)
            
            print(f"\nTotal unique papers across all sources: {len(self.all_papers)}")
        else:
            print("No papers found matching your search criteria.")
    
    def _deduplicate_papers(self, papers_df, source_name):
        """Deduplicate papers within a single source"""
        if papers_df.empty:
            return papers_df
        
        # Primary deduplication key depends on source
        if source_name == "pubmed" and 'pmid' in papers_df.columns:
            # Use PMID for PubMed
            valid_ids = papers_df['pmid'].notna() & (papers_df['pmid'] != '')
            if valid_ids.any():
                papers_df = papers_df.groupby('pmid', as_index=False).agg({
                    'title': 'first',
                    'authors': 'first',
                    'date': 'first',
                    'doi': 'first',
                    'abstract': 'first',
                    'journal': 'first',
                    'url': 'first',
                    'source': 'first',
                    'search_term': lambda x: ', '.join(sorted(set(x)))
                })
        elif 'doi' in papers_df.columns:
            # Use DOI if available
            valid_dois = papers_df['doi'].notna() & (papers_df['doi'] != '')
            if valid_dois.any():
                papers_df = papers_df.groupby('doi', as_index=False).agg({
                    'title': 'first',
                    'authors': 'first',
                    'date': 'first',
                    'abstract': 'first',
                    'journal': 'first',
                    'url': 'first',
                    'source': 'first',
                    'search_term': lambda x: ', '.join(sorted(set(x)))
                })
        else:
            # Fall back to title-based deduplication
            papers_df = papers_df.groupby('title', as_index=False).agg({
                'authors': 'first',
                'date': 'first',
                'doi': 'first',
                'abstract': 'first',
                'journal': 'first',
                'url': 'first',
                'source': 'first',
                'search_term': lambda x: ', '.join(sorted(set(x)))
            })
        
        return papers_df
    
    def _cross_source_deduplicate(self, papers_df):
        """Remove duplicates across sources (same paper published in both)"""
        if papers_df.empty:
            return papers_df
        
        # Create a similarity key for cross-source matching
        # Use first few words of title + first author last name
        def create_similarity_key(row):
            title_words = str(row['title']).lower().split()[:5]
            first_author = str(row['authors']).split(',')[0].split()[-1] if row['authors'] else ""
            return ' '.join(title_words) + '_' + first_author.lower()
        
        papers_df['similarity_key'] = papers_df.apply(create_similarity_key, axis=1)
        
        # Group by similarity key and prefer PubMed over bioRxiv for published papers
        def choose_best_version(group):
            if len(group) == 1:
                return group.iloc[0]
            
            # If we have both PubMed and bioRxiv, prefer PubMed (published version)
            pubmed_papers = group[group['source'] == 'PubMed']
            if not pubmed_papers.empty:
                return pubmed_papers.iloc[0]  # Take first PubMed version
            else:
                return group.iloc[0]  # Take first version if all same source
        
        # Apply deduplication
        deduplicated = papers_df.groupby('similarity_key').apply(choose_best_version).reset_index(drop=True)
        
        # Drop the similarity key
        deduplicated = deduplicated.drop('similarity_key', axis=1)
        
        return deduplicated
    
    def get_top_papers(self, limit=20, source=None):
        """Get top papers overall or from specific source"""
        if source:
            if source in self.source_papers:
                papers = self.source_papers[source]
            else:
                return pd.DataFrame()
        else:
            papers = self.all_papers
        
        if papers.empty:
            return pd.DataFrame()
        
        # Filter by minimum score and return top N
        filtered = papers[papers['score'] >= self.config["min_score"]]
        return filtered.head(limit)
    
    def print_results(self, limit=20, source=None):
        """Print formatted results"""
        if source:
            print(f"\n=== Top {source.upper()} Papers ===")
            top_papers = self.get_top_papers(limit, source)
        else:
            print(f"\n=== Top Papers Across All Sources ===")
            top_papers = self.get_top_papers(limit)
        
        if top_papers.empty:
            print(f"No papers met your minimum score of {self.config['min_score']}.")
            return
        
        print(f"\nShowing {len(top_papers)} papers (minimum score: {self.config['min_score']})\n")
        
        for i, (_, paper) in enumerate(top_papers.iterrows(), 1):
            print(f"{i}. [{paper['source']}] {paper['title']}")
            
            # Format authors
            authors = self._format_authors(paper['authors'])
            print(f"   Authors: {authors}")
            
            # Journal info
            if paper.get('journal'):
                print(f"   Journal: {paper['journal']} (Date: {paper['date']})")
            
            # IDs
            if paper.get('pmid'):
                print(f"   PMID: {paper['pmid']}")
            if paper.get('doi'):
                print(f"   DOI: {paper['doi']}")
            
            print(f"   Score: {paper['score']}")
            
            # Abstract preview
            if paper.get('abstract') and len(str(paper['abstract'])) > 10:
                abstract_preview = str(paper['abstract'])[:150] + "..." if len(str(paper['abstract'])) > 150 else str(paper['abstract'])
                print(f"   Abstract: {abstract_preview}")
            
            print()
    
    def _format_authors(self, authors_string):
        """Format authors to show only first and last few"""
        if not authors_string:
            return "No authors"
        
        # Clean up bioRxiv author formatting
        cleaned = str(authors_string).replace("View ORCID Profile", "").replace("  ", " ").strip()
        
        # Split authors
        authors = [a.strip() for a in re.split(r',|\sand\s', cleaned) if a.strip()]
        
        if len(authors) <= 4:
            return ", ".join(authors)
        else:
            return f"{authors[0]}, ... {', '.join(authors[-2:])}"
    
    def save_results(self, filename="literature_monitor.csv"):
        """Save all results to CSV files"""
        if self.all_papers.empty:
            print("No results to save.")
            return
        
        # Save combined results
        self._save_papers_csv(self.all_papers, filename)
        
        # Save source-specific results
        for source_name, papers in self.source_papers.items():
            source_filename = filename.replace('.csv', f'_{source_name}.csv')
            self._save_papers_csv(papers, source_filename)
    
    def _save_papers_csv(self, papers_df, filename):
        """Save papers DataFrame to CSV with formatting"""
        if papers_df.empty:
            return
        
        # Create formatted authors column
        papers_df['formatted_authors'] = papers_df['authors'].apply(self._format_authors)
        
        # Define columns to save
        columns_to_save = [
            'title', 'formatted_authors', 'journal', 'date', 
            'score', 'source', 'search_term', 'doi', 'url'
        ]
        
        # Add PMID if it exists
        if 'pmid' in papers_df.columns:
            columns_to_save.insert(-2, 'pmid')
        
        # Only include existing columns
        available_columns = [col for col in columns_to_save if col in papers_df.columns]
        
        # Save full version
        papers_df[available_columns].to_csv(filename, index=False)
        print(f"Results saved to {filename}")
        
        # Save simplified version
        simple_filename = filename.replace('.csv', '_simple.csv')
        simple_columns = ['title', 'formatted_authors', 'journal', 'date', 'score', 'source']
        simple_available = [col for col in simple_columns if col in papers_df.columns]
        
        papers_df[simple_available].to_csv(simple_filename, index=False)
        print(f"Simplified results saved to {simple_filename}")
    
    def save_config(self, filename="literature_monitor_config.yaml"):
        """Save current configuration"""
        with open(filename, 'w') as f:
            if filename.endswith('.yaml') or filename.endswith('.yml'):
                yaml.dump(self.config, f, default_flow_style=False, indent=2)
            else:
                # JSON fallback
                json.dump(self.config, f, indent=2)
        print(f"Configuration saved to {filename}")


def create_sample_config():
    """Create a sample YAML configuration file with explanations"""
    
    # Create the YAML content as a string with comments
    yaml_content = """# Literature Monitor Configuration File
# Customize this file for your research area and preferences

# Search terms - individual keywords to search for
search_terms:
  - microbiome
  - microbial community
  - metatranscriptomics
  - bacterial growth
  - bacterial metabolism

# Combined searches - more complex queries with AND/OR logic
combined_searches:
  - bacterial AND gene AND expression
  - microbiome AND transcriptomics
  - microbiota AND growth

# Important authors you want to follow (use "LastName Initial" format)
important_authors:
  - Knight R
  - Segata N
  - Huttenhower C
  - Gilbert J
  - Turnbaugh P

# Important journals (will get scoring bonus)
important_journals:
  - Nature Microbiology
  - Cell Host Microbe
  - Microbiome
  - ISME Journal
  - mSystems
  - PLoS Comput Biol

# Top-tier journals (get highest scoring bonus)
top_tier_journals:
  - Nature
  - Science
  - Cell

# Relevant biological systems you're interested in
relevant_systems:
  - soil microbiome
  - marine microbiome
  - gut microbiome
  - plant microbiome
  - oral microbiome

# Computational/methodological keywords
computational_keywords:
  - machine learning
  - deep learning
  - network analysis
  - modeling
  - computational
  - bioinformatics

# Scoring weights - adjust these to emphasize what's important to you
# Higher positive numbers = more important, negative = penalty
scoring:
  # PubMed scoring (published papers) - journal prestige matters here
  pubmed:
    system_relevance_title: 3      # System keywords in title
    system_relevance_abstract: 2   # System keywords in abstract
    computational_title: 2         # Computational keywords in title
    computational_abstract: 1      # Computational keywords in abstract
    top_tier_journal: 3           # Published in Nature/Science/Cell
    important_journal: 2          # Published in field-specific journals
    important_author: 1           # Paper by authors you follow
    dataset_available: 1          # Data/code availability mentioned
    review_penalty: -1            # Penalty for review papers
    comprehensive_review_bonus: 1  # Bonus for systematic reviews

  # bioRxiv scoring (preprints) - computational methods weighted higher
  biorxiv:
    system_relevance_title: 3      # System keywords in title
    system_relevance_abstract: 2   # System keywords in abstract
    computational_title: 3         # Computational keywords in title (higher than PubMed)
    computational_abstract: 2      # Computational keywords in abstract
    preprint_novelty_bonus: 1      # Bonus for being first-to-market
    important_author: 2           # Higher weight since no journal signal
    dataset_available: 2          # More important for preprints
    review_penalty: -1            # Penalty for review papers
    comprehensive_review_bonus: 1  # Bonus for systematic reviews

# Default settings
defaults:
  days_back: 3                    # Days to look back (can be overridden with --days)
  min_score: 3                    # Minimum score to include in results
  result_limit: 100               # Maximum papers to fetch per search term
  sources:                        # Which sources to search
    - pubmed
    - biorxiv

# API keys for faster access
api_keys:
  ncbi_api_key: ""               # Add your NCBI API key here for faster PubMed searches
"""
    
    with open("literature_config.yaml", 'w') as f:
        f.write(yaml_content)
    
    print("Sample YAML configuration created: literature_config.yaml")
    print("\nCustomize this file to:")
    print("• Add your research keywords to 'search_terms'")
    print("• List authors you follow in 'important_authors'")
    print("• Adjust scoring weights in 'scoring' section")
    print("• Set your NCBI API key for faster PubMed searches")
    print("\nThen use: --config literature_config.yaml")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Unified Literature Monitoring Tool')
    parser.add_argument('--config', help='Path to configuration file')
    parser.add_argument('--days', type=int, help='Number of days to look back')
    parser.add_argument('--limit', type=int, default=20, help='Maximum papers to return')
    parser.add_argument('--min-score', type=int, help='Minimum score threshold')
    parser.add_argument('--sources', nargs='+', choices=['pubmed', 'biorxiv'], 
                        help='Which sources to search (default: both)')
    parser.add_argument('--output', help='Output CSV filename')
    parser.add_argument('--save-config', action='store_true', help='Save current config to JSON')
    parser.add_argument('--create-config', action='store_true', help='Create sample configuration file')
    parser.add_argument('--api-key', help='NCBI API key for PubMed')
    
    args = parser.parse_args()
    
    # Create sample config if requested
    if args.create_config:
        create_sample_config()
        sys.exit(0)
    
    # Create monitor
    monitor = LiteratureMonitor(config_file=args.config)
    
    # Override config with command line arguments
    if args.days:
        monitor.config["days_back"] = args.days
    if args.min_score:
        monitor.config["min_score"] = args.min_score
    if args.sources:
        monitor.config["sources"] = args.sources
    if args.api_key:
        monitor.config["ncbi_api_key"] = args.api_key
    
    # Reinitialize sources with updated config
    monitor.sources = {}
    if "pubmed" in monitor.config["sources"]:
        monitor.sources["pubmed"] = PubMedSource(monitor.config)
    if "biorxiv" in monitor.config["sources"]:
        monitor.sources["biorxiv"] = BioRxivSource(monitor.config)
    
    # Run the monitoring
    monitor.run_searches()
    monitor.print_results(limit=args.limit)
    
    # Save results if requested
    if args.output:
        monitor.save_results(args.output)
    
    # Save config if requested
    if args.save_config:
        monitor.save_config()