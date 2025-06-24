#!/usr/bin/env python3
"""
PubMed Paper Scraper and Prioritizer

This script searches for and prioritizes papers from PubMed based on user-specified
criteria related to microbiome research, bacterial growth, and transcriptomics.
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

class PubMedScraper:
    """
    A class to search and prioritize papers from PubMed based on user-defined criteria.
    """
    
    def __init__(self, config_file=None):
        """Initialize with default config or from a file"""
        # Default configuration
        self.config = {
            "search_terms": [
                "microbiome",  "microbial community",
                "metatranscriptomics",
                "bacterial growth", "bacterial metabolism", 
            ],
            "combined_searches": [
                "bacterial AND gene AND expression"
                "microbiome AND transcriptomics",
                "microbiota AND growth",
                "bacterial AND transcriptomics AND computational"
            ],
            "important_authors": [
                # Add authors you follow closely
                "Knight R", "Segata N", "Huttenhower C",
                "Gilbert J", "Turnbaugh P", "Fischbach M"
            ],
            "important_journals": [
                "Nature Microbiology", "Cell Host Microbe", "Microbiome",
                "ISME Journal", "mSystems", "PLoS Comput Biol",
                "Nature", "Science", "Cell"
            ],
            "relevant_systems": [
                # Your specific systems of interest
                "soil microbiome", 
                "marine microbiome", "plant microbiome", "oral microbiome",
                "ocean microbiome", 'coral microbiome'
            ],
            "computational_keywords": [
                "machine learning", "deep learning", "network analysis",
                "bayesian", "simulation", "modeling", "computational",
                "metatranscriptomics"
                #"bioinformatics", "algorithm", "pipeline"
            ],
            "days_back": 3,  # PubMed often has a longer publication cycle
            "min_score": 3,
            "result_limit": 100,  # Number of results to retrieve per query
            "ncbi_api_key": ""  # Optional: add your NCBI API key here
        }
        
        # Load custom configuration if provided
        if config_file:
            with open(config_file, 'r') as f:
                custom_config = json.load(f)
                self.config.update(custom_config)
        
        # Calculate date range for searches
        self.today = datetime.now()
        self.start_date = self.today - timedelta(days=self.config["days_back"])
        
        # Format dates for PubMed
        self.date_from = self.start_date.strftime("%Y/%m/%d")
        self.date_to = self.today.strftime("%Y/%m/%d")
        
        # Prepare results storage
        self.all_papers = pd.DataFrame()
        
        # Base URLs for NCBI E-utilities
        self.base_esearch = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        self.base_efetch = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        self.base_esummary = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    
    def search_pubmed(self, query):
        """Search PubMed for a specific query within date range using NCBI E-utilities"""
        print(f"Searching PubMed for '{query}'...")
        
        # Create a properly formatted search query with date range
        search_query = f"{query} AND {self.date_from}:{self.date_to}[Date - Publication]"
        
        # Build the API parameters
        params = {
            "db": "pubmed",
            "term": search_query,
            "retmax": self.config["result_limit"],
            "usehistory": "y",
            "retmode": "json"
        }
        
        # Add API key if available
        if self.config["ncbi_api_key"]:
            params["api_key"] = self.config["ncbi_api_key"]
        
        # Step 1: Use ESearch to get IDs
        max_retries = 3
        for attempt in range(max_retries):
            try:
                response = requests.get(self.base_esearch, params=params, timeout=30)
                response.raise_for_status()
                
                search_data = response.json()
                
                # Check if we have results
                if 'esearchresult' in search_data and 'idlist' in search_data['esearchresult']:
                    id_list = search_data['esearchresult']['idlist']
                    
                    if not id_list:
                        print(f"No results found for '{query}'")
                        return pd.DataFrame()
                    
                    print(f"Found {len(id_list)} papers for '{query}'")
                    
                    # Get WebEnv and QueryKey from the search for efficient fetching
                    web_env = search_data['esearchresult']['webenv']
                    query_key = search_data['esearchresult']['querykey']
                    
                    # Step 2: Use EFetch to get full records
                    return self.fetch_paper_details(id_list, web_env, query_key, query)
                    
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
    
    def fetch_paper_details(self, id_list, web_env, query_key, original_query):
        """Fetch detailed paper information using EFetch"""
        # Prepare parameters for EFetch
        params = {
            "db": "pubmed",
            "WebEnv": web_env,
            "query_key": query_key,
            "retmode": "xml",
            "rettype": "abstract"
        }
        
        # Add API key if available
        if self.config["ncbi_api_key"]:
            params["api_key"] = self.config["ncbi_api_key"]
        
        # For larger result sets, fetch in batches to avoid timeouts
        batch_size = 50
        all_results = []
        
        for start in range(0, len(id_list), batch_size):
            end = min(start + batch_size, len(id_list))
            batch_ids = id_list[start:end]
            
            # Update parameters for this batch
            params["id"] = ",".join(batch_ids)
            
            max_retries = 3
            for attempt in range(max_retries):
                try:
                    response = requests.get(self.base_efetch, params=params, timeout=30)
                    response.raise_for_status()
                    
                    # Parse XML response
                    papers = self.parse_pubmed_xml(response.text, original_query)
                    all_results.extend(papers)
                    
                    # Be nice to the API
                    time.sleep(0.5)
                    break
                    
                except requests.exceptions.RequestException as e:
                    if attempt < max_retries - 1:
                        wait_time = 2 ** attempt
                        print(f"API request failed, retrying in {wait_time} seconds...")
                        time.sleep(wait_time)
                    else:
                        print(f"Error fetching paper details: {e}")
        
        # Convert to DataFrame
        if all_results:
            return pd.DataFrame(all_results)
        else:
            return pd.DataFrame()
    
    def parse_pubmed_xml(self, xml_text, original_query):
        """Parse PubMed XML response to extract paper details"""
        papers = []
        
        try:
            # Parse XML
            root = ET.fromstring(xml_text)
            
            # Find all PubmedArticle elements
            for article in root.findall(".//PubmedArticle"):
                try:
                    # Extract PMID
                    pmid = article.find(".//PMID").text if article.find(".//PMID") is not None else ""
                    
                    # Extract article title
                    title_element = article.find(".//ArticleTitle")
                    title = title_element.text if title_element is not None else "No title"
                    
                    # Extract abstract
                    abstract_elements = article.findall(".//AbstractText")
                    abstract = " ".join([elem.text for elem in abstract_elements if elem.text]) if abstract_elements else ""
                    
                    # Extract authors
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
                    
                    # Extract publication date
                    pub_date = article.find(".//PubDate")
                    year = pub_date.find("Year").text if pub_date.find("Year") is not None else ""
                    month = pub_date.find("Month").text if pub_date.find("Month") is not None else ""
                    day = pub_date.find("Day").text if pub_date.find("Day") is not None else ""
                    
                    date = f"{year}-{month}-{day}" if day else (f"{year}-{month}" if month else year)
                    
                    # Extract journal information
                    journal_element = article.find(".//Journal/Title")
                    journal = journal_element.text if journal_element is not None else ""
                    
                    # Extract DOI if available
                    doi_element = article.find(".//ArticleId[@IdType='doi']")
                    doi = doi_element.text if doi_element is not None else ""
                    
                    # Create a URL to the PubMed entry
                    url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else ""
                    
                    # Create paper dictionary
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
    
    def run_searches(self):
        """Run all configured searches and combine results"""
        print(f"Searching PubMed for papers from {self.date_from} to {self.date_to}")
        
        # Individual term searches
        term_dfs = []
        
        # Use ThreadPoolExecutor for parallel searches
        with ThreadPoolExecutor(max_workers=3) as executor:
            # Submit individual term searches
            term_futures = {executor.submit(self.search_pubmed, term): term 
                           for term in self.config["search_terms"]}
            
            # Submit combined searches
            combined_futures = {executor.submit(self.search_pubmed, query): query 
                               for query in self.config["combined_searches"]}
            
            # Collect results from individual terms
            for future in term_futures:
                term = term_futures[future]
                try:
                    df = future.result()
                    if not df.empty:
                        print(f"Processed {len(df)} papers for term '{term}'")
                        term_dfs.append(df)
                except Exception as e:
                    print(f"Error processing results for '{term}': {e}")
            
            # Collect results from combined searches
            for future in combined_futures:
                query = combined_futures[future]
                try:
                    df = future.result()
                    if not df.empty:
                        print(f"Processed {len(df)} papers for query '{query}'")
                        term_dfs.append(df)
                except Exception as e:
                    print(f"Error processing results for '{query}': {e}")
        
        # Combine all results
        if term_dfs:
            self.all_papers = pd.concat(term_dfs)
            
            # Group by PMID to consolidate search terms
            if 'pmid' in self.all_papers.columns and not self.all_papers['pmid'].isna().all():
                # First, make sure PMID is not empty
                valid_pmids = self.all_papers['pmid'].notna() & (self.all_papers['pmid'] != '')
                
                # Split dataframe into records with and without PMIDs
                papers_with_pmid = self.all_papers[valid_pmids]
                papers_without_pmid = self.all_papers[~valid_pmids]
                
                # Group by PMID and aggregate search terms
                if not papers_with_pmid.empty:
                    papers_with_pmid = papers_with_pmid.groupby('pmid', as_index=False).agg({
                        'title': 'first',
                        'authors': 'first',
                        'date': 'first',
                        'doi': 'first',
                        'abstract': 'first',
                        'journal': 'first',
                        'url': 'first',
                        'search_term': lambda x: ', '.join(sorted(set(x)))  # Combine all search terms
                    })
                
                # Group papers without PMID by DOI or title
                if not papers_without_pmid.empty:
                    if 'doi' in papers_without_pmid.columns and not papers_without_pmid['doi'].isna().all():
                        # Try grouping by DOI first
                        papers_without_pmid = papers_without_pmid.groupby('doi', as_index=False).agg({
                            'title': 'first',
                            'authors': 'first',
                            'date': 'first',
                            'pmid': 'first',
                            'abstract': 'first',
                            'journal': 'first',
                            'url': 'first',
                            'search_term': lambda x: ', '.join(sorted(set(x)))
                        })
                    else:
                        # If DOI not available, use title
                        papers_without_pmid = papers_without_pmid.groupby('title', as_index=False).agg({
                            'authors': 'first',
                            'date': 'first',
                            'pmid': 'first',
                            'doi': 'first',
                            'abstract': 'first',
                            'journal': 'first',
                            'url': 'first',
                            'search_term': lambda x: ', '.join(sorted(set(x)))
                        })
                
                # Combine the results back
                self.all_papers = pd.concat([papers_with_pmid, papers_without_pmid])
                
            else:
                # If PMID is not available, try DOI, then fall back to title
                if 'doi' in self.all_papers.columns and not self.all_papers['doi'].isna().all():
                    self.all_papers = self.all_papers.groupby('doi', as_index=False).agg({
                        'title': 'first',
                        'authors': 'first',
                        'date': 'first',
                        'pmid': 'first',
                        'abstract': 'first',
                        'journal': 'first',
                        'url': 'first',
                        'search_term': lambda x: ', '.join(sorted(set(x)))
                    })
                else:
                    self.all_papers = self.all_papers.groupby('title', as_index=False).agg({
                        'authors': 'first',
                        'date': 'first',
                        'pmid': 'first',
                        'doi': 'first',
                        'abstract': 'first',
                        'journal': 'first',
                        'url': 'first',
                        'search_term': lambda x: ', '.join(sorted(set(x)))
                    })
            
            print(f"Combined unique papers: {len(self.all_papers)}")
        else:
            self.all_papers = pd.DataFrame()
            print("No papers found matching your search criteria.")
    
    def score_papers(self):
        """Score papers based on relevance criteria"""
        if self.all_papers.empty:
            print("No papers to score.")
            return
        
        # Add a score column
        self.all_papers['score'] = 0
        
        # Make sure all string columns are actually strings to avoid errors
        string_columns = ['title', 'authors', 'abstract', 'journal']
        for col in string_columns:
            if col in self.all_papers.columns:
                self.all_papers[col] = self.all_papers[col].fillna('').astype(str)
            else:
                self.all_papers[col] = ''
                print(f"Warning: Column '{col}' not found, using empty strings instead")
        
        # Score based on title relevance to specific systems (most important: +3)
        system_pattern = '|'.join(self.config["relevant_systems"])
        self.all_papers['system_match'] = self.all_papers['title'].str.lower().str.contains(
            system_pattern, regex=True, case=False, na=False).astype(int) * 3
        # Also check abstract for system match
        self.all_papers['system_match_abstract'] = self.all_papers['abstract'].str.lower().str.contains(
            system_pattern, regex=True, case=False, na=False).astype(int) * 2
        # Use the higher of the two scores
        self.all_papers['system_match'] = self.all_papers[['system_match', 'system_match_abstract']].max(axis=1)
        self.all_papers['score'] += self.all_papers['system_match']
        
        # Score for computational methods (important: +2)
        comp_pattern = '|'.join(self.config["computational_keywords"])
        self.all_papers['computational_match'] = self.all_papers['title'].str.lower().str.contains(
            comp_pattern, regex=True, case=False, na=False).astype(int) * 2
        # Also check abstract for computational keywords
        self.all_papers['computational_match_abstract'] = self.all_papers['abstract'].str.lower().str.contains(
            comp_pattern, regex=True, case=False, na=False).astype(int) * 1
        # Use the higher of the two scores
        self.all_papers['computational_match'] = self.all_papers[['computational_match', 'computational_match_abstract']].max(axis=1)
        self.all_papers['score'] += self.all_papers['computational_match']
        
        # Score for important journals (useful: +2)
        journal_pattern = '|'.join(self.config["important_journals"])
        self.all_papers['journal_match'] = self.all_papers['journal'].str.lower().str.contains(
            journal_pattern, regex=True, case=False, na=False).astype(int) * 2
        self.all_papers['score'] += self.all_papers['journal_match']
        
        # Look for combined aspects (experimental + computational) (important: +2)
        # Only if abstract is available
        if not self.all_papers['abstract'].str.len().eq(0).all():
            self.all_papers['combined_approach'] = (
                self.all_papers['abstract'].str.lower().str.contains('experimental', regex=False, na=False) &
                self.all_papers['abstract'].str.lower().str.contains('computational', regex=False, na=False)
            ).astype(int) * 2
            self.all_papers['score'] += self.all_papers['combined_approach']
        
        # Check for important authors (useful but less critical: +1)
        author_pattern = '|'.join(self.config["important_authors"])
        self.all_papers['author_match'] = self.all_papers['authors'].str.contains(
            author_pattern, regex=True, case=False, na=False).astype(int)
        self.all_papers['score'] += self.all_papers['author_match']
        
        # Look for potential datasets (useful: +1)
        # Only if abstract is available
        if not self.all_papers['abstract'].str.len().eq(0).all():
            self.all_papers['dataset_available'] = self.all_papers['abstract'].str.lower().str.contains(
                'dataset|data available|github|available at|repository', regex=True, na=False).astype(int)
            self.all_papers['score'] += self.all_papers['dataset_available']
        
        # Penalize reviews slightly unless comprehensive (-1)
        self.all_papers['is_review'] = self.all_papers['title'].str.lower().str.contains(
            'review|overview|survey', regex=True, na=False).astype(int)
        self.all_papers['comprehensive_review'] = self.all_papers['title'].str.lower().str.contains(
            'comprehensive|systematic', regex=True, na=False).astype(int)
        self.all_papers['review_score'] = -1 * (self.all_papers['is_review'] & ~self.all_papers['comprehensive_review'])
        self.all_papers['score'] += self.all_papers['review_score']
        
        # Sort by score (descending)
        self.all_papers = self.all_papers.sort_values('score', ascending=False)
    
    def format_authors(self, authors_string):
        """Format authors to show only first and last three with ... in between"""
        # Split into individual authors
        authors = [a.strip() for a in authors_string.split(',') if a.strip()]
        
        if len(authors) <= 4:
            return ", ".join(authors)
        else:
            return f"{authors[0]}, ... {', '.join(authors[-3:])}"
    
    def get_top_papers(self, limit=20):
        """Get the top N papers by score"""
        if self.all_papers.empty:
            return pd.DataFrame()
        
        # Filter by minimum score
        filtered = self.all_papers[self.all_papers['score'] >= self.config["min_score"]]
        
        # Return top N papers
        return filtered.head(limit)
    
    def print_results(self, limit=20):
        """Print the top results in a formatted way"""
        top_papers = self.get_top_papers(limit)
        
        if top_papers.empty:
            print(f"No papers met your minimum score of {self.config['min_score']}.")
            return
        
        print(f"\n=== Top {len(top_papers)} Papers (minimum score: {self.config['min_score']}) ===\n")
        
        for i, (_, paper) in enumerate(top_papers.iterrows(), 1):
            print(f"{i}. {paper['title']}")
            
            # Format authors to show only first and last three
            formatted_authors = self.format_authors(paper['authors'])
            print(f"   Authors: {formatted_authors}")
            
            # Print journal and date
            print(f"   Journal: {paper['journal']} (Date: {paper['date']})")
            
            # Print PMID and DOI if available
            if 'pmid' in paper and paper['pmid']:
                print(f"   PMID: {paper['pmid']}")
            if 'doi' in paper and paper['doi']:
                print(f"   DOI: {paper['doi']}")
            if 'url' in paper and paper['url']:
                print(f"   URL: {paper['url']}")
                
            # Print search terms that found this paper
            if 'search_term' in paper and paper['search_term']:
                print(f"   Found via: {paper['search_term']}")
                
            print(f"   Score: {paper['score']}")
            
            # Print beginning of abstract if available
            if 'abstract' in paper and paper['abstract'] and len(paper['abstract']) > 10:
                abstract_preview = paper['abstract'][:150] + "..." if len(paper['abstract']) > 150 else paper['abstract']
                print(f"   Abstract: {abstract_preview}")
                
            print()
    
    def save_results(self, filename="pubmed_papers.csv"):
        """Save results to a CSV file with all relevant information"""
        if self.all_papers.empty:
            print("No results to save.")
            return
            
        # Create formatted authors column
        self.all_papers['formatted_authors'] = self.all_papers['authors'].apply(self.format_authors)
        
        # Add matching_terms column that lists terms found in title and abstract
        self.all_papers['matching_terms'] = ''
        
        # Check which system terms match
        for system in self.config["relevant_systems"]:
            mask = (
                self.all_papers['title'].str.lower().str.contains(system.lower(), regex=False, na=False) | 
                self.all_papers['abstract'].str.lower().str.contains(system.lower(), regex=False, na=False)
            )
            self.all_papers.loc[mask, 'matching_terms'] += system + ', '
        
        # Check which computational terms match
        for term in self.config["computational_keywords"]:
            mask = (
                self.all_papers['title'].str.lower().str.contains(term.lower(), regex=False, na=False) | 
                self.all_papers['abstract'].str.lower().str.contains(term.lower(), regex=False, na=False)
            )
            self.all_papers.loc[mask, 'matching_terms'] += term + ', '
        
        # Clean up the matching_terms column
        self.all_papers['matching_terms'] = self.all_papers['matching_terms'].str.rstrip(', ')
        
        # Make sure the search_term column exists (in case it wasn't populated)
        if 'search_term' not in self.all_papers.columns:
            self.all_papers['search_term'] = ''
        
        # Define columns to include in output
        columns_to_save = [
            'title', 
            'authors', 
            'formatted_authors',
            'journal',
            'date', 
            'pmid',
            'doi', 
            'score', 
            'search_term', 
            'matching_terms',
            'abstract',
            'url'
        ]
        
        # Only include columns that exist in the dataframe
        available_columns = [col for col in columns_to_save if col in self.all_papers.columns]
        
        # Save to CSV
        self.all_papers[available_columns].to_csv(filename, index=False)
        print(f"Results saved to {filename}")
        
        # Also save a simplified version with just the essential columns
        simple_filename = filename.replace('.csv', '_simple.csv')
        if simple_filename == filename:  # If no .csv extension was found
            simple_filename = filename + '_simple'
            
        simple_columns = [col for col in ['title', 'formatted_authors', 'journal', 'date', 'pmid', 'doi', 'score', 'search_term'] 
                          if col in self.all_papers.columns]
        
        self.all_papers[simple_columns].to_csv(simple_filename, index=False)
        print(f"Simplified results saved to {simple_filename}")
    
    def save_config(self, filename="pubmed_config.json"):
        """Save current configuration to a JSON file"""
        with open(filename, 'w') as f:
            json.dump(self.config, f, indent=2)
        print(f"Configuration saved to {filename}")


if __name__ == "__main__":
    # Set up command line arguments
    parser = argparse.ArgumentParser(description='Search and prioritize PubMed papers')
    parser.add_argument('--config', help='Path to configuration file')
    parser.add_argument('--days', type=int, help='Number of days to look back')
    parser.add_argument('--limit', type=int, default=20, help='Maximum papers to return')
    parser.add_argument('--min-score', type=int, help='Minimum score threshold')
    parser.add_argument('--save-config', action='store_true', help='Save current config to JSON')
    parser.add_argument('--output', help='Output CSV filename')
    parser.add_argument('--api-key', help='NCBI API key for higher rate limits')
    
    args = parser.parse_args()
    
    # Create scraper with optional config file
    scraper = PubMedScraper(config_file=args.config)
    
    # Override days_back if specified
    if args.days:
        scraper.config["days_back"] = args.days
    
    # Override min_score if specified
    if args.min_score:
        scraper.config["min_score"] = args.min_score
    
    # Override API key if specified
    if args.api_key:
        scraper.config["ncbi_api_key"] = args.api_key
    
    # Run the search and scoring
    scraper.run_searches()
    scraper.score_papers()
    scraper.print_results(limit=args.limit)
    
    # Save results if requested
    if args.output:
        scraper.save_results(args.output)
    
    # Save config if requested
    if args.save_config:
        scraper.save_config()