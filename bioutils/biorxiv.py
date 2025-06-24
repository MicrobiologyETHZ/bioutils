#!/usr/bin/env python3
"""
BioRxiv Paper Scraper and Prioritizer

This script searches for and prioritizes preprints from bioRxiv based on user-specified
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
from bs4 import BeautifulSoup  # Make sure to install this: pip install beautifulsoup4

class BioRxivScraper:
    """
    A class to search, scrape, and prioritize papers from bioRxiv based on user-defined criteria.
    """
    
    def __init__(self, config_file=None):
        """Initialize with default config or from a file"""
        # Check if BeautifulSoup is installed
        try:
            import bs4
        except ImportError:
            print("ERROR: BeautifulSoup4 is required for this script to work.")
            print("Please install it with: pip install beautifulsoup4")
            sys.exit(1)  # Exit with error code
            
        # Default configuration
        self.config = {
            "search_terms": [
                "microbiome", "microbial+community",
                "metatranscriptomics", 
                "bacterial+growth", "bacterial+metabolism", 
            ],
            "combined_searches": [
                "microbiome+transcriptomics",
                "microbiota+growth",
                "bacterial+transcriptomics+computational"
            ],
            "important_authors": [
                # Add authors you follow closely
                "Knight, Rob", "Segata, Nicola", "Huttenhower, Curtis",
                "Gilbert, Jack", "Turnbaugh, Peter", "Fischbach, Michael"
            ],
            "relevant_systems": [
                # Your specific systems of interest
                "soil+microbiome", 
                "marine+microbiome", "plant+microbiome", "oral+microbiome",
                "ocean+microbiome", 'coral+microbiome', "gut+microbiome"
            ],
            "computational_keywords": [
                "machine learning", "deep learning", "network analysis",
                "metabolic modeling", "computational", "niche partitioning", 
                "labour divisiont", "life strategies"
            ],
            "days_back": 7,
            "min_score": 3
        }
        
        # Load custom configuration if provided
        if config_file:
            with open(config_file, 'r') as f:
                custom_config = json.load(f)
                self.config.update(custom_config)
        
        # Calculate date range for searches
        self.today = datetime.now()
        self.start_date = self.today - timedelta(days=self.config["days_back"])
        
        # Prepare results storage
        self.all_papers = pd.DataFrame()
    
    def search_biorxiv(self, query):
        """Search bioRxiv preprints for a specific query within date range using direct scraping"""
        date_str = self.start_date.strftime("%Y-%m-%d")
        today_str = self.today.strftime("%Y-%m-%d")
        
        # Use modern browser headers
        headers = {
            "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36",
            "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,*/*;q=0.8",
            "Accept-Language": "en-US,en;q=0.5",
            "Connection": "keep-alive",
            "Upgrade-Insecure-Requests": "1",
            "Sec-Fetch-Dest": "document",
            "Sec-Fetch-Mode": "navigate",
            "Sec-Fetch-Site": "none",
            "Sec-Fetch-User": "?1",
            "DNT": "1"
        }
        
        # Direct HTML scraping from the bioRxiv website
        print(f"Scraping bioRxiv for '{query}'...")
        
        max_retries = 3
        for attempt in range(max_retries):
            try:
                # Properly encode the search URL
                # Format: https://www.biorxiv.org/search/microbiome%20jcode%3Abiorxiv%20limit_from%3A2025-03-01%20limit_to%3A2025-03-10
                search_url = f"https://www.biorxiv.org/search/{query}%20jcode%3Abiorxiv%20limit_from%3A{date_str}%20limit_to%3A{today_str}%20numresults%3A100%20sort%3Apublication-date"
                # Add a slight delay to avoid hammering the server
                time.sleep(1)
                
                response = requests.get(search_url, headers=headers, timeout=30)
                
                if response.status_code == 200:
                    soup = BeautifulSoup(response.text, 'html.parser')
                    
                    # Check if we have results or a "no results" message
                    no_results = soup.select('.highwire-search-no-results')
                    if no_results:
                        print(f"No results found for '{query}'")
                        return pd.DataFrame()
                    
                    # Find all paper entries
                    paper_entries = soup.select('.highwire-article-citation')
                    
                    if not paper_entries:
                        print(f"Warning: Found response page but no paper entries for '{query}'")
                        # Save the HTML for debugging if needed
                        with open(f"debug_{query.replace('+', '_')}.html", 'w', encoding='utf-8') as f:
                            f.write(response.text)
                        return pd.DataFrame()
                    
                    print(f"Found {len(paper_entries)} papers for '{query}'")
                    
                    results = []
                    for entry in paper_entries:
                        try:
                            # Extract key information
                            title_elem = entry.select_one('.highwire-cite-title')
                            authors_elem = entry.select_one('.highwire-citation-authors')
                            date_elem = entry.select_one('.highwire-cite-metadata-date')
                            doi_elem = entry.select_one('.highwire-cite-metadata-doi')
                            
                            # Get abstract if available (sometimes shown in search results)
                            abstract_elem = entry.select_one('.highwire-cite-snippet')
                            
                            # Get paper URL for fetching abstract later if needed
                            paper_url = None
                            if title_elem and title_elem.find('a'):
                                paper_url = title_elem.find('a').get('href')
                                if paper_url and not paper_url.startswith('http'):
                                    paper_url = 'https://www.biorxiv.org' + paper_url
                            
                            if title_elem and authors_elem:
                                paper = {
                                    'title': title_elem.get_text(strip=True),
                                    'authors': authors_elem.get_text(strip=True).replace('\n', ', '),
                                    'date': date_elem.get_text(strip=True) if date_elem else '',
                                    'doi': doi_elem.get_text(strip=True).replace('doi:', '').strip() if doi_elem else '',
                                    'abstract': abstract_elem.get_text(strip=True) if abstract_elem else '',
                                    'version': '1',  # Default version
                                    'url': paper_url,
                                    'search_term': query  # Add the search term that found this paper
                                }
                                results.append(paper)
                        except Exception as e:
                            print(f"Error parsing paper entry: {e}")
                            continue
                    
                    # Create DataFrame
                    df = pd.DataFrame(results)
                    return df
                else:
                    print(f"HTTP Error: {response.status_code} for '{query}'")
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
    
    def run_searches(self):
        """Run all configured searches and combine results"""
        print(f"Searching bioRxiv for papers from {self.start_date.strftime('%Y-%m-%d')} to {self.today.strftime('%Y-%m-%d')}")
        
        # Individual term searches
        term_dfs = []
        
        # Use ThreadPoolExecutor for parallel searches
        with ThreadPoolExecutor(max_workers=3) as executor:  # Reduced from 5 to 3 to avoid overloading server
            # Submit individual term searches
            term_futures = {executor.submit(self.search_biorxiv, term): term 
                           for term in self.config["search_terms"]}
            
            # Submit combined searches
            combined_futures = {executor.submit(self.search_biorxiv, query): query 
                               for query in self.config["combined_searches"]}
            
            # Collect results from individual terms
            for future in term_futures:
                term = term_futures[future]
                try:
                    df = future.result()
                    if not df.empty:
                        print(f"Found {len(df)} papers for term '{term}'")
                        term_dfs.append(df)
                except Exception as e:
                    print(f"Error processing results for '{term}': {e}")
            
            # Collect results from combined searches
            for future in combined_futures:
                query = combined_futures[future]
                try:
                    df = future.result()
                    if not df.empty:
                        print(f"Found {len(df)} papers for query '{query}'")
                        term_dfs.append(df)
                except Exception as e:
                    print(f"Error processing results for '{query}': {e}")
        
        # Combine all results
        if term_dfs:
            self.all_papers = pd.concat(term_dfs)
            
            # Group by DOI to consolidate search terms
            if 'doi' in self.all_papers.columns and not self.all_papers['doi'].isna().all():
                # First, make sure DOI is not empty
                valid_dois = self.all_papers['doi'].notna() & (self.all_papers['doi'] != '')
                
                # Split dataframe into records with and without DOIs
                papers_with_doi = self.all_papers[valid_dois]
                papers_without_doi = self.all_papers[~valid_dois]
                
                # Group by DOI and aggregate search terms
                if not papers_with_doi.empty:
                    papers_with_doi = papers_with_doi.groupby('doi', as_index=False).agg({
                        'title': 'first',
                        'authors': 'first',
                        'date': 'first',
                        'abstract': 'first',
                        'version': 'first',
                        'url': 'first',
                        'search_term': lambda x: ', '.join(sorted(set(x)))  # Combine all search terms
                    })
                
                # Group papers without DOI by title
                if not papers_without_doi.empty:
                    papers_without_doi = papers_without_doi.groupby('title', as_index=False).agg({
                        'authors': 'first',
                        'date': 'first',
                        'doi': 'first',
                        'abstract': 'first',
                        'version': 'first',
                        'url': 'first',
                        'search_term': lambda x: ', '.join(sorted(set(x)))  # Combine all search terms
                    })
                
                # Combine the results back
                self.all_papers = pd.concat([papers_with_doi, papers_without_doi])
                
            else:
                # If DOI is not available, use title for deduplication
                self.all_papers = self.all_papers.groupby('title', as_index=False).agg({
                    'authors': 'first',
                    'date': 'first',
                    'doi': 'first',
                    'abstract': 'first', 
                    'version': 'first',
                    'url': 'first',
                    'search_term': lambda x: ', '.join(sorted(set(x)))  # Combine all search terms
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
        string_columns = ['title', 'authors', 'abstract']
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
            comp_pattern, regex=True, case=False, na=False).astype(int) * 2
        # Use the higher of the two scores
        self.all_papers['computational_match'] = self.all_papers[['computational_match', 'computational_match_abstract']].max(axis=1)
        self.all_papers['score'] += self.all_papers['computational_match']
        
        
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
        
        # Add length of search terms to the score
        self.all_papers['search_term_score'] = self.all_papers['search_term'].apply(lambda x: len(x.split(","))-1)
        self.all_papers['score'] += self.all_papers['search_term_score']

        # Sort by score (descending)
        self.all_papers = self.all_papers.sort_values('score', ascending=False)
    
    def format_authors(self, authors_string):
        """Format authors to show only first and last three with ... in between"""
        # Remove "View ORCID Profile" text which is often present in bioRxiv author strings
        cleaned = authors_string.replace("View ORCID Profile", "").replace("  ", " ").strip()
        
        # Split into individual authors
        authors = [a.strip() for a in re.split(r',|\sand\s', cleaned) if a.strip()]
        
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
            
            # Print DOI if available, otherwise URL
            if paper['doi'] and paper['doi'].strip():
                print(f"   DOI: {paper['doi']} (Date: {paper['date']})")
            elif 'url' in paper and paper['url']:
                print(f"   URL: {paper['url']} (Date: {paper['date']})")
            else:
                print(f"   Date: {paper['date']}")
                
            # Print search terms that found this paper
            if 'search_term' in paper and paper['search_term']:
                print(f"   Found via: {paper['search_term']}")
                
            print(f"   Score: {paper['score']}")
            
            # Print beginning of abstract if available
            if 'abstract' in paper and paper['abstract'] and len(paper['abstract']) > 10:
                abstract_preview = paper['abstract'][:150] + "..." if len(paper['abstract']) > 150 else paper['abstract']
                print(f"   Abstract: {abstract_preview}")
                
            print()
    
    def save_results(self, filename="biorxiv_papers.csv"):
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
            #'authors', 
            'formatted_authors',
            'doi', 
            'date', 
            'score', 
            'search_term', 
            'matching_terms',
            #'abstract',
            #'url'
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
            
        simple_columns = [col for col in ['title', 'formatted_authors', 'doi', 'date', 'score', 'search_term'] 
                          if col in self.all_papers.columns]
        
        self.all_papers[simple_columns].to_csv(simple_filename, index=False)
        print(f"Simplified results saved to {simple_filename}")
    
    def save_config(self, filename="biorxiv_config.json"):
        """Save current configuration to a JSON file"""
        with open(filename, 'w') as f:
            json.dump(self.config, f, indent=2)
        print(f"Configuration saved to {filename}")


if __name__ == "__main__":
    # Set up command line arguments
    parser = argparse.ArgumentParser(description='Search and prioritize bioRxiv papers')
    parser.add_argument('--config', help='Path to configuration file')
    parser.add_argument('--days', type=int, help='Number of days to look back')
    parser.add_argument('--limit', type=int, default=20, help='Maximum papers to return')
    parser.add_argument('--min-score', type=int, help='Minimum score threshold')
    parser.add_argument('--save-config', action='store_true', help='Save current config to JSON')
    parser.add_argument('--output', help='Output CSV filename')
    
    args = parser.parse_args()
    
    # Create scraper with optional config file
    scraper = BioRxivScraper(config_file=args.config)
    
    # Override days_back if specified
    if args.days:
        scraper.config["days_back"] = args.days
    
    # Override min_score if specified
    if args.min_score:
        scraper.config["min_score"] = args.min_score
    
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