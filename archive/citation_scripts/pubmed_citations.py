#!/usr/bin/env python3
"""
PubMed Citations Collector for Awesome-Virome

This script collects academic citation data from PubMed for tools in the Awesome-Virome 
repository. It uses the NCBI E-utilities API to search for relevant publications, retrieve
metadata, and store formatted citations.

Usage:
    python pubmed_citations.py [--output OUTPUT] [--api-key API_KEY] [--email EMAIL]
    
Environment variables:
    NCBI_API_KEY: NCBI API key for higher rate limits
    CONTACT_EMAIL: Contact email for API rate limiting
"""

import os
import json
import time
import logging
import argparse
from datetime import datetime
from typing import Dict, List, Optional, Any, Tuple
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import re
import sys

# Add scripts directory to path to find the APIs module
scripts_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(scripts_dir)

from apis.citations_api import (
    PubMedAPI,
    CrossRefAPI,
    GitHubAPI,
    cache_manager
)

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("pubmed_citations.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Directory constants
METADATA_DIR = os.path.join("metadata", "pubmed_citations")
OUTPUT_FILE = os.path.join(METADATA_DIR, "pubmed_citations.json")
SUMMARY_FILE = os.path.join(METADATA_DIR, "summary.json")
CITATIONS_FILE = "citations.md"


class PubMedCitationsCollector:
    """Collects and manages PubMed citation data for tools."""
    
    def __init__(self, api_key: Optional[str] = None, email: Optional[str] = None,
                 github_token: Optional[str] = None, metadata_dir: str = METADATA_DIR):
        """Initialize the PubMed citations collector."""
        self.pubmed_api = PubMedAPI(api_key=api_key, email=email)
        self.crossref_api = CrossRefAPI(email=email)
        self.github_api = GitHubAPI(token=github_token)
        
        # Store directory paths
        self.metadata_dir = metadata_dir
        self.output_file = os.path.join(self.metadata_dir, "pubmed_citations.json")
        self.summary_file = os.path.join(self.metadata_dir, "summary.json")
        
        # Create directories
        os.makedirs(self.metadata_dir, exist_ok=True)
        
        # Load existing metadata if available
        self.existing_metadata = {}
        if os.path.exists(self.output_file):
            try:
                with open(self.output_file, 'r') as f:
                    self.existing_metadata = json.load(f)
            except json.JSONDecodeError:
                logger.warning(f"Could not parse existing metadata from {self.output_file}")
        
        # Statistics for progress tracking
        self.stats = {
            "tools_processed": 0,
            "tools_with_publications": 0,
            "total_publications": 0,
            "start_time": datetime.now().isoformat()
        }
    
    def get_tool_publication(self, tool: Dict[str, Any]) -> Dict[str, Any]:
        """Find the best publication for a tool using PubMed."""
        tool_name = tool.get('name', '')
        repo_url = tool.get('url', '')
        description = tool.get('description', '')
        
        logger.info(f"Finding publications for tool: {tool_name}")
        
        # Check if we already have recent data for this tool
        cache_key = f"{tool_name}_{repo_url}"
        if cache_key in self.existing_metadata:
            last_updated = self.existing_metadata[cache_key].get('last_updated', '')
            # If updated within the last 30 days, use cached data
            if last_updated and (datetime.now() - datetime.fromisoformat(last_updated)).days < 30:
                logger.info(f"Using cached citation data for {tool_name}")
                return self.existing_metadata[cache_key]
        
        # First try from the GitHub repo if available
        citation_info = {}
        if repo_url and 'github.com' in repo_url:
            citation_file = self.github_api.find_citation_file(repo_url)
            if citation_file:
                logger.info(f"Found citation file for {tool_name}")
                if citation_file['format'] == 'cff':
                    citation_data = self.github_api.parse_citation_cff(citation_file['content'])
                    if citation_data.get('doi'):
                        # Get citation information from DOI
                        publication = self.crossref_api.get_work_by_doi(citation_data['doi'])
                        if publication:
                            citation_info = {
                                'source': 'github_citation_file',
                                'publication': publication,
                                'formatted_citations': {
                                    'apa': self.crossref_api.format_citation(publication, 'apa'),
                                    'bibtex': self.crossref_api.format_citation(publication, 'bibtex')
                                }
                            }
                else:
                    # Extract DOI from other citation formats
                    doi = self.github_api.extract_doi_from_text(citation_file['content'])
                    if doi:
                        publication = self.crossref_api.get_work_by_doi(doi)
                        if publication:
                            citation_info = {
                                'source': 'github_citation_file',
                                'publication': publication,
                                'formatted_citations': {
                                    'apa': self.crossref_api.format_citation(publication, 'apa'),
                                    'bibtex': self.crossref_api.format_citation(publication, 'bibtex')
                                }
                            }
                        
        # If no citation file found, search PubMed
        if not citation_info:
            publication = self.pubmed_api.find_best_publication_for_tool(tool_name, repo_url)
            if publication:
                # Format citations
                citation_info = {
                    'source': 'pubmed',
                    'publication': publication,
                    'formatted_citations': {
                        'apa': self.pubmed_api.format_citation(publication, 'apa'),
                        'bibtex': self.pubmed_api.format_citation(publication, 'bibtex'),
                        'mla': self.pubmed_api.format_citation(publication, 'mla')
                    }
                }
        
        # Compile result with tool information
        result = {
            'name': tool_name,
            'url': repo_url,
            'citation_info': citation_info,
            'last_updated': datetime.now().isoformat()
        }
        
        # Save individual tool data
        tool_filename = re.sub(r'[^\w\-\.]', '_', tool_name)
        tool_file = os.path.join(self.metadata_dir, f"{tool_filename}.json")
        with open(tool_file, 'w') as f:
            json.dump(result, f, indent=2)
        
        return result
    
    def collect_pubmed_citations(self, tools: List[Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
        """Collect PubMed citation data for all tools."""
        results = {}
        
        # Process tools in parallel with a limited number of workers to avoid overloading the API
        with ThreadPoolExecutor(max_workers=5) as executor:
            future_to_tool = {executor.submit(self.get_tool_publication, tool): tool for tool in tools}
            
            for future in as_completed(future_to_tool):
                tool = future_to_tool[future]
                tool_name = tool.get('name', '')
                repo_url = tool.get('url', '')
                
                try:
                    result = future.result()
                    cache_key = f"{tool_name}_{repo_url}"
                    results[cache_key] = result
                    
                    # Update statistics
                    self.stats["tools_processed"] += 1
                    if result.get('citation_info'):
                        self.stats["tools_with_publications"] += 1
                        self.stats["total_publications"] += 1
                    
                    # Log progress periodically
                    if self.stats["tools_processed"] % 10 == 0:
                        logger.info(f"Processed {self.stats['tools_processed']} tools. " +
                                   f"Found publications for {self.stats['tools_with_publications']} tools.")
                
                except Exception as e:
                    logger.error(f"Error processing tool {tool_name}: {e}")
        
        # Save all results
        with open(self.output_file, 'w') as f:
            json.dump(results, f, indent=2)
        
        # Generate summary
        self.generate_summary(results)
        
        return results
    
    def generate_summary(self, results: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
        """Generate a summary of the citation data."""
        total_tools = len(results)
        tools_with_publications = sum(1 for data in results.values() if data.get('citation_info'))
        
        # Find most cited publications
        citation_counts = []
        for tool_key, data in results.items():
            citation_info = data.get('citation_info', {})
            publication = citation_info.get('publication', {})
            
            # Different sources might store citation count differently
            citation_count = 0
            if 'citation_count' in publication:
                citation_count = publication.get('citation_count', 0)
            
            if citation_count > 0:
                citation_counts.append((
                    data.get('name', ''),
                    data.get('url', ''),
                    citation_count,
                    publication.get('title', ''),
                    publication.get('journal', ''),
                    publication.get('year', '')
                ))
        
        # Sort by citation count (descending)
        most_cited = sorted(citation_counts, key=lambda x: x[2], reverse=True)[:20]
        
        # Count publications by year
        publications_by_year = {}
        for data in results.values():
            citation_info = data.get('citation_info', {})
            publication = citation_info.get('publication', {})
            year = publication.get('year', '')
            
            if year:
                year_str = str(year)
                if year_str in publications_by_year:
                    publications_by_year[year_str] += 1
                else:
                    publications_by_year[year_str] = 1
        
        # Generate final summary
        summary = {
            'total_tools': total_tools,
            'tools_with_publications': tools_with_publications,
            'most_cited': most_cited,
            'publications_by_year': dict(sorted(publications_by_year.items())),
            'generated': datetime.now().isoformat(),
            'processing_time': (datetime.now() - datetime.fromisoformat(self.stats['start_time'])).total_seconds()
        }
        
        # Save summary
        with open(self.summary_file, 'w') as f:
            json.dump(summary, f, indent=2)
        
        return summary
    
    def generate_citations_md(self, results: Dict[str, Dict[str, Any]]) -> str:
        """Generate a markdown file with formatted citations for each tool."""
        logger.info("Generating citations markdown file")
        
        # Sort tools alphabetically by name
        sorted_tools = sorted(results.values(), key=lambda x: x.get('name', '').lower())
        
        # Build markdown content
        md_content = "# Academic Citations for Awesome-Virome Tools\n\n"
        md_content += "This document contains academic citations for tools in the Awesome-Virome collection. "
        md_content += "These citations can be used to properly acknowledge the tool authors in academic publications.\n\n"
        
        # Table of contents
        md_content += "## Table of Contents\n\n"
        for tool_data in sorted_tools:
            tool_name = tool_data.get('name', '')
            if tool_data.get('citation_info'):
                anchor = tool_name.lower().replace(' ', '-').replace('(', '').replace(')', '')
                md_content += f"- [{tool_name}](#{anchor})\n"
        
        # Generate citation sections for each tool
        for tool_data in sorted_tools:
            tool_name = tool_data.get('name', '')
            tool_url = tool_data.get('url', '')
            citation_info = tool_data.get('citation_info', {})
            
            if citation_info:
                md_content += f"\n\n## {tool_name}\n\n"
                md_content += f"**Tool URL**: [{tool_url}]({tool_url})\n\n"
                
                publication = citation_info.get('publication', {})
                formatted_citations = citation_info.get('formatted_citations', {})
                
                # Show publication info
                if isinstance(publication, dict) and publication:
                    md_content += "### Publication Information\n\n"
                    
                    title = publication.get('title', '')
                    journal = publication.get('journal', '')
                    year = publication.get('year', '')
                    doi = publication.get('doi', '')
                    pmid = publication.get('pmid', '')
                    
                    if title:
                        md_content += f"**Title**: {title}\n\n"
                    if journal:
                        md_content += f"**Journal**: {journal}\n\n"
                    if year:
                        md_content += f"**Year**: {year}\n\n"
                    
                    # Add identifiers
                    if doi or pmid:
                        md_content += "**Identifiers**:\n"
                        if doi:
                            md_content += f"- DOI: [{doi}](https://doi.org/{doi})\n"
                        if pmid:
                            md_content += f"- PMID: [{pmid}](https://pubmed.ncbi.nlm.nih.gov/{pmid}/)\n"
                        md_content += "\n"
                
                # Add formatted citations
                if formatted_citations:
                    md_content += "### Cite this tool\n\n"
                    
                    # APA style
                    if 'apa' in formatted_citations and formatted_citations['apa']:
                        md_content += "**APA format**:\n\n```\n"
                        md_content += formatted_citations['apa']
                        md_content += "\n```\n\n"
                    
                    # BibTeX format (for LaTeX/academic writing)
                    if 'bibtex' in formatted_citations and formatted_citations['bibtex']:
                        md_content += "**BibTeX format**:\n\n```bibtex\n"
                        md_content += formatted_citations['bibtex']
                        md_content += "\n```\n"
        
        # Write to file
        with open(CITATIONS_FILE, 'w') as f:
            f.write(md_content)
        
        logger.info(f"Generated citations markdown file: {CITATIONS_FILE}")
        return md_content
    
    def update_data_json(self, data_json_path: str = "data.json") -> None:
        """Update data.json with citation information."""
        logger.info(f"Updating citation information in {data_json_path}")
        
        # Load data.json
        try:
            with open(data_json_path, 'r') as f:
                data = json.load(f)
        except (json.JSONDecodeError, FileNotFoundError) as e:
            logger.error(f"Error loading {data_json_path}: {e}")
            return
        
        # Update nodes with citation information
        nodes = data.get('nodes', [])
        updated_count = 0
        
        for node in nodes:
            tool_name = node.get('name', '')
            repo_url = node.get('url', '')
            
            if tool_name and repo_url:
                cache_key = f"{tool_name}_{repo_url}"
                if cache_key in self.existing_metadata:
                    citation_info = self.existing_metadata[cache_key].get('citation_info', {})
                    if citation_info:
                        # Add citation information to the node
                        publication = citation_info.get('publication', {})
                        
                        if not node.get('citations'):
                            node['citations'] = {}
                        
                        # Add relevant citation data
                        if isinstance(publication, dict):
                            if 'doi' in publication:
                                node['citations']['doi'] = publication['doi']
                            if 'pmid' in publication:
                                node['citations']['pmid'] = publication['pmid']
                            if 'title' in publication:
                                node['citations']['title'] = publication['title']
                            if 'journal' in publication:
                                node['citations']['journal'] = publication['journal']
                            if 'year' in publication:
                                node['citations']['year'] = publication['year']
                            if 'citation_count' in publication:
                                node['citations']['citation_count'] = publication['citation_count']
                        
                        updated_count += 1
        
        logger.info(f"Updated citation information for {updated_count} tools")
        
        # Save updated data.json
        with open(data_json_path, 'w') as f:
            json.dump(data, f, indent=2)
        
        logger.info(f"Saved updated {data_json_path}")


def load_tools_data(data_json_path: str = "data.json") -> List[Dict[str, Any]]:
    """Load tools data from data.json."""
    try:
        with open(data_json_path, 'r') as f:
            data = json.load(f)
            return [node for node in data.get('nodes', []) if node.get('type') not in ['category', 'subcategory']]
    except (json.JSONDecodeError, FileNotFoundError) as e:
        logger.error(f"Error loading data.json: {e}")
        return []


def main():
    """Main function to run the PubMed citations collection."""
    parser = argparse.ArgumentParser(description='Collect PubMed citation data for Awesome-Virome tools')
    parser.add_argument('--output', help='Output directory for citation data', default=METADATA_DIR)
    parser.add_argument('--api-key', help='NCBI API key')
    parser.add_argument('--email', help='Contact email for API rate limiting')
    parser.add_argument('--github-token', help='GitHub API token')
    args = parser.parse_args()
    
    # Get API keys from environment variables if not provided as arguments
    ncbi_api_key = args.api_key or os.environ.get('NCBI_API_KEY')
    contact_email = args.email or os.environ.get('CONTACT_EMAIL')
    github_token = args.github_token or os.environ.get('GITHUB_TOKEN')
    
    # Warn if no API key provided
    if not ncbi_api_key:
        logger.warning("No NCBI API key provided. Rate limits will be restricted.")
    
    # Set up directory paths
    metadata_dir = args.output
    
    # Ensure directories exist
    os.makedirs(metadata_dir, exist_ok=True)
    
    # Load tools data
    tools = load_tools_data()
    if not tools:
        logger.error("No tools found in data.json")
        return
    
    logger.info(f"Loaded {len(tools)} tools from data.json")
    
    # Initialize collector
    collector = PubMedCitationsCollector(
        api_key=ncbi_api_key,
        email=contact_email,
        github_token=github_token,
        metadata_dir=metadata_dir
    )
    
    # Collect citation data
    start_time = time.time()
    results = collector.collect_pubmed_citations(tools)
    
    # Generate markdown file with citations
    collector.generate_citations_md(results)
    
    # Update data.json with citation information
    collector.update_data_json()
    
    # Print summary
    elapsed_time = time.time() - start_time
    logger.info(f"Citation data collection complete in {elapsed_time:.2f}s")
    logger.info(f"Processed {len(results)} tools")
    logger.info(f"Found publications for {collector.stats['tools_with_publications']} tools")
    logger.info(f"Generated citations markdown file: {CITATIONS_FILE}")


if __name__ == "__main__":
    main()