#!/usr/bin/env python3
"""
Academic Impact Metadata Collection for the Awesome-Virome Repository

This script collects academic impact and citation information for tools
in the Awesome-Virome collection, including:
- DOI identification from various sources
- Citation counts and metrics
- Related papers in the virology field
- Standardized citation formats

Usage:
    python academic_impact.py [--output OUTPUT] [--token TOKEN]
"""

import os
import json
import argparse
import logging
from datetime import datetime
from typing import Dict, List, Optional, Any, Tuple
from concurrent.futures import ThreadPoolExecutor
import re

from apis.citations_api import (
    ZenodoAPI, 
    SemanticScholarAPI, 
    CrossRefAPI, 
    GitHubAPI
)

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Directory constants
METADATA_DIR = os.path.join("metadata", "academic_impact")
CACHE_DIR = os.path.join("metadata", "cache")
OUTPUT_FILE = os.path.join(METADATA_DIR, "academic_impact.json")
SUMMARY_FILE = os.path.join(METADATA_DIR, "summary.json")


class AcademicImpactCollector:
    """Collects academic impact metadata for bioinformatics tools."""
    
    def __init__(self, github_token: Optional[str] = None, 
                 semantic_scholar_key: Optional[str] = None,
                 contact_email: Optional[str] = None,
                 zenodo_token: Optional[str] = None,
                 metadata_dir: str = METADATA_DIR,
                 cache_dir: str = CACHE_DIR):
        """Initialize the academic impact collector."""
        self.github_api = GitHubAPI(token=github_token)
        self.zenodo_api = ZenodoAPI(token=zenodo_token)
        self.semantic_scholar_api = SemanticScholarAPI(api_key=semantic_scholar_key)
        self.crossref_api = CrossRefAPI(email=contact_email)
        
        # Store directory paths as instance variables
        self.metadata_dir = metadata_dir
        self.cache_dir = cache_dir
        self.output_file = os.path.join(self.metadata_dir, "academic_impact.json")
        self.summary_file = os.path.join(self.metadata_dir, "summary.json")
        
        # Create directories
        os.makedirs(self.metadata_dir, exist_ok=True)
        os.makedirs(self.cache_dir, exist_ok=True)
        
        # Load existing metadata if available
        self.existing_metadata = {}
        if os.path.exists(self.output_file):
            try:
                with open(self.output_file, 'r') as f:
                    self.existing_metadata = json.load(f)
            except json.JSONDecodeError:
                logger.warning(f"Could not parse existing metadata from {self.output_file}")
    
    def get_repo_citation_info(self, tool_name: str, repo_url: str) -> Dict[str, Any]:
        """Get citation information from repository files like CITATION.cff."""
        if not repo_url or 'github.com' not in repo_url:
            return {}
        
        # Check if we already have this data
        cache_key = f"{tool_name}_{repo_url}"
        if cache_key in self.existing_metadata:
            last_updated = self.existing_metadata[cache_key].get('last_updated', '')
            # If updated within the last 30 days, use cached data
            if last_updated and (datetime.now() - datetime.fromisoformat(last_updated)).days < 30:
                logger.info(f"Using cached citation info for {tool_name}")
                return self.existing_metadata[cache_key]
        
        logger.info(f"Fetching repository citation info for {tool_name}")
        citation_file = self.github_api.find_citation_file(repo_url)
        
        if not citation_file:
            logger.info(f"No citation file found for {tool_name}")
            return {}
        
        citation_data = {}
        
        # Parse according to file format
        if citation_file['format'] == 'cff':
            citation_data = self.github_api.parse_citation_cff(citation_file['content'])
        else:
            # For other formats, just extract DOI if present
            doi = self.github_api.extract_doi_from_text(citation_file['content'])
            if doi:
                citation_data['doi'] = doi
        
        return {
            'source': 'github_citation_file',
            'citation_file': citation_file['path'],
            'citation_format': citation_file['format'],
            'data': citation_data,
            'last_updated': datetime.now().isoformat()
        }
    
    def find_doi(self, tool_name: str, repo_url: str, description: str = "") -> Optional[str]:
        """Find DOI for a tool checking multiple sources."""
        # First check if we already have a DOI from repository citation files
        citation_info = self.get_repo_citation_info(tool_name, repo_url)
        doi = citation_info.get('data', {}).get('doi')
        
        if doi:
            logger.info(f"Found DOI from citation file for {tool_name}: {doi}")
            return doi
        
        # Try to find DOI in repository README or other documentation
        if repo_url and 'github.com' in repo_url:
            match = re.search(r'github\.com/([^/]+)/([^/]+)', repo_url)
            if match:
                owner, repo = match.groups()
                repo = repo.rstrip('.git')
                
                readme = self.github_api.get_file_content(owner, repo, "README.md")
                if readme:
                    doi = self.github_api.extract_doi_from_text(readme)
                    if doi:
                        logger.info(f"Found DOI in README for {tool_name}: {doi}")
                        return doi
        
        # Search Zenodo for DOI
        doi = self.zenodo_api.find_tool_doi(tool_name, repo_url)
        if doi:
            logger.info(f"Found DOI from Zenodo for {tool_name}: {doi}")
            return doi
        
        # Search CrossRef using tool name and description
        search_terms = [
            tool_name,
            f"{tool_name} software",
            f"{tool_name} bioinformatics",
            f"{tool_name} virology"
        ]
        
        for term in search_terms:
            results = self.crossref_api.search_works(term)
            for work in results:
                title = work.get('title', [''])[0].lower()
                if tool_name.lower() in title and (
                    'software' in title or
                    'tool' in title or
                    'bioinformatics' in title or
                    'pipeline' in title or
                    'workflow' in title
                ):
                    doi = work.get('DOI')
                    if doi:
                        logger.info(f"Found DOI from CrossRef for {tool_name}: {doi}")
                        return doi
        
        logger.info(f"No DOI found for {tool_name}")
        return None
    
    def get_citation_metrics(self, doi: str) -> Dict[str, Any]:
        """Get citation metrics for a DOI."""
        if not doi:
            return {}
        
        # Get publication details from CrossRef
        publication = self.crossref_api.get_work_by_doi(doi)
        
        # Get citation metrics from Semantic Scholar
        paper = self.semantic_scholar_api.get_paper_by_doi(doi)
        
        if not paper:
            logger.warning(f"Publication not found in Semantic Scholar for DOI: {doi}")
            return {
                'publication': publication,
                'metrics': {
                    'total_citations': 0,
                    'influential_citations': 0,
                    'citations_by_year': {}
                },
                'formatted_citations': {}
            }
        
        paper_id = paper.get('paperId')
        metrics = self.semantic_scholar_api.get_citation_metrics(paper_id)
        
        # Create formatted citations
        formatted_citations = {}
        if publication:
            formatted_citations['apa'] = self.crossref_api.format_citation(publication, 'apa')
            formatted_citations['bibtex'] = self.crossref_api.format_citation(publication, 'bibtex')
        
        return {
            'publication': publication,
            'metrics': metrics,
            'formatted_citations': formatted_citations,
            'paper_id': paper_id
        }
    
    def get_related_papers(self, paper_id: str, field: str = "virology") -> List[Dict[str, Any]]:
        """Get related papers in the virology field."""
        if not paper_id:
            return []
        
        return self.semantic_scholar_api.find_related_papers(paper_id, field)
    
    def process_tool(self, tool: Dict[str, Any]) -> Dict[str, Any]:
        """Process a single tool to collect academic impact metadata."""
        tool_name = tool.get('name', '')
        repo_url = tool.get('url', '')
        description = tool.get('description', '')
        
        logger.info(f"Processing academic impact data for {tool_name}")
        
        # Check if we already have recent data for this tool
        cache_key = f"{tool_name}_{repo_url}"
        if cache_key in self.existing_metadata:
            last_updated = self.existing_metadata[cache_key].get('last_updated', '')
            # If updated within the last 30 days, use cached data
            if last_updated and (datetime.now() - datetime.fromisoformat(last_updated)).days < 30:
                logger.info(f"Using cached academic impact data for {tool_name}")
                return self.existing_metadata[cache_key]
        
        # Find DOI for the tool
        doi = self.find_doi(tool_name, repo_url, description)
        
        # Get citation information
        citation_info = self.get_repo_citation_info(tool_name, repo_url)
        
        # Get citation metrics if DOI is available
        citation_metrics = {}
        related_papers = []
        if doi:
            citation_metrics = self.get_citation_metrics(doi)
            
            # Get related papers if we have a paper ID
            paper_id = citation_metrics.get('paper_id')
            if paper_id:
                related_papers = self.get_related_papers(paper_id)
        
        # Compile all academic impact data
        academic_impact = {
            'name': tool_name,
            'url': repo_url,
            'doi': doi,
            'citation_info': citation_info,
            'citation_metrics': citation_metrics,
            'related_papers': related_papers,
            'last_updated': datetime.now().isoformat()
        }
        
        # Save individual tool data
        tool_file = os.path.join(self.metadata_dir, f"{tool_name.replace('/', '_')}.json")
        with open(tool_file, 'w') as f:
            json.dump(academic_impact, f, indent=2)
        
        return academic_impact
    
    def collect_academic_impact(self, tools: List[Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
        """Collect academic impact data for all tools."""
        results = {}
        
        # Process tools in parallel
        with ThreadPoolExecutor(max_workers=5) as executor:
            for tool, result in zip(tools, executor.map(self.process_tool, tools)):
                tool_name = tool.get('name', '')
                repo_url = tool.get('url', '')
                cache_key = f"{tool_name}_{repo_url}"
                results[cache_key] = result
        
        # Save all results
        with open(self.output_file, 'w') as f:
            json.dump(results, f, indent=2)
        
        return results
    
    def generate_summary(self, results: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
        """Generate a summary of academic impact data."""
        total_tools = len(results)
        tools_with_doi = sum(1 for data in results.values() if data.get('doi'))
        tools_with_citations = sum(1 for data in results.values() 
                                if data.get('citation_metrics', {}).get('metrics', {}).get('total_citations', 0) > 0)
        
        # Calculate total citations
        total_citations = sum(data.get('citation_metrics', {}).get('metrics', {}).get('total_citations', 0) 
                            for data in results.values())
        
        # Find most cited tools
        most_cited = sorted(
            [(key, data.get('name', ''), data.get('citation_metrics', {}).get('metrics', {}).get('total_citations', 0)) 
             for key, data in results.items()],
            key=lambda x: x[2],
            reverse=True
        )[:10]
        
        # Citations by year
        citations_by_year = {}
        for data in results.values():
            yearly_citations = data.get('citation_metrics', {}).get('metrics', {}).get('citations_by_year', {})
            for year, count in yearly_citations.items():
                # Convert year to string to ensure consistent type for keys
                year_key = str(year)
                if year_key in citations_by_year:
                    citations_by_year[year_key] += count
                else:
                    citations_by_year[year_key] = count
        
        summary = {
            'total_tools': total_tools,
            'tools_with_doi': tools_with_doi,
            'tools_with_citations': tools_with_citations,
            'total_citations': total_citations,
            'most_cited_tools': most_cited,
            'citations_by_year': dict(sorted(citations_by_year.items())),
            'generated': datetime.now().isoformat()
        }
        
        # Save summary
        with open(self.summary_file, 'w') as f:
            json.dump(summary, f, indent=2)
        
        return summary


def load_tools_data() -> List[Dict[str, Any]]:
    """Load tools data from data.json."""
    try:
        with open('data.json', 'r') as f:
            data = json.load(f)
            return data.get('nodes', [])
    except (json.JSONDecodeError, FileNotFoundError) as e:
        logger.error(f"Error loading data.json: {e}")
        return []


def main():
    """Main function to run the academic impact data collection."""
    parser = argparse.ArgumentParser(description='Collect academic impact data for Awesome-Virome tools')
    parser.add_argument('--output', help='Output directory for academic impact data', default=METADATA_DIR)
    parser.add_argument('--github-token', help='GitHub API token')
    parser.add_argument('--semantic-scholar-key', help='Semantic Scholar API key')
    parser.add_argument('--contact-email', help='Contact email for API rate limiting')
    parser.add_argument('--zenodo-token', help='Zenodo API token')
    args = parser.parse_args()
    
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
    
    # Initialize collector with custom directory paths
    collector = AcademicImpactCollector(
        github_token=args.github_token,
        semantic_scholar_key=args.semantic_scholar_key,
        contact_email=args.contact_email,
        zenodo_token=args.zenodo_token,
        metadata_dir=metadata_dir
    )
    
    # Collect academic impact data
    results = collector.collect_academic_impact(tools)
    
    # Generate summary
    summary = collector.generate_summary(results)
    
    logger.info(f"Academic impact data collection complete. Processed {len(results)} tools.")
    logger.info(f"Found DOIs for {summary['tools_with_doi']} tools.")
    logger.info(f"Found citations for {summary['tools_with_citations']} tools.")
    logger.info(f"Total citations: {summary['total_citations']}")


if __name__ == "__main__":
    main()