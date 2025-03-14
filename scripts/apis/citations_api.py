#!/usr/bin/env python3
"""
Module for integrating with various citation and academic impact APIs.
"""

import os
import json
import time
import requests
import logging
import base64
import re
from datetime import datetime, timedelta
from typing import Dict, List, Optional, Any, Tuple, Union
from urllib.parse import quote

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class RateLimiter:
    """Simple rate limiter to prevent API rate limit issues."""
    
    def __init__(self, calls_per_minute: int = 10):
        self.calls_per_minute = calls_per_minute
        self.interval = 60 / calls_per_minute
        self.last_call = 0
    
    def wait(self):
        """Wait if necessary to respect the rate limit."""
        elapsed = time.time() - self.last_call
        if elapsed < self.interval:
            time.sleep(self.interval - elapsed)
        self.last_call = time.time()


class ZenodoAPI:
    """Integration with Zenodo for DOI information."""
    
    BASE_URL = "https://zenodo.org/api"
    CACHE_DIR = os.path.join("metadata", "cache", "zenodo")
    CACHE_EXPIRY = 30  # days
    
    def __init__(self, token: Optional[str] = None):
        self.token = token
        self.rate_limiter = RateLimiter(calls_per_minute=20)
        os.makedirs(self.CACHE_DIR, exist_ok=True)
    
    def _get_cache_path(self, query: str) -> str:
        """Get the path for a cached result."""
        safe_query = base64.urlsafe_b64encode(query.encode()).decode()
        return os.path.join(self.CACHE_DIR, f"{safe_query}.json")
    
    def _cache_valid(self, cache_path: str) -> bool:
        """Check if the cache is still valid."""
        if not os.path.exists(cache_path):
            return False
        
        try:
            with open(cache_path, 'r') as f:
                data = json.load(f)
                cache_date = datetime.fromisoformat(data.get('cache_date', '2000-01-01'))
                return datetime.now() - cache_date < timedelta(days=self.CACHE_EXPIRY)
        except (json.JSONDecodeError, ValueError, KeyError):
            return False
    
    def _save_to_cache(self, cache_path: str, data: Dict[str, Any]) -> None:
        """Save API response to cache."""
        cache_data = {
            'cache_date': datetime.now().isoformat(),
            'data': data
        }
        with open(cache_path, 'w') as f:
            json.dump(cache_data, f, indent=2)
    
    def _get_from_cache(self, cache_path: str) -> Dict[str, Any]:
        """Get data from cache."""
        with open(cache_path, 'r') as f:
            return json.load(f)['data']
    
    def search_records(self, query: str) -> List[Dict[str, Any]]:
        """Search Zenodo records for the given query."""
        cache_path = self._get_cache_path(f"search_{query}")
        
        if self._cache_valid(cache_path):
            logger.info(f"Using cached Zenodo results for '{query}'")
            return self._get_from_cache(cache_path)
        
        self.rate_limiter.wait()
        url = f"{self.BASE_URL}/records"
        params = {
            'q': query,
            'sort': 'mostrecent',
            'size': 10
        }
        
        headers = {}
        if self.token:
            headers['Authorization'] = f"Bearer {self.token}"
        
        try:
            response = requests.get(url, params=params, headers=headers)
            response.raise_for_status()
            results = response.json().get('hits', {}).get('hits', [])
            self._save_to_cache(cache_path, results)
            return results
        except requests.RequestException as e:
            logger.error(f"Error searching Zenodo: {e}")
            return []
    
    def get_doi_metadata(self, doi: str) -> Dict[str, Any]:
        """Get metadata for a specific DOI."""
        cache_path = self._get_cache_path(f"doi_{doi}")
        
        if self._cache_valid(cache_path):
            logger.info(f"Using cached Zenodo DOI metadata for '{doi}'")
            return self._get_from_cache(cache_path)
        
        self.rate_limiter.wait()
        url = f"{self.BASE_URL}/records"
        params = {'q': f"doi:\"{doi}\""} 
        
        headers = {}
        if self.token:
            headers['Authorization'] = f"Bearer {self.token}"
        
        try:
            response = requests.get(url, params=params, headers=headers)
            response.raise_for_status()
            results = response.json().get('hits', {}).get('hits', [])
            if results:
                self._save_to_cache(cache_path, results[0])
                return results[0]
            return {}
        except requests.RequestException as e:
            logger.error(f"Error getting DOI metadata from Zenodo: {e}")
            return {}
    
    def find_tool_doi(self, tool_name: str, repo_url: str) -> Optional[str]:
        """Find a DOI for a tool based on name and repository URL."""
        search_terms = [
            tool_name,
            tool_name.replace('-', ' '),
            f"{tool_name} software",
            f"{tool_name} bioinformatics"
        ]
        
        for term in search_terms:
            results = self.search_records(term)
            for record in results:
                metadata = record.get('metadata', {})
                # Check if this record matches our tool
                description = metadata.get('description', '').lower()
                related_ids = metadata.get('related_identifiers', [])
                
                # Look for GitHub URL in related identifiers
                for related in related_ids:
                    if related.get('relation', '') == 'isSupplementTo':
                        identifier = related.get('identifier', '')
                        if repo_url.lower() in identifier.lower():
                            return metadata.get('doi')
                
                # Check if tool name is in title
                title = metadata.get('title', '').lower()
                if tool_name.lower() in title and (
                    'software' in title or 
                    'tool' in title or 
                    'bioinformatics' in title or
                    'pipeline' in title or
                    'workflow' in title
                ):
                    return metadata.get('doi')
        
        return None


class SemanticScholarAPI:
    """Integration with Semantic Scholar API for citation metrics."""
    
    BASE_URL = "https://api.semanticscholar.org/graph/v1"
    CACHE_DIR = os.path.join("metadata", "cache", "semanticscholar")
    CACHE_EXPIRY = 30  # days
    
    def __init__(self, api_key: Optional[str] = None):
        self.api_key = api_key
        self.rate_limiter = RateLimiter(calls_per_minute=10)
        os.makedirs(self.CACHE_DIR, exist_ok=True)
    
    def _get_cache_path(self, query: str) -> str:
        """Get the path for a cached result."""
        safe_query = base64.urlsafe_b64encode(query.encode()).decode()
        return os.path.join(self.CACHE_DIR, f"{safe_query}.json")
    
    def _cache_valid(self, cache_path: str) -> bool:
        """Check if the cache is still valid."""
        if not os.path.exists(cache_path):
            return False
        
        try:
            with open(cache_path, 'r') as f:
                data = json.load(f)
                cache_date = datetime.fromisoformat(data.get('cache_date', '2000-01-01'))
                return datetime.now() - cache_date < timedelta(days=self.CACHE_EXPIRY)
        except (json.JSONDecodeError, ValueError, KeyError):
            return False
    
    def _save_to_cache(self, cache_path: str, data: Dict[str, Any]) -> None:
        """Save API response to cache."""
        cache_data = {
            'cache_date': datetime.now().isoformat(),
            'data': data
        }
        with open(cache_path, 'w') as f:
            json.dump(cache_data, f, indent=2)
    
    def _get_from_cache(self, cache_path: str) -> Dict[str, Any]:
        """Get data from cache."""
        with open(cache_path, 'r') as f:
            return json.load(f)['data']
    
    def search_paper(self, query: str) -> List[Dict[str, Any]]:
        """Search for papers matching a query."""
        cache_path = self._get_cache_path(f"search_{query}")
        
        if self._cache_valid(cache_path):
            logger.info(f"Using cached Semantic Scholar results for '{query}'")
            return self._get_from_cache(cache_path)
        
        self.rate_limiter.wait()
        url = f"{self.BASE_URL}/paper/search"
        params = {
            'query': query,
            'limit': 10,
            'fields': 'title,authors,year,venue,citationCount,externalIds,url,abstract'
        }
        
        headers = {}
        if self.api_key:
            headers['x-api-key'] = self.api_key
        
        try:
            response = requests.get(url, params=params, headers=headers)
            response.raise_for_status()
            results = response.json().get('data', [])
            self._save_to_cache(cache_path, results)
            return results
        except requests.RequestException as e:
            logger.error(f"Error searching Semantic Scholar: {e}")
            return []
    
    def get_paper_by_doi(self, doi: str) -> Dict[str, Any]:
        """Get paper information using DOI."""
        cache_path = self._get_cache_path(f"doi_{doi}")
        
        if self._cache_valid(cache_path):
            logger.info(f"Using cached Semantic Scholar DOI data for '{doi}'")
            return self._get_from_cache(cache_path)
        
        self.rate_limiter.wait()
        url = f"{self.BASE_URL}/paper/DOI:{doi}"
        params = {
            'fields': 'title,authors,year,venue,citationCount,referenceCount,influentialCitationCount,externalIds,url,abstract,references,citations'
        }
        
        headers = {}
        if self.api_key:
            headers['x-api-key'] = self.api_key
        
        try:
            response = requests.get(url, params=params, headers=headers)
            response.raise_for_status()
            data = response.json()
            self._save_to_cache(cache_path, data)
            return data
        except requests.RequestException as e:
            logger.error(f"Error getting paper by DOI from Semantic Scholar: {e}")
            return {}
    
    def get_paper_by_title_author(self, title: str, authors: List[str]) -> Dict[str, Any]:
        """Find a paper by its title and authors."""
        # Generate a search query with title and first author
        first_author = authors[0] if authors else ""
        query = f"{title} {first_author}"
        
        results = self.search_paper(query)
        
        # Try to find an exact title match
        for paper in results:
            if paper.get('title', '').lower() == title.lower():
                # Further validate with author matching
                paper_authors = [author.get('name', '') for author in paper.get('authors', [])]
                # Check if at least one author matches
                for author in authors:
                    if any(author.lower() in paper_author.lower() for paper_author in paper_authors):
                        return paper
        
        # If no exact match, return the top result if it seems relevant
        if results and title.lower() in results[0].get('title', '').lower():
            return results[0]
        
        return {}
    
    def get_citation_metrics(self, paper_id: str) -> Dict[str, Any]:
        """Get detailed citation metrics for a paper."""
        cache_path = self._get_cache_path(f"metrics_{paper_id}")
        
        if self._cache_valid(cache_path):
            logger.info(f"Using cached Semantic Scholar metrics for '{paper_id}'")
            return self._get_from_cache(cache_path)
        
        self.rate_limiter.wait()
        url = f"{self.BASE_URL}/paper/{paper_id}"
        params = {
            'fields': 'citationCount,influentialCitationCount,citations.year'
        }
        
        headers = {}
        if self.api_key:
            headers['x-api-key'] = self.api_key
        
        try:
            response = requests.get(url, params=params, headers=headers)
            response.raise_for_status()
            data = response.json()
            
            # Process citation data to get citations by year
            citations_by_year = {}
            for citation in data.get('citations', []):
                year = citation.get('year')
                if year:
                    citations_by_year[year] = citations_by_year.get(year, 0) + 1
            
            metrics = {
                'total_citations': data.get('citationCount', 0),
                'influential_citations': data.get('influentialCitationCount', 0),
                'citations_by_year': citations_by_year
            }
            
            self._save_to_cache(cache_path, metrics)
            return metrics
        except requests.RequestException as e:
            logger.error(f"Error getting citation metrics from Semantic Scholar: {e}")
            return {
                'total_citations': 0,
                'influential_citations': 0,
                'citations_by_year': {}
            }
    
    def find_related_papers(self, paper_id: str, field: str = "virology") -> List[Dict[str, Any]]:
        """Find related papers in a specific field."""
        cache_path = self._get_cache_path(f"related_{paper_id}_{field}")
        
        if self._cache_valid(cache_path):
            logger.info(f"Using cached related papers for '{paper_id}' in '{field}'")
            return self._get_from_cache(cache_path)
        
        self.rate_limiter.wait()
        url = f"{self.BASE_URL}/paper/{paper_id}/related"
        params = {
            'fields': 'title,authors,year,venue,citationCount,url,abstract',
            'limit': 20
        }
        
        headers = {}
        if self.api_key:
            headers['x-api-key'] = self.api_key
        
        try:
            response = requests.get(url, params=params, headers=headers)
            response.raise_for_status()
            all_related = response.json().get('data', [])
            
            # Filter for papers related to the specified field
            field_related = []
            for paper in all_related:
                abstract = paper.get('abstract', '').lower()
                title = paper.get('title', '').lower()
                venue = paper.get('venue', '').lower()
                
                if (field.lower() in abstract or
                    field.lower() in title or
                    'virus' in abstract or
                    'viral' in abstract or
                    'pathogen' in abstract):
                    field_related.append(paper)
            
            self._save_to_cache(cache_path, field_related[:10])  # Store top 10
            return field_related[:10]
        except requests.RequestException as e:
            logger.error(f"Error finding related papers from Semantic Scholar: {e}")
            return []


class CrossRefAPI:
    """Integration with CrossRef API for publication information."""
    
    BASE_URL = "https://api.crossref.org"
    CACHE_DIR = os.path.join("metadata", "cache", "crossref")
    CACHE_EXPIRY = 30  # days
    
    def __init__(self, email: Optional[str] = None):
        """Initialize CrossRef API client with optional email for Polite Pool."""
        self.email = email
        self.rate_limiter = RateLimiter(calls_per_minute=10)
        os.makedirs(self.CACHE_DIR, exist_ok=True)
    
    def _get_cache_path(self, query: str) -> str:
        """Get the path for a cached result."""
        safe_query = base64.urlsafe_b64encode(query.encode()).decode()
        return os.path.join(self.CACHE_DIR, f"{safe_query}.json")
    
    def _cache_valid(self, cache_path: str) -> bool:
        """Check if the cache is still valid."""
        if not os.path.exists(cache_path):
            return False
        
        try:
            with open(cache_path, 'r') as f:
                data = json.load(f)
                cache_date = datetime.fromisoformat(data.get('cache_date', '2000-01-01'))
                return datetime.now() - cache_date < timedelta(days=self.CACHE_EXPIRY)
        except (json.JSONDecodeError, ValueError, KeyError):
            return False
    
    def _save_to_cache(self, cache_path: str, data: Dict[str, Any]) -> None:
        """Save API response to cache."""
        cache_data = {
            'cache_date': datetime.now().isoformat(),
            'data': data
        }
        with open(cache_path, 'w') as f:
            json.dump(cache_data, f, indent=2)
    
    def _get_from_cache(self, cache_path: str) -> Dict[str, Any]:
        """Get data from cache."""
        with open(cache_path, 'r') as f:
            return json.load(f)['data']
    
    def search_works(self, query: str) -> List[Dict[str, Any]]:
        """Search CrossRef works with the given query."""
        cache_path = self._get_cache_path(f"search_{query}")
        
        if self._cache_valid(cache_path):
            logger.info(f"Using cached CrossRef results for '{query}'")
            return self._get_from_cache(cache_path)
        
        self.rate_limiter.wait()
        url = f"{self.BASE_URL}/works"
        params = {
            'query': query,
            'rows': 10,
            'sort': 'relevance'
        }
        
        headers = {}
        if self.email:
            headers['User-Agent'] = f"AwesomeVirome/1.0 (mailto:{self.email})"
        
        try:
            response = requests.get(url, params=params, headers=headers)
            response.raise_for_status()
            results = response.json().get('message', {}).get('items', [])
            self._save_to_cache(cache_path, results)
            return results
        except requests.RequestException as e:
            logger.error(f"Error searching CrossRef: {e}")
            return []
    
    def get_work_by_doi(self, doi: str) -> Dict[str, Any]:
        """Get detailed work information by DOI."""
        cache_path = self._get_cache_path(f"doi_{doi}")
        
        if self._cache_valid(cache_path):
            logger.info(f"Using cached CrossRef DOI data for '{doi}'")
            return self._get_from_cache(cache_path)
        
        self.rate_limiter.wait()
        url = f"{self.BASE_URL}/works/{quote(doi)}"
        
        headers = {}
        if self.email:
            headers['User-Agent'] = f"AwesomeVirome/1.0 (mailto:{self.email})"
        
        try:
            response = requests.get(url, headers=headers)
            response.raise_for_status()
            data = response.json().get('message', {})
            self._save_to_cache(cache_path, data)
            return data
        except requests.RequestException as e:
            logger.error(f"Error getting work by DOI from CrossRef: {e}")
            return {}
    
    def format_citation(self, work: Dict[str, Any], style: str = "apa") -> str:
        """Format a citation in the specified style."""
        if not work:
            return ""
        
        try:
            # Extract publication data
            title = work.get('title', [''])[0]
            authors = []
            for author in work.get('author', []):
                given = author.get('given', '')
                family = author.get('family', '')
                if given and family:
                    authors.append(f"{family}, {given[0]}.")
                elif family:
                    authors.append(family)
            
            year = work.get('published', {}).get('date-parts', [['']])[0][0]
            journal = work.get('container-title', [''])[0]
            volume = work.get('volume', '')
            issue = work.get('issue', '')
            pages = work.get('page', '')
            doi = work.get('DOI', '')
            
            if style.lower() == "apa":
                # Format authors
                if len(authors) == 0:
                    author_text = ""
                elif len(authors) == 1:
                    author_text = f"{authors[0]}"
                elif len(authors) < 8:
                    author_text = ", ".join(authors[:-1]) + f", & {authors[-1]}"
                else:
                    author_text = ", ".join(authors[:6]) + f", ... {authors[-1]}"
                
                # Format APA style citation
                citation = f"{author_text} ({year}). {title}."
                if journal:
                    citation += f" {journal}"
                    if volume:
                        citation += f", {volume}"
                        if issue:
                            citation += f"({issue})"
                    if pages:
                        citation += f", {pages}"
                citation += f". https://doi.org/{doi}" if doi else ""
                
                return citation
            
            elif style.lower() == "bibtex":
                # Create a BibTeX key from first author's last name and year
                first_author = work.get('author', [{}])[0].get('family', 'Unknown').lower()
                bibtex_key = f"{first_author}{year}"
                
                # Format BibTeX entry
                bibtex = f"@article{{{bibtex_key},\n"
                bibtex += f"  title = {{{title}}},\n"
                
                # Authors
                if authors:
                    bibtex += f"  author = {{{' and '.join(authors)}}},\n"
                
                # Other fields
                if year:
                    bibtex += f"  year = {{{year}}},\n"
                if journal:
                    bibtex += f"  journal = {{{journal}}},\n"
                if volume:
                    bibtex += f"  volume = {{{volume}}},\n"
                if issue:
                    bibtex += f"  number = {{{issue}}},\n"
                if pages:
                    bibtex += f"  pages = {{{pages}}},\n"
                if doi:
                    bibtex += f"  doi = {{{doi}}},\n"
                
                bibtex += "}"
                return bibtex
            
            # Default to a simple citation format
            return f"{', '.join(authors[:3])} {'et al. ' if len(authors) > 3 else ''}({year}). {title}. {journal}. DOI: {doi}"
        
        except Exception as e:
            logger.error(f"Error formatting citation: {e}")
            return ""


class GitHubAPI:
    """Integration with GitHub API to extract citation information."""
    
    BASE_URL = "https://api.github.com"
    CACHE_DIR = os.path.join("metadata", "cache", "github")
    CACHE_EXPIRY = 7  # days
    
    def __init__(self, token: Optional[str] = None):
        self.token = token
        self.rate_limiter = RateLimiter(calls_per_minute=30)
        os.makedirs(self.CACHE_DIR, exist_ok=True)
    
    def _get_cache_path(self, query: str) -> str:
        """Get the path for a cached result."""
        safe_query = base64.urlsafe_b64encode(query.encode()).decode()
        return os.path.join(self.CACHE_DIR, f"{safe_query}.json")
    
    def _cache_valid(self, cache_path: str) -> bool:
        """Check if the cache is still valid."""
        if not os.path.exists(cache_path):
            return False
        
        try:
            with open(cache_path, 'r') as f:
                data = json.load(f)
                cache_date = datetime.fromisoformat(data.get('cache_date', '2000-01-01'))
                return datetime.now() - cache_date < timedelta(days=self.CACHE_EXPIRY)
        except (json.JSONDecodeError, ValueError, KeyError):
            return False
    
    def _save_to_cache(self, cache_path: str, data: Dict[str, Any]) -> None:
        """Save API response to cache."""
        cache_data = {
            'cache_date': datetime.now().isoformat(),
            'data': data
        }
        with open(cache_path, 'w') as f:
            json.dump(cache_data, f, indent=2)
    
    def _get_from_cache(self, cache_path: str) -> Dict[str, Any]:
        """Get data from cache."""
        with open(cache_path, 'r') as f:
            return json.load(f)['data']
    
    def get_repo_contents(self, owner: str, repo: str, path: str = "") -> List[Dict[str, Any]]:
        """Get contents of a repository at the given path."""
        cache_path = self._get_cache_path(f"contents_{owner}_{repo}_{path}")
        
        if self._cache_valid(cache_path):
            logger.info(f"Using cached GitHub contents for '{owner}/{repo}/{path}'")
            return self._get_from_cache(cache_path)
        
        self.rate_limiter.wait()
        url = f"{self.BASE_URL}/repos/{owner}/{repo}/contents/{path}"
        
        headers = {
            'Accept': 'application/vnd.github.v3+json'
        }
        if self.token:
            headers['Authorization'] = f"token {self.token}"
        
        try:
            response = requests.get(url, headers=headers)
            response.raise_for_status()
            contents = response.json()
            self._save_to_cache(cache_path, contents)
            return contents if isinstance(contents, list) else [contents]
        except requests.RequestException as e:
            logger.error(f"Error getting repo contents from GitHub: {e}")
            return []
    
    def get_file_content(self, owner: str, repo: str, path: str) -> Optional[str]:
        """Get the content of a specific file."""
        cache_path = self._get_cache_path(f"file_{owner}_{repo}_{path}")
        
        if self._cache_valid(cache_path):
            logger.info(f"Using cached GitHub file for '{owner}/{repo}/{path}'")
            return self._get_from_cache(cache_path)
        
        self.rate_limiter.wait()
        url = f"{self.BASE_URL}/repos/{owner}/{repo}/contents/{path}"
        
        headers = {
            'Accept': 'application/vnd.github.v3+json'
        }
        if self.token:
            headers['Authorization'] = f"token {self.token}"
        
        try:
            response = requests.get(url, headers=headers)
            response.raise_for_status()
            data = response.json()
            
            if isinstance(data, dict) and 'content' in data:
                content = base64.b64decode(data['content']).decode('utf-8')
                self._save_to_cache(cache_path, content)
                return content
            return None
        except requests.RequestException as e:
            logger.error(f"Error getting file content from GitHub: {e}")
            return None
    
    def find_citation_file(self, repo_url: str) -> Optional[Dict[str, Any]]:
        """Find CITATION.cff or other citation files in a repository."""
        # Extract owner and repo from URL
        match = re.search(r'github\.com/([^/]+)/([^/]+)', repo_url)
        if not match:
            logger.error(f"Could not extract owner/repo from URL: {repo_url}")
            return None
        
        owner, repo = match.groups()
        repo = repo.rstrip('.git')
        
        # Common citation file paths to check
        citation_paths = [
            "CITATION.cff",
            "CITATION",
            "citation.cff",
            "citation.txt",
            "CITATION.txt",
            "CITATION.md",
            "citation.md",
            "CITATION.bib",
            "citation.bib"
        ]
        
        for path in citation_paths:
            content = self.get_file_content(owner, repo, path)
            if content:
                return {
                    'path': path,
                    'content': content,
                    'format': path.split('.')[-1] if '.' in path else 'txt'
                }
        
        # Look in .github directory
        for path in [f".github/{p}" for p in citation_paths]:
            content = self.get_file_content(owner, repo, path)
            if content:
                return {
                    'path': path,
                    'content': content,
                    'format': path.split('.')[-1] if '.' in path else 'txt'
                }
        
        # Check README for citation information
        readme_content = self.get_file_content(owner, repo, "README.md")
        if readme_content:
            # Look for citation section in README
            citation_section = None
            citation_patterns = [
                r'(?i)#+\s*citation\s*\n+(.*?)(?:\n#+\s|\Z)',
                r'(?i)#+\s*how\s+to\s+cite\s*\n+(.*?)(?:\n#+\s|\Z)',
                r'(?i)#+\s*citing\s+this\s+work\s*\n+(.*?)(?:\n#+\s|\Z)'
            ]
            
            for pattern in citation_patterns:
                match = re.search(pattern, readme_content, re.DOTALL)
                if match:
                    citation_section = match.group(1).strip()
                    break
            
            if citation_section:
                return {
                    'path': 'README.md',
                    'content': citation_section,
                    'format': 'md'
                }
        
        return None
    
    def parse_citation_cff(self, content: str) -> Dict[str, Any]:
        """Parse CITATION.cff file for citation information."""
        try:
            import yaml
            citation_data = yaml.safe_load(content)
            
            # Extract relevant citation information
            citation = {
                'title': citation_data.get('title', ''),
                'authors': [
                    f"{author.get('family-names', '')} {author.get('given-names', '')}"
                    for author in citation_data.get('authors', [])
                ],
                'year': citation_data.get('date-released', '').split('-')[0] if citation_data.get('date-released') else '',
                'doi': citation_data.get('doi', ''),
                'version': citation_data.get('version', ''),
                'url': citation_data.get('url', ''),
                'repository': citation_data.get('repository-code', '')
            }
            
            # Format according to specified styles
            formatted = {}
            if citation_data.get('preferred-citation'):
                preferred = citation_data.get('preferred-citation', {})
                if 'doi' in preferred:
                    # If the preferred citation has a DOI, we can get formatted citations later
                    citation['preferred_citation_doi'] = preferred.get('doi')
                else:
                    # Otherwise, build from the preferred citation data
                    authors = []
                    for author in preferred.get('authors', []):
                        given = author.get('given-names', '')
                        family = author.get('family-names', '')
                        if given and family:
                            authors.append(f"{family}, {given[0]}.")
                        elif family:
                            authors.append(family)
                    
                    journal = preferred.get('journal', '')
                    volume = preferred.get('volume', '')
                    issue = preferred.get('issue', '')
                    pages = preferred.get('pages', '')
                    year = preferred.get('year', '') or citation['year']
                    
                    # APA style
                    if authors:
                        if len(authors) == 1:
                            apa_authors = authors[0]
                        elif len(authors) < 8:
                            apa_authors = ", ".join(authors[:-1]) + f", & {authors[-1]}"
                        else:
                            apa_authors = ", ".join(authors[:6]) + f", ... {authors[-1]}"
                        
                        apa = f"{apa_authors} ({year}). {preferred.get('title', '')}."
                        if journal:
                            apa += f" {journal}"
                            if volume:
                                apa += f", {volume}"
                                if issue:
                                    apa += f"({issue})"
                            if pages:
                                apa += f", {pages}"
                        formatted['apa'] = apa
            
            citation['formatted'] = formatted
            return citation
        
        except Exception as e:
            logger.error(f"Error parsing CITATION.cff: {e}")
            return {}
    
    def extract_doi_from_text(self, text: str) -> Optional[str]:
        """Extract DOI from text content."""
        # DOI pattern
        doi_pattern = r'(?:doi:|DOI:|https?://doi.org/)(10\.\d+/[^\s\'"]+)'
        match = re.search(doi_pattern, text)
        if match:
            return match.group(1)
        return None