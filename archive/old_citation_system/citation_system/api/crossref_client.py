#!/usr/bin/env python3
"""
CrossRef API client for retrieving citation data.
"""

import datetime
import logging
import re
from typing import Any, Dict, List, Optional, Tuple, Union

from ..config import CROSSREF_API_URL, CROSSREF_RATE_LIMIT
from .base_client import BaseAPIClient

logger = logging.getLogger(__name__)


class CrossRefClient(BaseAPIClient):
    """Client for retrieving citation data from CrossRef."""
    
    def __init__(self):
        """Initialize the CrossRef API client."""
        super().__init__(
            base_url=CROSSREF_API_URL,
            rate_limit=CROSSREF_RATE_LIMIT,
            cache_namespace="crossref"
        )
        
        # Set up headers with proper User-Agent
        self.headers = {
            "User-Agent": "AwesomeVirome/1.0 (https://github.com/shandley/awesome-virome; mailto:your.email@example.com)",
            "Accept": "application/json"
        }
    
    def get_publication_data(
        self,
        doi: str,
        use_cache: bool = True
    ) -> Tuple[bool, Union[Dict[str, Any], str]]:
        """
        Get publication metadata from CrossRef.
        
        Args:
            doi: DOI to look up
            use_cache: Whether to use cached data
        
        Returns:
            Tuple containing success status and publication data or error message
        """
        # Clean up DOI - remove any URL prefix
        doi = doi.strip().lower()
        if doi.startswith('https://doi.org/'):
            doi = doi[len('https://doi.org/'):]
        
        logger.info(f"Fetching publication data for DOI: {doi}")
        
        # Make API request
        success, response = self._make_request(
            endpoint=doi,
            headers=self.headers,
            use_cache=use_cache
        )
        
        if not success:
            return False, f"Failed to get publication data: {response}"
        
        # Check for message
        if 'message' not in response:
            return False, "Invalid response format from CrossRef"
        
        return True, response['message']
    
    def get_citation_count(
        self,
        doi: str,
        use_cache: bool = True
    ) -> Tuple[bool, Union[Dict[str, Any], str]]:
        """
        Get citation count for a DOI.
        
        Args:
            doi: DOI to look up
            use_cache: Whether to use cached data
        
        Returns:
            Tuple containing success status and citation data or error message
        """
        # Get publication metadata
        success, pub_data = self.get_publication_data(doi, use_cache)
        
        if not success:
            return False, pub_data
        
        # CrossRef doesn't directly provide citation counts, but we can get
        # some useful metadata for our citation system
        
        # Extract relevant information with improved error handling
        try:
            # Safely extract title
            title = ''
            if 'title' in pub_data and isinstance(pub_data['title'], list) and pub_data['title']:
                title = pub_data['title'][0]
            
            # Get publication date
            published_date = None
            if 'published' in pub_data and 'date-parts' in pub_data['published']:
                date_parts = pub_data['published']['date-parts']
                if date_parts and isinstance(date_parts, list) and date_parts[0]:
                    date_part = date_parts[0]
                    if len(date_part) >= 3:
                        try:
                            published_date = datetime.date(date_part[0], date_part[1], date_part[2]).isoformat()
                        except (ValueError, TypeError):
                            published_date = str(date_part[0]) if date_part else None
                    elif len(date_part) >= 1:
                        published_date = str(date_part[0])
            
            # Get authors
            authors = []
            if 'author' in pub_data and isinstance(pub_data['author'], list):
                for author in pub_data['author']:
                    if isinstance(author, dict):
                        if 'family' in author and 'given' in author:
                            authors.append(f"{author['family']}, {author['given']}")
                        elif 'family' in author:
                            authors.append(author['family'])
            
            # Get journal/container info
            journal = ''
            if 'container-title' in pub_data and isinstance(pub_data['container-title'], list) and pub_data['container-title']:
                journal = pub_data['container-title'][0]
            
            # Structure the citation data
            citation_data = {
                'doi': doi,
                'source': 'crossref',
                'total_citations': 0,  # CrossRef doesn't provide this directly
                'timestamp': datetime.datetime.now().isoformat(),
                'metadata': {
                    'title': title,
                    'authors': authors,
                    'journal': journal,
                    'publication_date': published_date,
                    'type': pub_data.get('type', ''),
                    'publisher': pub_data.get('publisher', '')
                }
            }
            
            logger.info(f"Successfully retrieved metadata for DOI: {doi}")
            return True, citation_data
            
        except Exception as e:
            error_msg = f"Error parsing CrossRef data: {str(e)}"
            logger.error(error_msg)
            return False, error_msg
    
    def search_by_title(
        self,
        title: str,
        limit: int = 5,
        use_cache: bool = True
    ) -> Tuple[bool, Union[List[Dict[str, Any]], str]]:
        """
        Search for publications by title.
        
        Args:
            title: Publication title to search for
            limit: Maximum number of results to return
            use_cache: Whether to use cached data
        
        Returns:
            Tuple containing success status and list of publications or error message
        """
        logger.info(f"Searching for publication with title: {title}")
        
        # Make API request to the works endpoint with query
        try:
            success, response = self._make_request(
                endpoint="works",  # Updated endpoint to correct path
                params={
                    "query.title": title,
                    "rows": limit,
                    "sort": "relevance"
                },
                headers=self.headers,
                use_cache=use_cache
            )
            
            if not success:
                return False, f"Failed to search for title: {response}"
            
            # Check for message and items
            if 'message' not in response or 'items' not in response['message']:
                return False, "Invalid response format from CrossRef search"
            
            # Extract and return results
            results = []
            for item in response['message']['items']:
                # Only include items with DOIs
                if 'DOI' in item:
                    title_text = ''
                    if 'title' in item and isinstance(item['title'], list) and item['title']:
                        title_text = item['title'][0]
                    
                    results.append({
                        'doi': item['DOI'],
                        'title': title_text,
                        'score': item.get('score', 0)
                    })
            
            return True, results
        except Exception as e:
            error_msg = f"Error in search_by_title: {str(e)}"
            logger.error(error_msg)
            return False, error_msg