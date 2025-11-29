#!/usr/bin/env python3
"""
CrossRef API client for retrieving bibliographic data.
"""

import datetime
import json
import logging
import urllib.parse
from typing import Any, Dict, List, Optional, Tuple, Union

from ..base_client import BaseAPIClient

logger = logging.getLogger(__name__)


class CrossRefClient(BaseAPIClient):
    """Client for retrieving bibliographic data from CrossRef."""
    
    def __init__(self, api_url: str, rate_limit: float):
        """
        Initialize the CrossRef API client.
        
        Args:
            api_url: Base URL for the CrossRef API
            rate_limit: Maximum requests per second
        """
        super().__init__(
            base_url=api_url,
            rate_limit=rate_limit,
            cache_namespace="crossref"
        )
        
        # User agent headers are required by CrossRef
        self.headers = {
            "User-Agent": "AwesomeVirome/1.0 (https://github.com/shandley/awesome-virome; mailto:your-email@example.com)",
            "Accept": "application/json"
        }
    
    def get_citation_data(self, doi: str, use_cache: bool = True) -> Tuple[bool, Dict[str, Any]]:
        """
        Get citation data for a given DOI from CrossRef.
        
        Note: CrossRef does not directly provide citation counts.
        This method returns basic metadata that can be used for validation.
        
        Args:
            doi: DOI to look up
            use_cache: Whether to use cached responses
        
        Returns:
            Tuple containing success status and bibliographic data or error message
        """
        # Clean the DOI and encode for URL
        clean_doi = doi.strip().lower()
        encoded_doi = urllib.parse.quote_plus(clean_doi)
        
        # Construct the CrossRef API endpoint
        endpoint = f"works/{encoded_doi}"
        
        # Make the API request
        success, response = self._make_request(
            endpoint=endpoint,
            headers=self.headers,
            use_cache=use_cache
        )
        
        if not success:
            return False, f"Failed to get CrossRef data: {response}"
        
        # CrossRef returns data in a specific format
        if "message" not in response:
            return False, f"Invalid CrossRef response for DOI {doi}"
        
        message = response["message"]
        
        try:
            # CrossRef doesn't provide citation counts directly
            # We're just getting metadata
            
            # Extract basic information
            title = message.get("title", [""])[0] if message.get("title") else ""
            
            # Extract authors
            authors = []
            if "author" in message:
                for author in message["author"]:
                    author_name = []
                    if "given" in author:
                        author_name.append(author["given"])
                    if "family" in author:
                        author_name.append(author["family"])
                    
                    if author_name:
                        authors.append(" ".join(author_name))
            
            # Extract publication date
            year = None
            if "published" in message and "date-parts" in message["published"]:
                date_parts = message["published"]["date-parts"]
                if date_parts and len(date_parts[0]) > 0:
                    year = date_parts[0][0]
            
            # Extract journal/publisher information
            journal = message.get("container-title", [""])[0] if message.get("container-title") else ""
            publisher = message.get("publisher", "")
            
            # Create bibliographic data (without citation count)
            bibliographic_data = {
                "doi": doi,
                "source": "crossref",
                "total_citations": 0,  # CrossRef doesn't provide this directly
                "citations_by_year": {},
                "timestamp": datetime.datetime.now().isoformat(),
                "metadata": {
                    "title": title,
                    "authors": authors,
                    "journal": journal,
                    "publisher": publisher,
                    "year": year,
                    "doi": message.get("DOI", doi)
                }
            }
            
            logger.info(f"Successfully retrieved CrossRef data for DOI: {doi}")
            return True, bibliographic_data
            
        except Exception as e:
            error_msg = f"Error parsing CrossRef data: {str(e)}"
            logger.error(error_msg)
            return False, error_msg
    
    def search_works(self, query: str, limit: int = 10, use_cache: bool = True) -> Tuple[bool, Dict[str, Any]]:
        """
        Search for works in CrossRef.
        
        Args:
            query: Search query
            limit: Maximum number of results
            use_cache: Whether to use cached responses
        
        Returns:
            Tuple containing success status and search results or error message
        """
        # Construct the search query
        encoded_query = urllib.parse.quote_plus(query)
        
        # Construct the CrossRef API endpoint with parameters
        params = {
            "query": query,
            "rows": limit,
            "sort": "relevance"
        }
        
        # Make the API request
        success, response = self._make_request(
            endpoint="works",
            params=params,
            headers=self.headers,
            use_cache=use_cache
        )
        
        if not success:
            return False, f"Failed to search CrossRef: {response}"
        
        # Check for valid response
        if "message" not in response or "items" not in response["message"]:
            return False, "Invalid CrossRef search response"
        
        # Extract search results
        items = response["message"]["items"]
        
        # Process and return results
        results = []
        for item in items:
            try:
                # Extract basic information
                doi = item.get("DOI", "")
                title = item.get("title", [""])[0] if item.get("title") else ""
                
                # Extract authors
                authors = []
                if "author" in item:
                    for author in item["author"]:
                        author_name = []
                        if "given" in author:
                            author_name.append(author["given"])
                        if "family" in author:
                            author_name.append(author["family"])
                        
                        if author_name:
                            authors.append(" ".join(author_name))
                
                # Extract publication date
                year = None
                if "published" in item and "date-parts" in item["published"]:
                    date_parts = item["published"]["date-parts"]
                    if date_parts and len(date_parts[0]) > 0:
                        year = date_parts[0][0]
                
                # Extract journal/publisher information
                journal = item.get("container-title", [""])[0] if item.get("container-title") else ""
                publisher = item.get("publisher", "")
                
                # Add to results
                results.append({
                    "doi": doi,
                    "title": title,
                    "authors": authors,
                    "journal": journal,
                    "publisher": publisher,
                    "year": year
                })
                
            except Exception as e:
                logger.warning(f"Error processing search result: {e}")
        
        return True, {
            "source": "crossref",
            "query": query,
            "total_results": response["message"].get("total-results", 0),
            "items": results
        }