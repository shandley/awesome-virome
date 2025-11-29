#!/usr/bin/env python3
"""
Web of Science API client for retrieving citation data from Clarivate.
"""

import datetime
import json
import logging
import urllib.parse
import requests
from typing import Any, Dict, List, Optional, Tuple, Union

from ..base_client import BaseAPIClient

logger = logging.getLogger(__name__)


class WebOfScienceClient(BaseAPIClient):
    """Client for retrieving citation data from Web of Science."""
    
    def __init__(self, api_url: str, api_key: str, api_secret: str = "", rate_limit: float = 2.0):
        """
        Initialize the Web of Science API client.
        
        Args:
            api_url: Base URL for the Web of Science API
            api_key: Web of Science API key for authentication
            api_secret: Web of Science API secret for authentication
            rate_limit: Maximum requests per second
        """
        super().__init__(
            base_url=api_url,
            rate_limit=rate_limit,
            cache_namespace="wos"
        )
        
        self.api_key = api_key
        self.api_secret = api_secret
        self.auth_token = None
        self.auth_token_expiry = None
        
        # Track API credential status
        self.has_credentials = bool(api_key and api_secret)
    
    def _ensure_auth_token(self) -> bool:
        """
        Ensure that we have a valid authentication token.
        
        Returns:
            True if a valid token is available, False otherwise
        """
        # Check if we already have a valid token
        if self.auth_token and self.auth_token_expiry and datetime.datetime.now() < self.auth_token_expiry:
            return True
        
        # If no credentials, can't get a token
        if not self.has_credentials:
            logger.warning("Cannot authenticate with Web of Science API: missing credentials")
            return False
        
        # Request a new token
        try:
            auth_url = "https://api.clarivate.com/apis/wos-starter/auth/v1"
            
            headers = {
                "Accept": "application/json",
                "Content-Type": "application/json"
            }
            
            payload = {
                "client_id": self.api_key,
                "client_secret": self.api_secret
            }
            
            response = requests.post(auth_url, headers=headers, json=payload, timeout=30)
            response.raise_for_status()
            
            auth_data = response.json()
            
            # Extract token and expiry
            self.auth_token = auth_data.get("access_token")
            expires_in = auth_data.get("expires_in", 3600)  # Default to 1 hour if not specified
            
            # Calculate expiry time with a safety margin (subtract 5 minutes)
            self.auth_token_expiry = datetime.datetime.now() + datetime.timedelta(seconds=expires_in - 300)
            
            logger.info("Successfully authenticated with Web of Science API")
            return True
            
        except Exception as e:
            logger.error(f"Failed to authenticate with Web of Science API: {e}")
            self.auth_token = None
            self.auth_token_expiry = None
            return False
    
    def get_citation_data(self, doi: str, use_cache: bool = True) -> Tuple[bool, Dict[str, Any]]:
        """
        Get citation data for a given DOI from Web of Science.
        
        Args:
            doi: DOI to look up
            use_cache: Whether to use cached responses
        
        Returns:
            Tuple containing success status and citation data or error message
        """
        # Clean the DOI and encode for URL
        clean_doi = doi.strip().lower()
        
        # Skip Zenodo DOIs as they're generally not in Web of Science
        if clean_doi.startswith('10.5281/zenodo'):
            return False, f"DOI {doi} is a Zenodo DOI which is typically not indexed in Web of Science"
        
        # Ensure we have a valid authentication token
        if not self._ensure_auth_token():
            return False, "Could not authenticate with Web of Science API"
        
        # Search for the publication by DOI
        encoded_doi = urllib.parse.quote_plus(clean_doi)
        endpoint = f"document/doi/{encoded_doi}"
        
        headers = {
            "Accept": "application/json",
            "Authorization": f"Bearer {self.auth_token}"
        }
        
        success, response = self._make_request(
            endpoint=endpoint,
            headers=headers,
            use_cache=use_cache
        )
        
        if not success:
            logger.warning(f"Failed to get Web of Science data for DOI {doi}: {response}")
            return False, f"Failed to get Web of Science data: {response}"
        
        try:
            # Extract data from response
            data = response.get("data", [])
            if not data or len(data) == 0:
                return False, f"DOI {doi} not found in Web of Science"
            
            # Get the first matching document
            document = data[0]
            
            # Extract citation count
            citation_count = 0
            metrics = document.get("metrics", {})
            if "timesCited" in metrics:
                citation_count = int(metrics["timesCited"])
            
            # Extract metadata
            title = document.get("title", {}).get("value", "")
            
            # Extract authors
            authors = []
            author_list = document.get("authors", [])
            for author in author_list:
                name = author.get("name", "")
                if name:
                    authors.append(name)
            
            # Extract journal information
            journal = ""
            source = document.get("source", {})
            if source:
                journal = source.get("sourceTitle", "")
            
            # Extract publication year
            year = None
            pub_info = document.get("publicationDate", {})
            if pub_info:
                year_str = pub_info.get("year", "")
                if year_str:
                    try:
                        year = int(year_str)
                    except (ValueError, TypeError):
                        logger.warning(f"Invalid year value for DOI {doi}")
            
            # Get yearly citation breakdown if available
            citations_by_year = self._get_yearly_citations(doi, document.get("uid", ""), use_cache)
            
            # Create the citation data object
            citation_data = {
                "doi": doi,
                "source": "wos",
                "total_citations": citation_count,
                "citations_by_year": citations_by_year,
                "year": year,
                "timestamp": datetime.datetime.now().isoformat(),
                "metadata": {
                    "title": title,
                    "authors": authors,
                    "journal": journal,
                    "year": year,
                    "doi": doi,
                    "wos_id": document.get("uid", "")
                }
            }
            
            logger.info(f"Successfully retrieved Web of Science data for DOI: {doi} - Citations: {citation_count}")
            return True, citation_data
            
        except Exception as e:
            error_msg = f"Error parsing Web of Science data: {str(e)}"
            logger.error(error_msg)
            return False, error_msg
    
    def _get_yearly_citations(self, doi: str, wos_id: str, use_cache: bool = True) -> Dict[str, int]:
        """
        Get yearly citation breakdown for a publication.
        
        Args:
            doi: DOI of the publication
            wos_id: Web of Science ID (UID) of the publication
            use_cache: Whether to use cached responses
        
        Returns:
            Dictionary mapping years to citation counts
        """
        # If no WoS ID, we can't get yearly citations
        if not wos_id:
            return {}
        
        # Ensure we have a valid authentication token
        if not self._ensure_auth_token():
            logger.warning(f"Could not authenticate with Web of Science API to get yearly citations for {doi}")
            return {}
        
        endpoint = f"document/{wos_id}/citations/year"
        
        headers = {
            "Accept": "application/json",
            "Authorization": f"Bearer {self.auth_token}"
        }
        
        success, response = self._make_request(
            endpoint=endpoint,
            headers=headers,
            use_cache=use_cache
        )
        
        citations_by_year = {}
        
        if not success:
            logger.warning(f"Failed to get yearly citation breakdown for DOI {doi}: {response}")
            return citations_by_year
        
        try:
            # Extract yearly data from the response
            yearly_data = response.get("data", [])
            
            for year_data in yearly_data:
                year = year_data.get("year")
                count = year_data.get("count")
                
                if year and count and int(count) > 0:
                    citations_by_year[year] = int(count)
            
            return citations_by_year
            
        except Exception as e:
            logger.error(f"Error extracting yearly citation data for DOI {doi}: {e}")
            return citations_by_year
    
    def batch_get_citation_data(self, dois: List[str], use_cache: bool = True) -> Dict[str, Dict[str, Any]]:
        """
        Get citation data for multiple DOIs in batch.
        
        Args:
            dois: List of DOIs to look up
            use_cache: Whether to use cached responses
        
        Returns:
            Dictionary mapping DOIs to citation data
        """
        results = {}
        
        # Filter out Zenodo DOIs
        filtered_dois = [doi for doi in dois if not doi.strip().lower().startswith('10.5281/zenodo')]
        
        if len(filtered_dois) < len(dois):
            logger.info(f"Filtered out {len(dois) - len(filtered_dois)} Zenodo DOIs which are typically not indexed in Web of Science")
        
        # Ensure we have a valid authentication token for the batch
        if not self._ensure_auth_token() and filtered_dois:
            logger.warning("Could not authenticate with Web of Science API for batch request")
            error_msg = "Could not authenticate with Web of Science API"
            
            # Return error for all DOIs
            for doi in dois:
                results[doi] = {
                    "error": error_msg,
                    "source": "wos",
                    "doi": doi
                }
            
            return results
        
        # Process each DOI individually
        for doi in filtered_dois:
            success, data = self.get_citation_data(doi, use_cache)
            if success:
                results[doi] = data
            else:
                logger.warning(f"Failed to get Web of Science data for DOI {doi}: {data}")
                results[doi] = {"error": str(data), "source": "wos", "doi": doi}
        
        # Add skipped Zenodo DOIs to results
        for doi in dois:
            if doi not in results and doi.strip().lower().startswith('10.5281/zenodo'):
                results[doi] = {
                    "error": "Zenodo DOI not indexed in Web of Science",
                    "source": "wos",
                    "doi": doi
                }
        
        return results