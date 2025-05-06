#!/usr/bin/env python3
"""
Scopus API client for retrieving citation data from Elsevier Scopus.
"""

import datetime
import logging
import urllib.parse
import requests
from typing import Any, Dict, List, Optional, Tuple, Union

from ..base_client import BaseAPIClient

logger = logging.getLogger(__name__)


class ScopusClient(BaseAPIClient):
    """Client for retrieving citation data from Scopus."""
    
    def __init__(self, api_url: str, api_key: str, institutional_token: str = "", rate_limit: float = 5.0):
        """
        Initialize the Scopus API client.
        
        Args:
            api_url: Base URL for the Scopus API
            api_key: Scopus API key for authentication
            institutional_token: Institutional token for additional access (optional)
            rate_limit: Maximum requests per second
        """
        super().__init__(
            base_url=api_url,
            rate_limit=rate_limit,
            cache_namespace="scopus"
        )
        
        # Set up authentication headers
        self.headers = {
            "Accept": "application/json",
            "X-ELS-APIKey": api_key
        }
        
        # Add institutional token if provided
        if institutional_token:
            self.headers["X-ELS-Insttoken"] = institutional_token
    
    def get_citation_data(self, doi: str, use_cache: bool = True) -> Tuple[bool, Dict[str, Any]]:
        """
        Get citation data for a given DOI from Scopus.
        
        Args:
            doi: DOI to look up
            use_cache: Whether to use cached responses
        
        Returns:
            Tuple containing success status and citation data or error message
        """
        # Clean the DOI and encode for URL
        clean_doi = doi.strip().lower()
        encoded_doi = urllib.parse.quote_plus(clean_doi)
        
        # First, get the Scopus ID for this DOI
        scopus_id = self._get_scopus_id(clean_doi, use_cache)
        if not scopus_id:
            return False, f"DOI {doi} not found in Scopus"
        
        # Now get citation data using the Scopus ID
        return self._get_citations_by_scopus_id(scopus_id, doi, use_cache)
    
    def _get_scopus_id(self, doi: str, use_cache: bool = True) -> Optional[str]:
        """
        Get Scopus ID for a DOI.
        
        Args:
            doi: DOI to look up
            use_cache: Whether to use cached responses
        
        Returns:
            Scopus ID as string if found, None otherwise
        """
        encoded_doi = urllib.parse.quote_plus(doi)
        endpoint = f"content/search/scopus?query=DOI({encoded_doi})&view=STANDARD"
        
        success, response = self._make_request(
            endpoint=endpoint,
            headers=self.headers,
            use_cache=use_cache
        )
        
        if not success:
            logger.warning(f"Failed to get Scopus ID for DOI {doi}: {response}")
            return None
        
        # Extract Scopus ID from response
        try:
            results = response.get("search-results", {}).get("entry", [])
            if not results:
                logger.debug(f"No Scopus records found for DOI {doi}")
                return None
            
            # Get the first matching result
            return results[0].get("dc:identifier", "").replace("SCOPUS_ID:", "")
        except Exception as e:
            logger.error(f"Error extracting Scopus ID for DOI {doi}: {e}")
            return None
    
    def _get_citations_by_scopus_id(self, scopus_id: str, doi: str, use_cache: bool = True) -> Tuple[bool, Dict[str, Any]]:
        """
        Get citation data for a Scopus ID.
        
        Args:
            scopus_id: Scopus ID to look up
            doi: Original DOI (for reference)
            use_cache: Whether to use cached responses
        
        Returns:
            Tuple containing success status and citation data or error message
        """
        # Get the abstract and citation data
        endpoint = f"content/abstract/scopus_id/{scopus_id}"
        
        success, response = self._make_request(
            endpoint=endpoint,
            headers=self.headers,
            use_cache=use_cache
        )
        
        if not success:
            return False, f"Failed to get Scopus data: {response}"
        
        try:
            # Extract abstract data
            abstract_data = response.get("abstracts-retrieval-response", {})
            if not abstract_data:
                return False, f"No abstract data found for Scopus ID {scopus_id}"
            
            # Extract citation count - Scopus has detailed citation info
            citation_info = abstract_data.get("coredata", {})
            citation_count = int(citation_info.get("citedby-count", "0"))
            
            # Extract bibliographic data
            biblio_data = abstract_data.get("item", {}).get("bibrecord", {}).get("head", {})
            title = ""
            journal = ""
            authors = []
            year = None
            
            # Extract title
            source_data = biblio_data.get("citation", {}).get("source", {})
            if "sourcetitle" in source_data:
                journal = source_data["sourcetitle"]
            
            # Extract publication year
            citation_info = biblio_data.get("citation", {})
            if "publicationyear" in citation_info:
                try:
                    year = int(citation_info["publicationyear"].get("@first", "0"))
                except (ValueError, TypeError):
                    logger.warning(f"Invalid year value for Scopus ID {scopus_id}")
            
            # Extract title
            title_elements = biblio_data.get("citation", {}).get("titletext", {}).get("titletext", "")
            if title_elements:
                title = title_elements
            
            # Extract authors
            author_group = biblio_data.get("author-group", [])
            if isinstance(author_group, dict):
                author_group = [author_group]
            
            for group in author_group:
                author_list = group.get("author", [])
                if isinstance(author_list, dict):
                    author_list = [author_list]
                
                for author in author_list:
                    given_name = author.get("ce:given-name", "")
                    surname = author.get("ce:surname", "")
                    if given_name and surname:
                        authors.append(f"{surname}, {given_name}")
                    elif surname:
                        authors.append(surname)
            
            # Now get yearly citation breakdown if available
            citations_by_year = self._get_yearly_citations(scopus_id, use_cache)
            
            # Create the citation data object
            citation_data = {
                "doi": doi,
                "source": "scopus",
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
                    "scopus_id": scopus_id
                }
            }
            
            logger.info(f"Successfully retrieved Scopus data for DOI: {doi} - Citations: {citation_count}")
            return True, citation_data
            
        except Exception as e:
            error_msg = f"Error parsing Scopus data: {str(e)}"
            logger.error(error_msg)
            return False, error_msg
    
    def _get_yearly_citations(self, scopus_id: str, use_cache: bool = True) -> Dict[str, int]:
        """
        Get yearly citation breakdown for a Scopus ID.
        
        Args:
            scopus_id: Scopus ID to look up
            use_cache: Whether to use cached responses
        
        Returns:
            Dictionary mapping years to citation counts
        """
        # Get citing documents by year
        endpoint = f"content/abstract/citations?scopus_id={scopus_id}&date=2000-{datetime.datetime.now().year}&view=STANDARD"
        
        success, response = self._make_request(
            endpoint=endpoint,
            headers=self.headers,
            use_cache=use_cache
        )
        
        citations_by_year = {}
        
        if not success:
            logger.warning(f"Failed to get yearly citation breakdown for Scopus ID {scopus_id}")
            return citations_by_year
        
        try:
            # Extract yearly data from the response
            yearly_data = response.get("abstract-citations-response", {}).get("citeInfoMatrix", {}).get("citeInfoMatrixXML", {}).get("citationMatrix", {}).get("citeInfo", [])
            
            if not yearly_data:
                return citations_by_year
            
            # Scopus returns data as an array of yearly counts
            for year_data in yearly_data:
                year = year_data.get("year")
                count = year_data.get("valueByYear")
                
                if year and count and int(count) > 0:
                    citations_by_year[year] = int(count)
            
            return citations_by_year
            
        except Exception as e:
            logger.error(f"Error extracting yearly citation data for Scopus ID {scopus_id}: {e}")
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
        
        for doi in dois:
            success, data = self.get_citation_data(doi, use_cache)
            if success:
                results[doi] = data
            else:
                logger.warning(f"Failed to get Scopus data for DOI {doi}: {data}")
                results[doi] = {"error": str(data)}
        
        return results