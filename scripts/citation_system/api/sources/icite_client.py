#!/usr/bin/env python3
"""
iCite API client for retrieving NIH citation data.
"""

import datetime
import logging
import urllib.parse
from typing import Any, Dict, List, Optional, Tuple, Union

from ..base_client import BaseAPIClient

logger = logging.getLogger(__name__)


class ICiteClient(BaseAPIClient):
    """Client for retrieving citation data from iCite."""
    
    def __init__(self, api_url: str, rate_limit: float):
        """
        Initialize the iCite API client.
        
        Args:
            api_url: Base URL for the iCite API
            rate_limit: Maximum requests per second
        """
        super().__init__(
            base_url=api_url,
            rate_limit=rate_limit,
            cache_namespace="icite"
        )
        
        # iCite is open access, so no auth is needed
        self.headers = {
            "Accept": "application/json"
        }
    
    def get_citation_data(self, doi: str, use_cache: bool = True) -> Tuple[bool, Dict[str, Any]]:
        """
        Get citation data for a given DOI from iCite.
        
        Args:
            doi: DOI to look up
            use_cache: Whether to use cached responses
        
        Returns:
            Tuple containing success status and citation data or error message
        """
        # Clean the DOI and encode for URL
        clean_doi = doi.strip().lower()
        encoded_doi = urllib.parse.quote_plus(clean_doi)
        
        # Construct the iCite API endpoint
        endpoint = f"pubs/v1/doi/{encoded_doi}"
        
        # Make the API request
        success, response = self._make_request(
            endpoint=endpoint,
            headers=self.headers,
            use_cache=use_cache
        )
        
        if not success:
            return False, f"Failed to get iCite data: {response}"
        
        # iCite returns data in a specific format
        if "data" not in response or not response["data"]:
            return False, f"DOI {doi} not found in iCite"
        
        # Extract the paper data
        paper_data = response["data"][0]
        
        try:
            # Extract citation count
            citation_count = int(paper_data.get("citation_count", 0))
            
            # Extract relative citation ratio (RCR) - an NIH-specific metric
            rcr = float(paper_data.get("relative_citation_ratio", 0))
            
            # Extract other metrics
            expected_citations = float(paper_data.get("expected_citations_per_year", 0))
            field_citation_rate = float(paper_data.get("field_citation_rate", 0))
            
            # Extract year info
            year = paper_data.get("year", None)
            
            # Create citations by year (iCite doesn't provide this directly)
            citations_by_year = {}
            if year and citation_count > 0:
                # We don't have actual yearly breakdown
                # This is just a placeholder for the aggregation system
                citations_by_year = {}
            
            # Extract title and authors
            title = paper_data.get("title", "")
            authors = []
            if "authors" in paper_data:
                authors = paper_data["authors"].split("|")
            
            # Extract journal information
            journal = paper_data.get("journal", "")
            
            # Create the citation data object
            citation_data = {
                "doi": doi,
                "source": "icite",
                "total_citations": citation_count,
                "rcr": rcr,  # Relative Citation Ratio - unique to iCite
                "expected_citations": expected_citations,
                "field_citation_rate": field_citation_rate,
                "citations_by_year": citations_by_year,
                "is_clinical": paper_data.get("is_clinical", False),
                "year": year,
                "timestamp": datetime.datetime.now().isoformat(),
                "metadata": {
                    "title": title,
                    "authors": authors,
                    "journal": journal,
                    "pmid": paper_data.get("pmid", ""),
                    "doi": paper_data.get("doi", ""),
                    "year": year
                }
            }
            
            logger.info(f"Successfully retrieved iCite data for DOI: {doi} - Citations: {citation_count}")
            return True, citation_data
            
        except Exception as e:
            error_msg = f"Error parsing iCite data: {str(e)}"
            logger.error(error_msg)
            return False, error_msg
    
    def search_by_pmid(self, pmid: str, use_cache: bool = True) -> Tuple[bool, Dict[str, Any]]:
        """
        Get citation data for a given PubMed ID from iCite.
        
        Args:
            pmid: PubMed ID to look up
            use_cache: Whether to use cached responses
        
        Returns:
            Tuple containing success status and citation data or error message
        """
        # Construct the iCite API endpoint
        endpoint = f"pubs/v1/pmid/{pmid}"
        
        # Make the API request
        success, response = self._make_request(
            endpoint=endpoint,
            headers=self.headers,
            use_cache=use_cache
        )
        
        if not success:
            return False, f"Failed to get iCite data: {response}"
        
        # iCite returns data in a specific format
        if "data" not in response or not response["data"]:
            return False, f"PMID {pmid} not found in iCite"
        
        # Extract the paper data
        paper_data = response["data"][0]
        
        # Get the DOI from the paper data
        doi = paper_data.get("doi", "")
        
        # If DOI is available, use it to get citation data
        if doi:
            return self.get_citation_data(doi, use_cache)
        else:
            return False, f"No DOI found for PMID {pmid}"
    
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
                logger.warning(f"Failed to get iCite data for DOI {doi}: {data}")
                results[doi] = {"error": str(data)}
        
        return results