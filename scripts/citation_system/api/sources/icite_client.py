#!/usr/bin/env python3
"""
iCite API client for retrieving NIH citation data.
"""

import datetime
import logging
import time
import urllib.parse
import requests
from typing import Any, Dict, List, Optional, Tuple, Union

from ..base_client import BaseAPIClient

logger = logging.getLogger(__name__)


class ICiteClient(BaseAPIClient):
    """Client for retrieving citation data from iCite."""
    
    def __init__(self, api_url: str, rate_limit: float, pubmed_api_url: str = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"):
        """
        Initialize the iCite API client.
        
        Args:
            api_url: Base URL for the iCite API
            rate_limit: Maximum requests per second
            pubmed_api_url: Base URL for the PubMed API (used for DOI to PMID conversion)
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
        
        self.pubmed_api_url = pubmed_api_url.rstrip('/')
    
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
        
        # First, try to find the PMID for this DOI using PubMed API
        pmid = self._get_pmid_for_doi(clean_doi)
        if pmid:
            return self.search_by_pmid(pmid, use_cache)
        
        # Fallback to directly querying with DOI
        # The current API doesn't support direct DOI lookup, so we'll use the pmids endpoint
        # and filter by DOI in the results

        # Construct the iCite API endpoint to search for recent papers (as a fallback)
        endpoint = f"pubs?format=json&limit=1&doi={urllib.parse.quote_plus(clean_doi)}"
        
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
            citation_count_val = paper_data.get("citation_count")
            citation_count = int(citation_count_val) if citation_count_val is not None else 0
            
            # Extract relative citation ratio (RCR) - an NIH-specific metric
            rcr_val = paper_data.get("relative_citation_ratio")
            rcr = float(rcr_val) if rcr_val is not None else 0
            
            # Extract other metrics
            expected_val = paper_data.get("expected_citations_per_year")
            expected_citations = float(expected_val) if expected_val is not None else 0
            field_rate_val = paper_data.get("field_citation_rate")
            field_citation_rate = float(field_rate_val) if field_rate_val is not None else 0
            
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
    
    def _get_pmid_for_doi(self, doi: str) -> Optional[str]:
        """
        Convert a DOI to a PubMed ID (PMID) using the PubMed API.
        
        Args:
            doi: DOI to convert
            
        Returns:
            PMID as string if found, None otherwise
        """
        # Skip Zenodo DOIs as they're generally not in PubMed
        if doi.startswith('10.5281/zenodo'):
            return None
            
        try:
            encoded_doi = urllib.parse.quote_plus(doi)
            url = f"{self.pubmed_api_url}/esearch.fcgi?db=pubmed&term={encoded_doi}[DOI]&retmode=json&tool=awesome-virome&email=email@example.com"
            
            # Add delay to respect rate limits
            time.sleep(0.5)  # Wait 500ms between requests
            
            response = requests.get(url, timeout=10)
            response.raise_for_status()
            
            data = response.json()
            id_list = data.get('esearchresult', {}).get('idlist', [])
            
            if id_list and len(id_list) > 0:
                return id_list[0]
                
        except requests.exceptions.HTTPError as e:
            if '429' in str(e):
                # If rate limited, just skip this DOI lookup
                logger.debug(f"Rate limited by PubMed API for DOI: {doi}")
            else:
                logger.warning(f"Failed to convert DOI to PMID: {str(e)}")
        except Exception as e:
            logger.warning(f"Failed to convert DOI to PMID: {str(e)}")
            
        return None
            
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
        endpoint = f"pubs?pmids={pmid}"
        
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
        
        try:
            # Extract citation count
            citation_count_val = paper_data.get("citation_count")
            citation_count = int(citation_count_val) if citation_count_val is not None else 0
            
            # Extract relative citation ratio (RCR) - an NIH-specific metric
            rcr_val = paper_data.get("relative_citation_ratio")
            rcr = float(rcr_val) if rcr_val is not None else 0
            
            # Extract other metrics
            expected_val = paper_data.get("expected_citations_per_year")
            expected_citations = float(expected_val) if expected_val is not None else 0
            field_rate_val = paper_data.get("field_citation_rate")
            field_citation_rate = float(field_rate_val) if field_rate_val is not None else 0
            
            # Extract year info
            year = paper_data.get("year", None)
            
            # Create citations by year (iCite doesn't provide this directly)
            citations_by_year = {}
            if year and citation_count > 0:
                # We don't have actual yearly breakdown
                citations_by_year = {}
            
            # Extract title and authors
            title = paper_data.get("title", "")
            authors = []
            if "authors" in paper_data:
                authors = paper_data["authors"].split(", ")
            
            # Extract journal information
            journal = paper_data.get("journal", "")
            
            # Get the DOI from the paper data (if available)
            doi = paper_data.get("doi", "")
            
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
                    "doi": doi,
                    "year": year
                }
            }
            
            logger.info(f"Successfully retrieved iCite data for PMID: {pmid} - Citations: {citation_count}")
            return True, citation_data
            
        except Exception as e:
            error_msg = f"Error parsing iCite data: {str(e)}"
            logger.error(error_msg)
            return False, error_msg
    
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