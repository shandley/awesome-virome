#!/usr/bin/env python3
"""
iCite API Client for Awesome-Virome

This module provides a clean, focused client for the NIH iCite API to fetch
citation data for academic publications by DOI. It includes:

- Fetch citation data for a single DOI
- Batch processing of multiple DOIs
- Rate limiting and retry mechanisms
- Structured data validation

The NIH iCite API provides citation information for biomedical publications
and is the primary source of citation data for the Awesome-Virome project.
"""

import json
import logging
import time
import urllib.parse
from typing import Dict, List, Any, Optional, Union
import requests

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger("icite_api")

class ICiteClient:
    """Client for the NIH iCite API."""
    
    # Base URL for the iCite API
    BASE_URL = "https://icite.od.nih.gov/api"
    
    # Endpoints
    PMID_ENDPOINT = "/pubs"
    DOI_ENDPOINT = "/pubs/doi"
    
    def __init__(self, rate_limit: float = 3.0, max_retries: int = 3):
        """
        Initialize the iCite API client.
        
        Args:
            rate_limit: Minimum seconds between API requests
            max_retries: Maximum number of retry attempts for failed requests
        """
        self.rate_limit = rate_limit
        self.max_retries = max_retries
        self.last_request_time = 0
        
        logger.info(f"Initialized iCite API client (rate_limit={rate_limit}s, max_retries={max_retries})")
    
    def _rate_limit_request(self):
        """Apply rate limiting to API requests."""
        current_time = time.time()
        elapsed = current_time - self.last_request_time
        
        if elapsed < self.rate_limit:
            sleep_time = self.rate_limit - elapsed
            logger.debug(f"Rate limiting: sleeping for {sleep_time:.2f}s")
            time.sleep(sleep_time)
        
        self.last_request_time = time.time()
    
    def _make_request(self, url: str) -> Optional[Dict[str, Any]]:
        """
        Make a request to the iCite API with retry logic.
        
        Args:
            url: The full URL to request
            
        Returns:
            Response data as a dictionary, or None if the request failed
        """
        self._rate_limit_request()
        
        for attempt in range(1, self.max_retries + 1):
            try:
                logger.debug(f"Making request to {url} (attempt {attempt}/{self.max_retries})")
                response = requests.get(url, timeout=30)
                
                # Check for successful response
                if response.status_code == 200:
                    return response.json()
                
                # Handle common errors
                if response.status_code == 404:
                    logger.warning(f"Resource not found: {url}")
                    return None
                elif response.status_code == 429:
                    wait_time = min(2 ** attempt, 60)  # Exponential backoff
                    logger.warning(f"Rate limit exceeded, waiting {wait_time}s before retry")
                    time.sleep(wait_time)
                else:
                    logger.error(f"HTTP error {response.status_code} for {url}")
                    
                    # Only retry on potentially temporary errors
                    if response.status_code >= 500 or response.status_code == 429:
                        continue
                    return None
                    
            except (requests.RequestException, json.JSONDecodeError) as e:
                logger.error(f"Request error: {str(e)}")
                
                # Wait before retry
                if attempt < self.max_retries:
                    wait_time = 2 ** attempt
                    logger.info(f"Retrying in {wait_time}s...")
                    time.sleep(wait_time)
        
        logger.error(f"Failed to retrieve data from {url} after {self.max_retries} attempts")
        return None
    
    def get_citation_data_by_doi(self, doi: str) -> Optional[Dict[str, Any]]:
        """
        Fetch citation data for a publication by DOI.
        
        Args:
            doi: Digital Object Identifier for the publication
            
        Returns:
            Dictionary containing citation data, or None if not found
        """
        # Clean and encode the DOI
        doi = doi.strip()
        encoded_doi = urllib.parse.quote(doi, safe='')
        
        # Build the request URL
        url = f"{self.BASE_URL}{self.DOI_ENDPOINT}/{encoded_doi}"
        
        logger.info(f"Fetching citation data for DOI: {doi}")
        response_data = self._make_request(url)
        
        if not response_data:
            logger.warning(f"No citation data found for DOI: {doi}")
            return None
        
        # Extract the actual record from the response
        if "data" in response_data and len(response_data["data"]) > 0:
            citation_data = response_data["data"][0]
            logger.info(f"Successfully retrieved citation data for DOI: {doi}")
            
            # Add source attribution
            citation_data["citation_source"] = "icite"
            
            return self._process_citation_data(citation_data)
        
        logger.warning(f"Empty or invalid response for DOI: {doi}")
        return None
    
    def get_citation_data_by_pmid(self, pmid: str) -> Optional[Dict[str, Any]]:
        """
        Fetch citation data for a publication by PubMed ID.
        
        Args:
            pmid: PubMed ID for the publication
            
        Returns:
            Dictionary containing citation data, or None if not found
        """
        # Clean and validate the PMID
        pmid = pmid.strip()
        if not pmid.isdigit():
            logger.warning(f"Invalid PMID format: {pmid}")
            return None
        
        # Build the request URL
        url = f"{self.BASE_URL}{self.PMID_ENDPOINT}/{pmid}"
        
        logger.info(f"Fetching citation data for PMID: {pmid}")
        response_data = self._make_request(url)
        
        if not response_data:
            logger.warning(f"No citation data found for PMID: {pmid}")
            return None
        
        # Extract the data
        if "data" in response_data:
            citation_data = response_data["data"]
            logger.info(f"Successfully retrieved citation data for PMID: {pmid}")
            
            # Add source attribution
            citation_data["citation_source"] = "icite"
            
            return self._process_citation_data(citation_data)
        
        logger.warning(f"Empty or invalid response for PMID: {pmid}")
        return None
    
    def _process_citation_data(self, data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Process and normalize raw citation data from iCite.
        
        Args:
            data: Raw citation data from iCite API
            
        Returns:
            Processed citation data with consistent structure
        """
        # Extract the core citation metrics
        processed = {
            "doi": data.get("doi", ""),
            "pmid": str(data.get("pmid", "")),
            "year": data.get("year", 0),
            "title": data.get("title", ""),
            "authors": data.get("authors", ""),
            "journal": data.get("journal", ""),
            "citation_count": data.get("citation_count", 0),
            "relative_citation_ratio": data.get("relative_citation_ratio", 0),
            "field_citation_rate": data.get("field_citation_rate", 0),
            "citation_source": "icite"
        }
        
        # Extract citations by year if available
        citations_by_year = {}
        if "citations_by_year" in data and isinstance(data["citations_by_year"], list):
            for year_data in data["citations_by_year"]:
                year = str(year_data.get("year", ""))
                count = year_data.get("citations", 0)
                if year and count > 0:
                    citations_by_year[year] = count
        
        processed["citations_by_year"] = citations_by_year
        
        return processed
    
    def batch_get_citation_data(self, dois: List[str], chunk_size: int = 10) -> Dict[str, Dict[str, Any]]:
        """
        Fetch citation data for multiple DOIs in batches.
        
        Args:
            dois: List of DOIs to fetch citation data for
            chunk_size: Number of DOIs to process in each batch
            
        Returns:
            Dictionary mapping DOIs to their citation data
        """
        results = {}
        total_dois = len(dois)
        
        logger.info(f"Batch processing {total_dois} DOIs in chunks of {chunk_size}")
        
        # Process DOIs in chunks to avoid overwhelming the API
        for i in range(0, total_dois, chunk_size):
            chunk = dois[i:i+chunk_size]
            logger.info(f"Processing chunk {i//chunk_size + 1}/{(total_dois+chunk_size-1)//chunk_size} ({len(chunk)} DOIs)")
            
            for doi in chunk:
                citation_data = self.get_citation_data_by_doi(doi)
                if citation_data:
                    results[doi] = citation_data
            
            # Log progress
            processed = i + len(chunk)
            logger.info(f"Processed {processed}/{total_dois} DOIs ({len(results)} with citation data)")
        
        return results


def test_icite_client():
    """Simple test function for the iCite client."""
    client = ICiteClient()
    
    # Test with a known DOI
    test_doi = "10.1093/nar/gkaa939"  # CheckV paper
    
    logger.info(f"Testing iCite client with DOI: {test_doi}")
    data = client.get_citation_data_by_doi(test_doi)
    
    if data:
        logger.info(f"Success! Found {data['citation_count']} citations for '{data['title']}'")
        if 'citations_by_year' in data:
            logger.info(f"Citations by year: {data['citations_by_year']}")
        return True
    else:
        logger.error("Test failed - no data returned")
        return False


if __name__ == "__main__":
    # Run a simple test if executed directly
    test_icite_client()