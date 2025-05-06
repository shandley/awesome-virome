#!/usr/bin/env python3
"""
Integration tests for the Web of Science citation source.

These tests verify that the Web of Science client can retrieve citation data
and that it integrates properly with the citation aggregator.

Note: These tests require valid Web of Science API credentials to run.
"""

import os
import unittest
import logging
from pathlib import Path
import sys

# Add parent directory to path to allow imports
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent.parent))

from scripts.citation_system.api.sources.wos_client import WebOfScienceClient
from scripts.citation_system.api.citation_registry import get_citation_source, get_available_sources
from scripts.citation_system.collectors.citation_aggregator import CitationAggregator
from scripts.citation_system.config import WOS_API_URL, WOS_API_KEY, WOS_API_SECRET, WOS_ENABLED

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger("test_wos_integration")


class TestWebOfScienceIntegration(unittest.TestCase):
    """Integration tests for Web of Science citation source."""

    @classmethod
    def setUpClass(cls):
        """Set up test fixtures that are reused across test methods."""
        # Check if Web of Science is enabled and API credentials are available
        cls.api_key = os.environ.get("WOS_API_KEY", WOS_API_KEY)
        cls.api_secret = os.environ.get("WOS_API_SECRET", WOS_API_SECRET)
        cls.wos_enabled = WOS_ENABLED and cls.api_key and cls.api_secret
        
        if not cls.wos_enabled:
            logger.warning("Web of Science tests skipped: WOS_ENABLED=false or no API credentials available")
            return
        
        # Initialize the Web of Science client
        cls.client = WebOfScienceClient(
            api_url=WOS_API_URL,
            api_key=cls.api_key,
            api_secret=cls.api_secret,
            rate_limit=2.0
        )
        
        # Test DOIs - these should be stable, well-known papers
        cls.test_dois = [
            "10.1038/s41586-020-2008-3",  # COVID-19 paper (high citations)
            "10.1093/bioinformatics/btab213"  # Bioinformatics paper
        ]
    
    def setUp(self):
        """Set up test fixtures."""
        if not self.wos_enabled:
            self.skipTest("Web of Science tests skipped: WOS_ENABLED=false or no API credentials available")
    
    def test_wos_client_direct(self):
        """Test direct use of Web of Science client."""
        for doi in self.test_dois:
            success, data = self.client.get_citation_data(doi)
            
            self.assertTrue(success, f"Failed to get Web of Science data for DOI {doi}")
            self.assertIn("total_citations", data, "Citation count missing from response")
            self.assertIn("metadata", data, "Metadata missing from response")
            self.assertIn("source", data, "Source identifier missing from response")
            self.assertEqual(data["source"], "wos", "Incorrect source in response")
            
            # Check that citation count is reasonable (greater than zero for these known papers)
            self.assertGreater(data["total_citations"], 0, f"Expected citations > 0 for {doi}")
            
            logger.info(f"Web of Science citation count for {doi}: {data['total_citations']}")
    
    def test_wos_client_batch(self):
        """Test batch citation retrieval."""
        results = self.client.batch_get_citation_data(self.test_dois)
        
        self.assertEqual(len(results), len(self.test_dois), "Not all DOIs processed")
        
        for doi in self.test_dois:
            self.assertIn(doi, results, f"DOI {doi} missing from batch results")
            data = results[doi]
            self.assertNotIn("error", data, f"Error in batch result for {doi}: {data.get('error', '')}")
            self.assertIn("total_citations", data, "Citation count missing from response")
            
            logger.info(f"Batch Web of Science citation count for {doi}: {data['total_citations']}")
    
    def test_citation_registry_integration(self):
        """Test that Web of Science is properly registered in the citation registry."""
        # Check if Web of Science is in available sources
        sources = get_available_sources()
        self.assertIn("wos", sources, "Web of Science not registered in citation sources")
        
        # Get Web of Science client instance from registry
        wos = get_citation_source("wos")
        self.assertIsNotNone(wos, "Failed to get Web of Science client from registry")
        self.assertIsInstance(wos, WebOfScienceClient, "Registry returned wrong client type")
    
    def test_citation_aggregator_with_wos(self):
        """Test citation aggregation with Web of Science as a source."""
        aggregator = CitationAggregator()
        
        for doi in self.test_dois:
            aggregated_data = aggregator.collect_citations(doi)
            
            self.assertNotIn("error", aggregated_data, 
                            f"Error in aggregated data for {doi}: {aggregated_data.get('error', '')}")
            self.assertIn("total_citations", aggregated_data, "Citation count missing from aggregated data")
            self.assertIn("sources_used", aggregated_data, "Sources used missing from aggregated data")
            self.assertIn("wos", aggregated_data["sources_used"], 
                         "Web of Science not listed in sources used for aggregation")
            
            # Check if Web of Science was selected as a data source
            # It might not be the primary source if Scopus has higher priority
            if aggregated_data.get("primary_source") == "wos":
                logger.info(f"Web of Science selected as primary source for {doi}")
            else:
                logger.info(f"Primary source for {doi}: {aggregated_data.get('primary_source', 'unknown')}")
            
            logger.info(f"Aggregated citation count for {doi}: {aggregated_data['total_citations']}")
            
            # Check yearly citation data source
            if "yearly_data_source" in aggregated_data:
                logger.info(f"Yearly citation data source for {doi}: {aggregated_data['yearly_data_source']}")
                if aggregated_data["yearly_data_source"] == "wos":
                    self.assertIn("citations_by_year", aggregated_data)
                    self.assertTrue(len(aggregated_data["citations_by_year"]) > 0)
    
    def test_zenodo_doi_handling(self):
        """Test handling of Zenodo DOIs."""
        # Zenodo DOIs should be skipped by Web of Science
        zenodo_doi = "10.5281/zenodo.1234567"
        success, data = self.client.get_citation_data(zenodo_doi)
        
        self.assertFalse(success)
        self.assertTrue("Zenodo DOI" in data)
        
        # Test in batch
        results = self.client.batch_get_citation_data([zenodo_doi])
        self.assertIn(zenodo_doi, results)
        self.assertIn("error", results[zenodo_doi])
        self.assertTrue("Zenodo DOI" in results[zenodo_doi]["error"])


if __name__ == "__main__":
    unittest.main()