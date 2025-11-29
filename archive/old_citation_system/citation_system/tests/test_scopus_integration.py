#!/usr/bin/env python3
"""
Integration tests for the Scopus citation source.

These tests verify that the Scopus client can retrieve citation data
and that it integrates properly with the citation aggregator.

Note: These tests require a valid Scopus API key to run.
"""

import os
import unittest
import logging
from pathlib import Path
import sys

# Add parent directory to path to allow imports
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent.parent))

from scripts.citation_system.api.sources.scopus_client import ScopusClient
from scripts.citation_system.api.citation_registry import get_citation_source, get_available_sources
from scripts.citation_system.collectors.citation_aggregator import CitationAggregator
from scripts.citation_system.config import SCOPUS_API_URL, SCOPUS_API_KEY, SCOPUS_ENABLED

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger("test_scopus_integration")


class TestScopusIntegration(unittest.TestCase):
    """Integration tests for Scopus citation source."""

    @classmethod
    def setUpClass(cls):
        """Set up test fixtures that are reused across test methods."""
        # Check if Scopus is enabled and API key is available
        cls.api_key = os.environ.get("SCOPUS_API_KEY", SCOPUS_API_KEY)
        cls.scopus_enabled = SCOPUS_ENABLED and cls.api_key
        
        if not cls.scopus_enabled:
            logger.warning("Scopus tests skipped: SCOPUS_ENABLED=false or no API key available")
            return
        
        # Initialize the Scopus client
        cls.client = ScopusClient(
            api_url=SCOPUS_API_URL,
            api_key=cls.api_key,
            institutional_token=os.environ.get("SCOPUS_INSTITUTIONAL_TOKEN", ""),
            rate_limit=5.0
        )
        
        # Test DOIs - these should be stable, well-known papers
        cls.test_dois = [
            "10.1038/s41586-020-2008-3",  # COVID-19 paper (high citations)
            "10.1093/bioinformatics/btab213"  # Bioinformatics paper
        ]
    
    def setUp(self):
        """Set up test fixtures."""
        if not self.scopus_enabled:
            self.skipTest("Scopus tests skipped: SCOPUS_ENABLED=false or no API key available")
    
    def test_scopus_client_direct(self):
        """Test direct use of Scopus client."""
        for doi in self.test_dois:
            success, data = self.client.get_citation_data(doi)
            
            self.assertTrue(success, f"Failed to get Scopus data for DOI {doi}")
            self.assertIn("total_citations", data, "Citation count missing from response")
            self.assertIn("metadata", data, "Metadata missing from response")
            self.assertIn("source", data, "Source identifier missing from response")
            self.assertEqual(data["source"], "scopus", "Incorrect source in response")
            
            # Check that citation count is reasonable (greater than zero for these known papers)
            self.assertGreater(data["total_citations"], 0, f"Expected citations > 0 for {doi}")
            
            logger.info(f"Scopus citation count for {doi}: {data['total_citations']}")
    
    def test_scopus_client_batch(self):
        """Test batch citation retrieval."""
        results = self.client.batch_get_citation_data(self.test_dois)
        
        self.assertEqual(len(results), len(self.test_dois), "Not all DOIs processed")
        
        for doi in self.test_dois:
            self.assertIn(doi, results, f"DOI {doi} missing from batch results")
            data = results[doi]
            self.assertNotIn("error", data, f"Error in batch result for {doi}: {data.get('error', '')}")
            self.assertIn("total_citations", data, "Citation count missing from response")
            
            logger.info(f"Batch Scopus citation count for {doi}: {data['total_citations']}")
    
    def test_citation_registry_integration(self):
        """Test that Scopus is properly registered in the citation registry."""
        # Check if Scopus is in available sources
        sources = get_available_sources()
        self.assertIn("scopus", sources, "Scopus not registered in citation sources")
        
        # Get Scopus client instance from registry
        scopus = get_citation_source("scopus")
        self.assertIsNotNone(scopus, "Failed to get Scopus client from registry")
        self.assertIsInstance(scopus, ScopusClient, "Registry returned wrong client type")
    
    def test_citation_aggregator_with_scopus(self):
        """Test citation aggregation with Scopus as a source."""
        aggregator = CitationAggregator()
        
        for doi in self.test_dois:
            aggregated_data = aggregator.collect_citations(doi)
            
            self.assertNotIn("error", aggregated_data, 
                            f"Error in aggregated data for {doi}: {aggregated_data.get('error', '')}")
            self.assertIn("total_citations", aggregated_data, "Citation count missing from aggregated data")
            self.assertIn("sources_used", aggregated_data, "Sources used missing from aggregated data")
            self.assertIn("scopus", aggregated_data["sources_used"], 
                         "Scopus not listed in sources used for aggregation")
            
            # Scopus should be primary source since it has highest priority
            self.assertEqual(aggregated_data.get("primary_source"), "scopus", 
                            "Scopus not selected as primary source despite having highest priority")
            
            logger.info(f"Aggregated citation count for {doi}: {aggregated_data['total_citations']}")


if __name__ == "__main__":
    unittest.main()