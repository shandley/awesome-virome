#!/usr/bin/env python3
"""
Unit tests for the iCite client.
"""

import json
import unittest
from unittest.mock import patch

import requests

from ..api.sources.icite_client import ICiteClient


class TestICiteClient(unittest.TestCase):
    """Tests for the iCite API client."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.client = ICiteClient(
            api_url="https://icite.od.nih.gov/api",
            rate_limit=10.0
        )
        
        # Example DOI for testing
        self.test_doi = "10.1093/bioinformatics/btab213"  # This is a real DOI for a bioinformatics paper
        
        # Example PMID for testing
        self.test_pmid = "33760067"  # PMID corresponding to the DOI above
    
    @patch("..api.sources.icite_client.BaseAPIClient._make_request")
    def test_get_citation_data_success(self, mock_make_request):
        """Test successful citation data retrieval."""
        # Mock response data
        mock_response = {
            "data": [
                {
                    "doi": "10.1093/bioinformatics/btab213",
                    "title": "CheckV: assessing the quality of metagenome-assembled viral genomes",
                    "authors": "Nayfach S|Camargo AP|Schulz F|Eloe-Fadrosh E|Roux S|Kyrpides NC",
                    "year": 2021,
                    "journal": "Bioinformatics",
                    "citation_count": 123,
                    "relative_citation_ratio": 2.5,
                    "expected_citations_per_year": 10.0,
                    "field_citation_rate": 4.0,
                    "is_clinical": False,
                    "pmid": "33760067"
                }
            ]
        }
        
        # Configure the mock to return the test data
        mock_make_request.return_value = (True, mock_response)
        
        # Call the method being tested
        success, data = self.client.get_citation_data(self.test_doi)
        
        # Assert that the request was made correctly
        mock_make_request.assert_called_once()
        
        # Check the results
        self.assertTrue(success)
        self.assertEqual(data["doi"], self.test_doi)
        self.assertEqual(data["total_citations"], 123)
        self.assertEqual(data["rcr"], 2.5)
        self.assertEqual(data["metadata"]["title"], "CheckV: assessing the quality of metagenome-assembled viral genomes")
        self.assertEqual(len(data["metadata"]["authors"]), 6)
    
    @patch("..api.sources.icite_client.BaseAPIClient._make_request")
    def test_get_citation_data_not_found(self, mock_make_request):
        """Test when a DOI is not found in iCite."""
        # Mock response with empty data
        mock_response = {"data": []}
        
        # Configure the mock
        mock_make_request.return_value = (True, mock_response)
        
        # Call the method being tested
        success, data = self.client.get_citation_data("10.1234/nonexistent")
        
        # Assert that the request was made
        mock_make_request.assert_called_once()
        
        # Check the results
        self.assertFalse(success)
        self.assertTrue("not found" in data)
    
    @patch("..api.sources.icite_client.BaseAPIClient._make_request")
    def test_get_citation_data_api_error(self, mock_make_request):
        """Test handling of API errors."""
        # Configure mock to simulate an API error
        mock_make_request.return_value = (False, "API error message")
        
        # Call the method being tested
        success, data = self.client.get_citation_data(self.test_doi)
        
        # Assert that the request was made
        mock_make_request.assert_called_once()
        
        # Check the results
        self.assertFalse(success)
        self.assertTrue("Failed to get iCite data" in data)
    
    @patch("..api.sources.icite_client.BaseAPIClient._make_request")
    def test_search_by_pmid(self, mock_make_request):
        """Test searching by PMID."""
        # Mock response data
        mock_response = {
            "data": [
                {
                    "doi": self.test_doi,
                    "pmid": self.test_pmid,
                    "title": "CheckV: assessing the quality of metagenome-assembled viral genomes",
                    "citation_count": 123
                }
            ]
        }
        
        # Configure mock for first call (PMID search)
        mock_make_request.return_value = (True, mock_response)
        
        # For the actual test, we'll patch the get_citation_data method
        with patch.object(self.client, "get_citation_data") as mock_get_citation_data:
            # Set up the return value for get_citation_data
            mock_get_citation_data.return_value = (True, {"doi": self.test_doi, "total_citations": 123})
            
            # Call the method being tested
            success, data = self.client.search_by_pmid(self.test_pmid)
            
            # Check that the PMID request was made
            mock_make_request.assert_called_once()
            
            # Check that get_citation_data was called with the correct DOI
            mock_get_citation_data.assert_called_once_with(self.test_doi, True)
            
            # Check the results
            self.assertTrue(success)
            self.assertEqual(data["doi"], self.test_doi)
    
    @patch("..api.sources.icite_client.BaseAPIClient._make_request")
    def test_batch_get_citation_data(self, mock_make_request):
        """Test batch processing of multiple DOIs."""
        # Create a list of test DOIs
        test_dois = [
            "10.1093/bioinformatics/btab213",  # First test DOI
            "10.1038/s41587-021-00980-x"       # Another real DOI
        ]
        
        # For simplicity, we'll patch the get_citation_data method directly
        with patch.object(self.client, "get_citation_data") as mock_get_citation_data:
            # Set up the return values
            mock_get_citation_data.side_effect = [
                (True, {"doi": test_dois[0], "total_citations": 123}),
                (True, {"doi": test_dois[1], "total_citations": 45})
            ]
            
            # Call the method being tested
            results = self.client.batch_get_citation_data(test_dois)
            
            # Check that get_citation_data was called twice
            self.assertEqual(mock_get_citation_data.call_count, 2)
            
            # Check the results
            self.assertEqual(len(results), 2)
            self.assertEqual(results[test_dois[0]]["total_citations"], 123)
            self.assertEqual(results[test_dois[1]]["total_citations"], 45)


if __name__ == "__main__":
    unittest.main()