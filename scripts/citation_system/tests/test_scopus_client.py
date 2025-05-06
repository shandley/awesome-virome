#!/usr/bin/env python3
"""
Tests for the Scopus API client.
"""

import os
import unittest
from unittest.mock import MagicMock, patch

import requests

from scripts.citation_system.api.sources.scopus_client import ScopusClient


class TestScopusClient(unittest.TestCase):
    """Tests for the Scopus API client."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.api_key = "test_api_key"
        self.client = ScopusClient(
            api_url="https://api.elsevier.com",
            api_key=self.api_key,
            institutional_token="test_institutional_token",
            rate_limit=10.0
        )
        
        # Mock common methods for all tests
        self.client._make_request = MagicMock()
    
    def test_get_scopus_id(self):
        """Test getting Scopus ID from DOI."""
        # Mock Scopus API response
        mock_response = {
            "search-results": {
                "entry": [
                    {
                        "dc:identifier": "SCOPUS_ID:85123456789"
                    }
                ]
            }
        }
        self.client._make_request.return_value = (True, mock_response)
        
        # Test successful ID retrieval
        scopus_id = self.client._get_scopus_id("10.1234/test", use_cache=False)
        
        # Verify the request was made correctly
        self.client._make_request.assert_called_once()
        args, kwargs = self.client._make_request.call_args
        self.assertTrue("DOI(10.1234/test)" in kwargs.get("endpoint", ""))
        
        # Verify the result
        self.assertEqual(scopus_id, "85123456789")
        
        # Test with empty response
        self.client._make_request.reset_mock()
        self.client._make_request.return_value = (True, {"search-results": {"entry": []}})
        
        scopus_id = self.client._get_scopus_id("10.1234/nonexistent", use_cache=False)
        self.assertIsNone(scopus_id)
        
        # Test with failed request
        self.client._make_request.reset_mock()
        self.client._make_request.return_value = (False, "API error")
        
        scopus_id = self.client._get_scopus_id("10.1234/error", use_cache=False)
        self.assertIsNone(scopus_id)
    
    def test_get_citations_by_scopus_id(self):
        """Test getting citation data using Scopus ID."""
        # Mock Scopus API response
        mock_response = {
            "abstracts-retrieval-response": {
                "coredata": {
                    "citedby-count": "42"
                },
                "item": {
                    "bibrecord": {
                        "head": {
                            "citation": {
                                "source": {
                                    "sourcetitle": "Test Journal"
                                },
                                "publicationyear": {"@first": "2020"},
                                "titletext": {"titletext": "Test Title"}
                            },
                            "author-group": {
                                "author": [
                                    {
                                        "ce:given-name": "John",
                                        "ce:surname": "Doe"
                                    },
                                    {
                                        "ce:given-name": "Jane",
                                        "ce:surname": "Smith"
                                    }
                                ]
                            }
                        }
                    }
                }
            }
        }
        self.client._make_request.return_value = (True, mock_response)
        
        # Mock yearly citations
        self.client._get_yearly_citations = MagicMock(return_value={"2020": 10, "2021": 15, "2022": 17})
        
        # Test successful citation data retrieval
        success, data = self.client._get_citations_by_scopus_id("85123456789", "10.1234/test", use_cache=False)
        
        # Verify the request was made correctly
        self.client._make_request.assert_called_once()
        args, kwargs = self.client._make_request.call_args
        self.assertTrue("85123456789" in kwargs.get("endpoint", ""))
        
        # Verify the result
        self.assertTrue(success)
        self.assertEqual(data["total_citations"], 42)
        self.assertEqual(data["year"], 2020)
        self.assertEqual(data["metadata"]["title"], "Test Title")
        self.assertEqual(data["metadata"]["journal"], "Test Journal")
        self.assertEqual(len(data["metadata"]["authors"]), 2)
        self.assertEqual(data["metadata"]["authors"][0], "Doe, John")
        
        # Test with empty response
        self.client._make_request.reset_mock()
        self.client._make_request.return_value = (True, {})
        
        success, data = self.client._get_citations_by_scopus_id("85123456789", "10.1234/nonexistent", use_cache=False)
        self.assertFalse(success)
        
        # Test with failed request
        self.client._make_request.reset_mock()
        self.client._make_request.return_value = (False, "API error")
        
        success, data = self.client._get_citations_by_scopus_id("85123456789", "10.1234/error", use_cache=False)
        self.assertFalse(success)
    
    def test_get_yearly_citations(self):
        """Test getting yearly citation breakdown."""
        # Mock Scopus API response
        mock_response = {
            "abstract-citations-response": {
                "citeInfoMatrix": {
                    "citeInfoMatrixXML": {
                        "citationMatrix": {
                            "citeInfo": [
                                {"year": "2020", "valueByYear": "10"},
                                {"year": "2021", "valueByYear": "15"},
                                {"year": "2022", "valueByYear": "17"}
                            ]
                        }
                    }
                }
            }
        }
        self.client._make_request.return_value = (True, mock_response)
        
        # Test successful yearly citation retrieval
        citations_by_year = self.client._get_yearly_citations("85123456789", use_cache=False)
        
        # Verify the request was made correctly
        self.client._make_request.assert_called_once()
        args, kwargs = self.client._make_request.call_args
        self.assertTrue("85123456789" in kwargs.get("endpoint", ""))
        
        # Verify the result
        self.assertEqual(len(citations_by_year), 3)
        self.assertEqual(citations_by_year["2020"], 10)
        self.assertEqual(citations_by_year["2021"], 15)
        self.assertEqual(citations_by_year["2022"], 17)
        
        # Test with empty response
        self.client._make_request.reset_mock()
        self.client._make_request.return_value = (True, {})
        
        citations_by_year = self.client._get_yearly_citations("85123456789", use_cache=False)
        self.assertEqual(len(citations_by_year), 0)
        
        # Test with failed request
        self.client._make_request.reset_mock()
        self.client._make_request.return_value = (False, "API error")
        
        citations_by_year = self.client._get_yearly_citations("85123456789", use_cache=False)
        self.assertEqual(len(citations_by_year), 0)
    
    def test_get_citation_data(self):
        """Test the full citation data retrieval workflow."""
        # Mock Scopus ID lookup
        self.client._get_scopus_id = MagicMock(return_value="85123456789")
        
        # Mock citation data retrieval
        mock_citation_data = {
            "doi": "10.1234/test",
            "source": "scopus",
            "total_citations": 42,
            "year": 2020,
            "metadata": {
                "title": "Test Title",
                "authors": ["Doe, John", "Smith, Jane"],
                "journal": "Test Journal"
            }
        }
        self.client._get_citations_by_scopus_id = MagicMock(return_value=(True, mock_citation_data))
        
        # Test successful full workflow
        success, data = self.client.get_citation_data("10.1234/test", use_cache=False)
        
        # Verify both methods were called
        self.client._get_scopus_id.assert_called_once_with("10.1234/test", False)
        self.client._get_citations_by_scopus_id.assert_called_once_with("85123456789", "10.1234/test", False)
        
        # Verify the result
        self.assertTrue(success)
        self.assertEqual(data, mock_citation_data)
        
        # Test with no Scopus ID found
        self.client._get_scopus_id.reset_mock()
        self.client._get_citations_by_scopus_id.reset_mock()
        self.client._get_scopus_id.return_value = None
        
        success, data = self.client.get_citation_data("10.1234/nonexistent", use_cache=False)
        
        # Verify only ID lookup was called
        self.client._get_scopus_id.assert_called_once_with("10.1234/nonexistent", False)
        self.client._get_citations_by_scopus_id.assert_not_called()
        
        # Verify the result
        self.assertFalse(success)
        self.assertTrue("not found in Scopus" in data)
    
    def test_batch_get_citation_data(self):
        """Test batch citation data retrieval."""
        # Mock get_citation_data to return different results
        def mock_get_citation_data(doi, use_cache):
            if doi == "10.1234/success":
                return True, {"doi": doi, "total_citations": 42}
            else:
                return False, f"Error for {doi}"
        
        self.client.get_citation_data = MagicMock(side_effect=mock_get_citation_data)
        
        # Test batch retrieval
        dois = ["10.1234/success", "10.1234/error"]
        results = self.client.batch_get_citation_data(dois, use_cache=False)
        
        # Verify calls were made for both DOIs
        self.assertEqual(self.client.get_citation_data.call_count, 2)
        
        # Verify results structure
        self.assertEqual(len(results), 2)
        self.assertIn("10.1234/success", results)
        self.assertIn("10.1234/error", results)
        
        # Verify successful result
        self.assertEqual(results["10.1234/success"]["total_citations"], 42)
        
        # Verify error result
        self.assertIn("error", results["10.1234/error"])


if __name__ == "__main__":
    unittest.main()