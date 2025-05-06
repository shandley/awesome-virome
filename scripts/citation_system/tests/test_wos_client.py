#!/usr/bin/env python3
"""
Tests for the Web of Science API client.
"""

import datetime
import unittest
from unittest.mock import MagicMock, patch

import requests

from scripts.citation_system.api.sources.wos_client import WebOfScienceClient


class TestWebOfScienceClient(unittest.TestCase):
    """Tests for the Web of Science API client."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.api_key = "test_api_key"
        self.api_secret = "test_api_secret"
        self.client = WebOfScienceClient(
            api_url="https://api.clarivate.com/apis/wos-starter/v1",
            api_key=self.api_key,
            api_secret=self.api_secret,
            rate_limit=10.0
        )
        
        # Mock authentication for all tests
        self.client.auth_token = "mock_token"
        self.client.auth_token_expiry = datetime.datetime.now() + datetime.timedelta(hours=1)
        
        # Mock common methods for all tests
        self.client._make_request = MagicMock()
    
    def test_ensure_auth_token_with_valid_token(self):
        """Test that _ensure_auth_token returns True when a valid token exists."""
        # Token is already set up in setUp
        self.assertTrue(self.client._ensure_auth_token())
    
    @patch('requests.post')
    def test_ensure_auth_token_with_expired_token(self, mock_post):
        """Test that _ensure_auth_token refreshes an expired token."""
        # Set up an expired token
        self.client.auth_token_expiry = datetime.datetime.now() - datetime.timedelta(minutes=5)
        
        # Mock the API response
        mock_response = MagicMock()
        mock_response.json.return_value = {
            "access_token": "new_mock_token",
            "expires_in": 3600
        }
        mock_response.raise_for_status = MagicMock()
        mock_post.return_value = mock_response
        
        # Call the method
        result = self.client._ensure_auth_token()
        
        # Verify the result
        self.assertTrue(result)
        self.assertEqual(self.client.auth_token, "new_mock_token")
        self.assertTrue(self.client.auth_token_expiry > datetime.datetime.now())
        
        # Verify the API call
        mock_post.assert_called_once()
        args, kwargs = mock_post.call_args
        self.assertEqual(args[0], "https://api.clarivate.com/apis/wos-starter/auth/v1")
        self.assertEqual(kwargs['json']['client_id'], self.api_key)
        self.assertEqual(kwargs['json']['client_secret'], self.api_secret)
    
    @patch('requests.post')
    def test_ensure_auth_token_with_error(self, mock_post):
        """Test that _ensure_auth_token handles errors correctly."""
        # Set up an expired token
        self.client.auth_token = None
        self.client.auth_token_expiry = None
        
        # Mock the API response
        mock_post.side_effect = requests.exceptions.RequestException("API error")
        
        # Call the method
        result = self.client._ensure_auth_token()
        
        # Verify the result
        self.assertFalse(result)
        self.assertIsNone(self.client.auth_token)
        self.assertIsNone(self.client.auth_token_expiry)
        
        # Verify the API call
        mock_post.assert_called_once()
    
    def test_get_citation_data(self):
        """Test the citation data retrieval workflow."""
        # Mock authentication
        self.client._ensure_auth_token = MagicMock(return_value=True)
        
        # Mock API response
        mock_response = {
            "data": [
                {
                    "uid": "WOS:000123456789",
                    "title": {"value": "Test Title"},
                    "authors": [
                        {"name": "Doe, John"},
                        {"name": "Smith, Jane"}
                    ],
                    "source": {"sourceTitle": "Test Journal"},
                    "publicationDate": {"year": "2020"},
                    "metrics": {"timesCited": 42}
                }
            ]
        }
        self.client._make_request.return_value = (True, mock_response)
        
        # Mock yearly citations
        self.client._get_yearly_citations = MagicMock(return_value={"2020": 10, "2021": 15, "2022": 17})
        
        # Test successful citation data retrieval
        success, data = self.client.get_citation_data("10.1234/test", use_cache=False)
        
        # Verify the API call
        self.client._make_request.assert_called_once()
        args, kwargs = self.client._make_request.call_args
        self.assertIn("doi/10.1234%2Ftest", kwargs.get("endpoint", ""))
        self.assertEqual(kwargs.get("headers", {}).get("Authorization"), "Bearer mock_token")
        
        # Verify the result
        self.assertTrue(success)
        self.assertEqual(data["total_citations"], 42)
        self.assertEqual(data["year"], 2020)
        self.assertEqual(data["metadata"]["title"], "Test Title")
        self.assertEqual(data["metadata"]["journal"], "Test Journal")
        self.assertEqual(len(data["metadata"]["authors"]), 2)
        self.assertEqual(data["metadata"]["authors"][0], "Doe, John")
        
        # Verify yearly citations were fetched
        self.client._get_yearly_citations.assert_called_once_with("10.1234/test", "WOS:000123456789", False)
    
    def test_get_citation_data_with_no_results(self):
        """Test citation data retrieval with no results."""
        # Mock authentication
        self.client._ensure_auth_token = MagicMock(return_value=True)
        
        # Mock API response with no data
        self.client._make_request.return_value = (True, {"data": []})
        
        # Test retrieval with no results
        success, data = self.client.get_citation_data("10.1234/nonexistent", use_cache=False)
        
        # Verify the result
        self.assertFalse(success)
        self.assertTrue("not found in Web of Science" in data)
    
    def test_get_citation_data_with_auth_failure(self):
        """Test citation data retrieval with authentication failure."""
        # Mock authentication failure
        self.client._ensure_auth_token = MagicMock(return_value=False)
        
        # Test retrieval with auth failure
        success, data = self.client.get_citation_data("10.1234/test", use_cache=False)
        
        # Verify the result
        self.assertFalse(success)
        self.assertTrue("Could not authenticate" in data)
        
        # Verify no API call was made
        self.client._make_request.assert_not_called()
    
    def test_get_citation_data_with_api_error(self):
        """Test citation data retrieval with API error."""
        # Mock authentication
        self.client._ensure_auth_token = MagicMock(return_value=True)
        
        # Mock API error
        self.client._make_request.return_value = (False, "API error")
        
        # Test retrieval with API error
        success, data = self.client.get_citation_data("10.1234/test", use_cache=False)
        
        # Verify the result
        self.assertFalse(success)
        self.assertTrue("Failed to get Web of Science data" in data)
    
    def test_get_yearly_citations(self):
        """Test yearly citation data retrieval."""
        # Mock authentication
        self.client._ensure_auth_token = MagicMock(return_value=True)
        
        # Mock API response
        mock_response = {
            "data": [
                {"year": "2020", "count": 10},
                {"year": "2021", "count": 15},
                {"year": "2022", "count": 17}
            ]
        }
        self.client._make_request.return_value = (True, mock_response)
        
        # Test successful yearly citation retrieval
        citations_by_year = self.client._get_yearly_citations("10.1234/test", "WOS:000123456789", use_cache=False)
        
        # Verify the API call
        self.client._make_request.assert_called_once()
        args, kwargs = self.client._make_request.call_args
        self.assertIn("document/WOS:000123456789/citations/year", kwargs.get("endpoint", ""))
        
        # Verify the result
        self.assertEqual(len(citations_by_year), 3)
        self.assertEqual(citations_by_year["2020"], 10)
        self.assertEqual(citations_by_year["2021"], 15)
        self.assertEqual(citations_by_year["2022"], 17)
    
    def test_get_yearly_citations_with_error(self):
        """Test yearly citation data retrieval with error."""
        # Mock authentication
        self.client._ensure_auth_token = MagicMock(return_value=True)
        
        # Mock API error
        self.client._make_request.return_value = (False, "API error")
        
        # Test retrieval with API error
        citations_by_year = self.client._get_yearly_citations("10.1234/test", "WOS:000123456789", use_cache=False)
        
        # Verify the result is an empty dictionary
        self.assertEqual(len(citations_by_year), 0)
    
    def test_batch_get_citation_data(self):
        """Test batch citation data retrieval."""
        # Mock get_citation_data to return different results
        def mock_get_citation_data(doi, use_cache):
            if doi == "10.1234/success":
                return True, {"doi": doi, "total_citations": 42}
            else:
                return False, f"Error for {doi}"
        
        self.client.get_citation_data = MagicMock(side_effect=mock_get_citation_data)
        self.client._ensure_auth_token = MagicMock(return_value=True)
        
        # Test batch retrieval
        dois = ["10.1234/success", "10.1234/error", "10.5281/zenodo.12345"]
        results = self.client.batch_get_citation_data(dois, use_cache=False)
        
        # Verify authentication was checked once
        self.client._ensure_auth_token.assert_called_once()
        
        # Verify get_citation_data was called for non-Zenodo DOIs
        self.assertEqual(self.client.get_citation_data.call_count, 2)
        
        # Verify results structure
        self.assertEqual(len(results), 3)
        self.assertIn("10.1234/success", results)
        self.assertIn("10.1234/error", results)
        self.assertIn("10.5281/zenodo.12345", results)
        
        # Verify successful result
        self.assertEqual(results["10.1234/success"]["total_citations"], 42)
        
        # Verify error result
        self.assertIn("error", results["10.1234/error"])
        
        # Verify Zenodo result
        self.assertIn("error", results["10.5281/zenodo.12345"])
        self.assertTrue("Zenodo DOI" in results["10.5281/zenodo.12345"]["error"])


if __name__ == "__main__":
    unittest.main()