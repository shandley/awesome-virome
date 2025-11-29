#!/usr/bin/env python3
"""
Tests for citation source attribution in the citation aggregator.
"""

import pytest
from unittest.mock import MagicMock, patch
import logging
import json
from datetime import datetime

from ..collectors.citation_aggregator import CitationAggregator

# Configure logging
logging.basicConfig(level=logging.DEBUG)


class TestCitationSourceAttribution:
    """Tests for citation source attribution in the citation aggregator."""
    
    def setup_method(self):
        """Set up test environment."""
        self.aggregator = CitationAggregator()
    
    def test_source_attribution_in_aggregation(self):
        """Test that citation sources are properly attributed in the aggregated data."""
        # Create sample source data
        source_data = {
            "scopus": {
                "total_citations": 100,
                "citations_by_year": {"2020": 30, "2021": 40, "2022": 30}
            },
            "wos": {
                "total_citations": 90,
                "citations_by_year": {"2020": 25, "2021": 35, "2022": 30}
            },
            "icite": {
                "total_citations": 95,
                "citations_by_year": {"2020": 28, "2021": 37, "2022": 30}
            }
        }
        
        # Aggregate the data
        aggregated = self.aggregator._aggregate_citation_data(source_data, "10.1000/test.123")
        
        # Check for source attribution fields
        assert "citation_source" in aggregated
        assert "yearly_citation_source" in aggregated
        
        # With the given data, scopus should be the citation source (highest priority)
        assert aggregated["citation_source"] == "scopus"
        assert aggregated["yearly_citation_source"] == "scopus"
        
        # Verify other metadata was preserved
        assert aggregated["total_citations"] == 100
        assert aggregated["citations_by_year"] == {"2020": 30, "2021": 40, "2022": 30}
    
    def test_different_sources_for_citations_and_yearly_data(self):
        """Test when different sources are used for citation count and yearly data."""
        # Create sample source data where one source has total citations 
        # but no yearly data, and another has yearly data
        source_data = {
            "scopus": {
                "total_citations": 100,
                "citations_by_year": {}  # No yearly data
            },
            "icite": {
                "total_citations": 90,
                "citations_by_year": {"2020": 25, "2021": 35, "2022": 30}
            }
        }
        
        # Aggregate the data
        aggregated = self.aggregator._aggregate_citation_data(source_data, "10.1000/test.123")
        
        # Check for source attribution fields
        assert "citation_source" in aggregated
        assert "yearly_citation_source" in aggregated
        
        # Scopus should be the citation source, but iCite should provide yearly data
        assert aggregated["citation_source"] == "scopus"
        assert aggregated["yearly_citation_source"] == "icite"
        
        # Verify data integrity
        assert aggregated["total_citations"] == 100  # From scopus
        assert aggregated["citations_by_year"] == {"2020": 25, "2021": 35, "2022": 30}  # From icite
    
    def test_backward_compatibility(self):
        """Test that backward compatibility with existing field names is maintained."""
        # Create sample source data
        source_data = {
            "scopus": {
                "total_citations": 100,
                "citations_by_year": {"2020": 30, "2021": 40, "2022": 30}
            }
        }
        
        # Aggregate the data
        aggregated = self.aggregator._aggregate_citation_data(source_data, "10.1000/test.123")
        
        # Check for both old and new field names
        assert "citation_source" in aggregated
        assert "primary_source" in aggregated
        assert "yearly_citation_source" in aggregated
        assert "yearly_data_source" in aggregated
        
        # Fields should have the same values
        assert aggregated["citation_source"] == aggregated["primary_source"]
        assert aggregated["yearly_citation_source"] == aggregated["yearly_data_source"]