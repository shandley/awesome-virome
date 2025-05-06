#!/usr/bin/env python3
"""
Tests for the citation source heatmap generator.
"""

import os
import sys
import json
import pytest
import tempfile
from pathlib import Path

# Add the scripts directory to the path
scripts_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(scripts_dir)

# Import the module
from citation_source_heatmap import (
    load_impact_data,
    get_source_data,
    generate_heatmap,
    generate_source_distribution_chart
)


class TestCitationSourceHeatmap:
    """Tests for the citation source heatmap generator."""
    
    def setup_method(self):
        """Set up test environment."""
        # Create a temporary file for the impact data
        self.temp_dir = tempfile.TemporaryDirectory()
        self.temp_path = Path(self.temp_dir.name)
        
        # Create a sample impact data file
        self.impact_data = {
            "tools": [
                {
                    "name": "Tool A",
                    "total_citations": 100,
                    "citation_source": "scopus",
                    "yearly_citation_source": "scopus",
                    "category": "Category 1"
                },
                {
                    "name": "Tool B",
                    "total_citations": 50,
                    "citation_source": "wos",
                    "yearly_citation_source": "wos",
                    "category": "Category 1"
                },
                {
                    "name": "Tool C",
                    "total_citations": 25,
                    "citation_source": "icite",
                    "yearly_citation_source": None,
                    "category": "Category 2"
                }
            ]
        }
        
        self.impact_data_path = self.temp_path / "impact_data.json"
        with open(self.impact_data_path, 'w') as f:
            json.dump(self.impact_data, f)
    
    def teardown_method(self):
        """Clean up after tests."""
        self.temp_dir.cleanup()
    
    def test_get_source_data(self):
        """Test extraction of source data from impact data."""
        df, categories = get_source_data(self.impact_data)
        
        # Check DataFrame structure
        assert not df.empty
        assert 'Tool' in df.columns
        assert 'Citation Source' in df.columns
        assert 'Yearly Citation Source' in df.columns
        assert 'Category' in df.columns
        assert 'Total Citations' in df.columns
        
        # Check data contents
        assert len(df) == 3
        assert set(df['Tool']) == {"Tool A", "Tool B", "Tool C"}
        assert set(df['Citation Source']) == {"scopus", "wos", "icite"}
        
        # Check categories
        assert "Category 1" in categories
        assert "Category 2" in categories
        assert set(categories["Category 1"]) == {"Tool A", "Tool B"}
        assert set(categories["Category 2"]) == {"Tool C"}
    
    def test_generate_heatmap(self):
        """Test heatmap generation."""
        df, _ = get_source_data(self.impact_data)
        output_path = self.temp_path / "test_heatmap.png"
        
        # Generate the heatmap
        generate_heatmap(df, str(output_path))
        
        # Check if file was created
        assert output_path.exists()
        assert output_path.stat().st_size > 0
    
    def test_generate_source_distribution_chart(self):
        """Test source distribution chart generation."""
        df, _ = get_source_data(self.impact_data)
        output_path = self.temp_path / "test_distribution.png"
        
        # Generate the chart
        generate_source_distribution_chart(df, str(output_path))
        
        # Check if file was created
        assert output_path.exists()
        assert output_path.stat().st_size > 0


if __name__ == "__main__":
    pytest.main(["-v", __file__])