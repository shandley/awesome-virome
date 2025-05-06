#!/usr/bin/env python3
"""
Test script for comprehensive_citation_data.py

This script tests that comprehensive_citation_data.py does not generate synthetic data
and only uses real citation data.
"""

import os
import sys
import json
import unittest
from unittest.mock import patch, MagicMock
import tempfile
from pathlib import Path

# Add parent directory to path so we can import the script
parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_dir)

# Import the module to test
import comprehensive_citation_data

class TestComprehensiveCitationData(unittest.TestCase):
    """Test class for comprehensive_citation_data.py"""
    
    def setUp(self):
        """Set up test environment"""
        # Create temporary directory for test files
        self.temp_dir = tempfile.TemporaryDirectory()
        self.test_dir = Path(self.temp_dir.name)
        
        # Create test directories
        self.reports_dir = self.test_dir / "reports" / "citations"
        self.academic_dir = self.test_dir / "metadata" / "academic_impact"
        self.bioinformatics_dir = self.test_dir / "metadata" / "bioinformatics"
        
        self.reports_dir.mkdir(parents=True, exist_ok=True)
        self.academic_dir.mkdir(parents=True, exist_ok=True)
        self.bioinformatics_dir.mkdir(parents=True, exist_ok=True)
        
        # Create sample data files
        self.create_sample_data()
        
        # Mock the paths in the module
        self.original_base_dir = comprehensive_citation_data.BASE_DIR
        self.original_reports_dir = comprehensive_citation_data.REPORTS_DIR
        self.original_academic_dir = comprehensive_citation_data.ACADEMIC_IMPACT_DIR
        self.original_bioinformatics_dir = comprehensive_citation_data.BIOINFORMATICS_DIR
        self.original_output_file = comprehensive_citation_data.OUTPUT_FILE
        
        comprehensive_citation_data.BASE_DIR = self.test_dir
        comprehensive_citation_data.REPORTS_DIR = self.reports_dir
        comprehensive_citation_data.ACADEMIC_IMPACT_DIR = self.academic_dir
        comprehensive_citation_data.BIOINFORMATICS_DIR = self.bioinformatics_dir
        comprehensive_citation_data.OUTPUT_FILE = self.test_dir / "impact_data.json"
    
    def tearDown(self):
        """Tear down test environment"""
        # Restore original paths
        comprehensive_citation_data.BASE_DIR = self.original_base_dir
        comprehensive_citation_data.REPORTS_DIR = self.original_reports_dir
        comprehensive_citation_data.ACADEMIC_IMPACT_DIR = self.original_academic_dir
        comprehensive_citation_data.BIOINFORMATICS_DIR = self.original_bioinformatics_dir
        comprehensive_citation_data.OUTPUT_FILE = self.original_output_file
        
        # Clean up temporary directory
        self.temp_dir.cleanup()
    
    def create_sample_data(self):
        """Create sample data files for testing"""
        # Create a sample most_cited_report.json
        most_cited_data = {
            "most_cited_tools": [
                {
                    "name": "Tool1",
                    "url": "https://example.com/tool1",
                    "doi": "10.1234/tool1",
                    "total_citations": 100,
                    "influential_citations": 25
                },
                {
                    "name": "Tool2",
                    "url": "https://example.com/tool2",
                    "doi": "10.1234/tool2",
                    "total_citations": 75,
                    "influential_citations": 15
                }
            ],
            "average_citations_per_tool": 87.5
        }
        with open(self.reports_dir / "most_cited_report.json", "w") as f:
            json.dump(most_cited_data, f)
        
        # Create a sample citation_trends_report.json
        trends_data = {
            "tool_yearly_data": {
                "Tool1": {
                    "2020": 20,
                    "2021": 30,
                    "2022": 50
                },
                "Tool2": {
                    "2020": 15,
                    "2021": 25,
                    "2022": 35
                }
            }
        }
        with open(self.reports_dir / "citation_trends_report.json", "w") as f:
            json.dump(trends_data, f)
        
        # Create a sample citation_reports_summary.json
        summary_data = {
            "total_tools_analyzed": 2,
            "total_tools_with_citations": 2,
            "total_citations": 175
        }
        with open(self.reports_dir / "citation_reports_summary.json", "w") as f:
            json.dump(summary_data, f)
        
        # Create sample academic impact data
        academic_tool1 = {
            "name": "Tool1",
            "url": "https://example.com/tool1",
            "doi": "10.1234/tool1",
            "citation_metrics": {
                "total_citations": 100,
                "influential_citations": 25,
                "citations_by_year": {
                    "2020": 20,
                    "2021": 30,
                    "2022": 50
                }
            }
        }
        with open(self.academic_dir / "Tool1.json", "w") as f:
            json.dump(academic_tool1, f)
        
        # Create sample bioinformatics data
        bioinformatics_tool2 = {
            "name": "Tool2",
            "url": "https://example.com/tool2",
            "academic_impact": {
                "doi": "10.1234/tool2",
                "total_citations": 75,
                "influential_citations": 15,
                "citations_by_year": {
                    "2020": 15,
                    "2021": 25,
                    "2022": 35
                }
            }
        }
        with open(self.bioinformatics_dir / "Tool2.json", "w") as f:
            json.dump(bioinformatics_tool2, f)
    
    def test_get_tool_category_does_not_use_heuristics(self):
        """Test that get_tool_category doesn't use heuristics"""
        # Test with various tool names that could trigger heuristic categorization
        tool_names = [
            "PhageFinder",  # Contains "phage"
            "HostPredictor",  # Contains "host"
            "VirusTaxonomy",  # Contains "tax"
            "VirusIdentifier",  # Contains "ident"
            "ViralAssembler",  # Contains "assembl"
            "GeneAnnotator",  # Contains "annot"
            "MetaVirome",  # Contains "meta"
            "StructureAnalyzer"  # Contains "struct"
        ]
        
        # All should return "Other Tools" if no tool_data provided
        for name in tool_names:
            category = comprehensive_citation_data.get_tool_category(name)
            self.assertEqual(category, "Other Tools", 
                            f"get_tool_category is using heuristics for '{name}'")
        
        # Should use category from tool_data if provided
        tool_data = {"category": "Test Category"}
        category = comprehensive_citation_data.get_tool_category("TestTool", tool_data)
        self.assertEqual(category, "Test Category", 
                        "get_tool_category is not using provided category")
    
    def test_collect_all_citation_data(self):
        """Test that collect_all_citation_data collects real data only"""
        # Run the function to collect data
        tools_data = comprehensive_citation_data.collect_all_citation_data()
        
        # Verify that the function collected the expected data
        self.assertIn("Tool1", tools_data, "Tool1 not found in collected data")
        self.assertIn("Tool2", tools_data, "Tool2 not found in collected data")
        
        # Verify citation counts match the input data (no modification)
        self.assertEqual(tools_data["Tool1"].get("total_citations"), 100, 
                        "Tool1 citation count was modified")
        self.assertEqual(tools_data["Tool2"].get("total_citations"), 75, 
                        "Tool2 citation count was modified")
        
        # Verify yearly citations match the input data (no synthetic distribution)
        self.assertEqual(tools_data["Tool1"].get("citations_by_year", {}).get("2020"), 20, 
                        "Tool1 yearly citations were modified")
        self.assertEqual(tools_data["Tool2"].get("citations_by_year", {}).get("2020"), 15, 
                        "Tool2 yearly citations were modified")
    
    def test_main_function(self):
        """Test the main function that generates impact_data.json"""
        # Run the main function
        comprehensive_citation_data.main()
        
        # Check if impact_data.json was created
        impact_data_path = self.test_dir / "impact_data.json"
        self.assertTrue(impact_data_path.exists(), "impact_data.json was not created")
        
        # Load and verify the content
        with open(impact_data_path, "r") as f:
            impact_data = json.load(f)
        
        # Verify that tools are included
        self.assertTrue(isinstance(impact_data.get("tools"), list), 
                       "impact_data.json does not contain a tools list")
        
        # Verify tool count
        self.assertEqual(len(impact_data.get("tools", [])), 2, 
                        "impact_data.json does not have the expected number of tools")
        
        # Verify citation totals
        self.assertEqual(impact_data.get("total_citations"), 175, 
                        "impact_data.json does not have the expected total_citations")
        
        # Verify no synthetic values were added
        for tool in impact_data.get("tools", []):
            # Check that citation counts weren't altered
            if tool["name"] == "Tool1":
                self.assertEqual(tool.get("total_citations"), 100, 
                                "Tool1 citation count was modified")
            elif tool["name"] == "Tool2":
                self.assertEqual(tool.get("total_citations"), 75, 
                                "Tool2 citation count was modified")

if __name__ == "__main__":
    unittest.main()