#!/usr/bin/env python3
"""
Test script for auto_fix_dois.py

This script tests that auto_fix_dois.py correctly fixes DOI formatting issues.
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
import auto_fix_dois

class TestAutoFixDOIs(unittest.TestCase):
    """Test class for auto_fix_dois.py"""
    
    def setUp(self):
        """Set up test environment"""
        # Create temporary directory for test files
        self.temp_dir = tempfile.TemporaryDirectory()
        self.test_dir = Path(self.temp_dir.name)
        
        # Create test directories
        self.metadata_dir = self.test_dir / "metadata"
        self.academic_dir = self.metadata_dir / "academic_impact"
        self.bioinformatics_dir = self.metadata_dir / "bioinformatics"
        self.pubmed_dir = self.metadata_dir / "pubmed_citations"
        
        self.metadata_dir.mkdir(parents=True, exist_ok=True)
        self.academic_dir.mkdir(parents=True, exist_ok=True)
        self.bioinformatics_dir.mkdir(parents=True, exist_ok=True)
        self.pubmed_dir.mkdir(parents=True, exist_ok=True)
        
        # Create sample data files with DOI issues
        self.create_sample_data()
        
        # Mock the paths in the module
        self.original_script_dir = auto_fix_dois.SCRIPT_DIR
        self.original_root_dir = auto_fix_dois.ROOT_DIR
        self.original_metadata_dir = auto_fix_dois.METADATA_DIR
        self.original_bioinformatics_dir = auto_fix_dois.BIOINFORMATICS_DIR
        self.original_academic_impact_dir = auto_fix_dois.ACADEMIC_IMPACT_DIR
        
        auto_fix_dois.SCRIPT_DIR = str(self.test_dir)
        auto_fix_dois.ROOT_DIR = str(self.test_dir)
        auto_fix_dois.METADATA_DIR = str(self.metadata_dir)
        auto_fix_dois.BIOINFORMATICS_DIR = str(self.bioinformatics_dir)
        auto_fix_dois.ACADEMIC_IMPACT_DIR = str(self.academic_dir)
    
    def tearDown(self):
        """Tear down test environment"""
        # Restore original paths
        auto_fix_dois.SCRIPT_DIR = self.original_script_dir
        auto_fix_dois.ROOT_DIR = self.original_root_dir
        auto_fix_dois.METADATA_DIR = self.original_metadata_dir
        auto_fix_dois.BIOINFORMATICS_DIR = self.original_bioinformatics_dir
        auto_fix_dois.ACADEMIC_IMPACT_DIR = self.original_academic_impact_dir
        
        # Clean up temporary directory
        self.temp_dir.cleanup()
    
    def create_sample_data(self):
        """Create sample data files with DOI issues for testing"""
        # File with trailing period in DOI
        tool1 = {
            "name": "Tool1",
            "doi": "10.1234/tool1."
        }
        with open(self.metadata_dir / "Tool1.json", "w") as f:
            json.dump(tool1, f)
        
        # File with DOI in markdown link format
        tool2 = {
            "name": "Tool2",
            "doi": "10.1234/tool2](https://doi.org/10.1234/tool2)"
        }
        with open(self.metadata_dir / "Tool2.json", "w") as f:
            json.dump(tool2, f)
        
        # File with doi.org prefix
        tool3 = {
            "name": "Tool3",
            "doi": "https://doi.org/10.1234/tool3"
        }
        with open(self.academic_dir / "Tool3.json", "w") as f:
            json.dump(tool3, f)
        
        # File with DOI in academic_impact field
        tool4 = {
            "name": "Tool4",
            "academic_impact": {
                "doi": "DOI: 10.1234/tool4"
            }
        }
        with open(self.bioinformatics_dir / "Tool4.json", "w") as f:
            json.dump(tool4, f)
        
        # File with DOI in citation field
        tool5 = {
            "name": "Tool5",
            "citation": {
                "doi": "10.1234/tool5)"
            }
        }
        with open(self.bioinformatics_dir / "Tool5.json", "w") as f:
            json.dump(tool5, f)
        
        # File with DOI in citation_info.publication field
        tool6 = {
            "name": "Tool6",
            "citation_info": {
                "publication": {
                    "doi": "doi:10.1234/tool6"
                }
            }
        }
        with open(self.pubmed_dir / "Tool6.json", "w") as f:
            json.dump(tool6, f)
    
    def test_fix_doi(self):
        """Test the fix_doi function with various DOI formats"""
        test_cases = [
            # (input_doi, expected_doi, should_change)
            ("10.1234/abc.", "10.1234/abc", True),  # Trailing period
            ("10.1234/abc](https://doi.org/10.1234/abc)", "10.1234/abc", True),  # Markdown link
            ("https://doi.org/10.1234/abc", "10.1234/abc", True),  # doi.org prefix
            ("DOI: 10.1234/abc", "10.1234/abc", True),  # DOI prefix
            ("doi:10.1234/abc", "10.1234/abc", True),  # doi prefix lowercase
            ("10.1234/abc)", "10.1234/abc", True),  # Trailing parenthesis
            ("10.1234/abc", "10.1234/abc", False),  # Already correctly formatted
            ("10.1234/ABC", "10.1234/abc", True),  # Uppercase to lowercase
            (" 10.1234/abc ", "10.1234/abc", True),  # Whitespace
            (None, None, False),  # None input
            ("", "", False),  # Empty string
            ("not a doi", "not a doi", False),  # Invalid DOI
            ("[10.1234/abc]", "10.1234/abc", True),  # Square brackets
        ]
        
        for input_doi, expected_doi, should_change in test_cases:
            fixed_doi, changed = auto_fix_dois.fix_doi(input_doi)
            self.assertEqual(fixed_doi, expected_doi, 
                           f"fix_doi failed for '{input_doi}', got '{fixed_doi}'")
            self.assertEqual(changed, should_change, 
                           f"fix_doi change detection failed for '{input_doi}'")
    
    def test_find_metadata_files(self):
        """Test the find_metadata_files function"""
        # Run the function
        metadata_files = auto_fix_dois.find_metadata_files()
        
        # Convert paths to strings for easier comparison
        metadata_files = [str(p) for p in metadata_files]
        
        # Verify all expected files are found
        expected_files = [
            str(self.metadata_dir / "Tool1.json"),
            str(self.metadata_dir / "Tool2.json"),
            str(self.academic_dir / "Tool3.json"),
            str(self.bioinformatics_dir / "Tool4.json"),
            str(self.bioinformatics_dir / "Tool5.json"),
            str(self.pubmed_dir / "Tool6.json")
        ]
        
        for expected_file in expected_files:
            self.assertIn(expected_file, metadata_files, 
                         f"Expected file {expected_file} not found")
    
    def test_process_metadata_file(self):
        """Test the process_metadata_file function"""
        # Process each test file and verify the results
        
        # Tool1: Direct DOI with trailing period
        result = auto_fix_dois.process_metadata_file(str(self.metadata_dir / "Tool1.json"))
        self.assertTrue(result, "process_metadata_file failed to fix Tool1.json")
        with open(self.metadata_dir / "Tool1.json", "r") as f:
            tool1 = json.load(f)
        self.assertEqual(tool1["doi"], "10.1234/tool1", "DOI not fixed in Tool1.json")
        
        # Tool2: DOI in markdown link format
        result = auto_fix_dois.process_metadata_file(str(self.metadata_dir / "Tool2.json"))
        self.assertTrue(result, "process_metadata_file failed to fix Tool2.json")
        with open(self.metadata_dir / "Tool2.json", "r") as f:
            tool2 = json.load(f)
        self.assertEqual(tool2["doi"], "10.1234/tool2", "DOI not fixed in Tool2.json")
        
        # Tool3: DOI with doi.org prefix
        result = auto_fix_dois.process_metadata_file(str(self.academic_dir / "Tool3.json"))
        self.assertTrue(result, "process_metadata_file failed to fix Tool3.json")
        with open(self.academic_dir / "Tool3.json", "r") as f:
            tool3 = json.load(f)
        self.assertEqual(tool3["doi"], "10.1234/tool3", "DOI not fixed in Tool3.json")
        
        # Tool4: DOI in academic_impact field
        result = auto_fix_dois.process_metadata_file(str(self.bioinformatics_dir / "Tool4.json"))
        self.assertTrue(result, "process_metadata_file failed to fix Tool4.json")
        with open(self.bioinformatics_dir / "Tool4.json", "r") as f:
            tool4 = json.load(f)
        self.assertEqual(tool4["academic_impact"]["doi"], "10.1234/tool4", "DOI not fixed in Tool4.json")
        
        # Tool5: DOI in citation field
        result = auto_fix_dois.process_metadata_file(str(self.bioinformatics_dir / "Tool5.json"))
        self.assertTrue(result, "process_metadata_file failed to fix Tool5.json")
        with open(self.bioinformatics_dir / "Tool5.json", "r") as f:
            tool5 = json.load(f)
        self.assertEqual(tool5["citation"]["doi"], "10.1234/tool5", "DOI not fixed in Tool5.json")
        
        # Tool6: DOI in citation_info.publication field
        result = auto_fix_dois.process_metadata_file(str(self.pubmed_dir / "Tool6.json"))
        self.assertTrue(result, "process_metadata_file failed to fix Tool6.json")
        with open(self.pubmed_dir / "Tool6.json", "r") as f:
            tool6 = json.load(f)
        self.assertEqual(tool6["citation_info"]["publication"]["doi"], "10.1234/tool6", 
                        "DOI not fixed in Tool6.json")
    
    @patch('argparse.ArgumentParser.parse_args')
    def test_main_function(self, mock_parse_args):
        """Test the main function"""
        # Mock the command line arguments
        mock_args = MagicMock()
        mock_args.dry_run = False
        mock_parse_args.return_value = mock_args
        
        # Run the main function
        with patch('sys.stdout'): # suppress print output
            auto_fix_dois.main()
        
        # Verify all DOIs were fixed
        with open(self.metadata_dir / "Tool1.json", "r") as f:
            tool1 = json.load(f)
        self.assertEqual(tool1["doi"], "10.1234/tool1", "DOI not fixed in Tool1.json")
        
        with open(self.metadata_dir / "Tool2.json", "r") as f:
            tool2 = json.load(f)
        self.assertEqual(tool2["doi"], "10.1234/tool2", "DOI not fixed in Tool2.json")
        
        with open(self.academic_dir / "Tool3.json", "r") as f:
            tool3 = json.load(f)
        self.assertEqual(tool3["doi"], "10.1234/tool3", "DOI not fixed in Tool3.json")
        
        with open(self.bioinformatics_dir / "Tool4.json", "r") as f:
            tool4 = json.load(f)
        self.assertEqual(tool4["academic_impact"]["doi"], "10.1234/tool4", "DOI not fixed in Tool4.json")
        
        with open(self.bioinformatics_dir / "Tool5.json", "r") as f:
            tool5 = json.load(f)
        self.assertEqual(tool5["citation"]["doi"], "10.1234/tool5", "DOI not fixed in Tool5.json")
        
        with open(self.pubmed_dir / "Tool6.json", "r") as f:
            tool6 = json.load(f)
        self.assertEqual(tool6["citation_info"]["publication"]["doi"], "10.1234/tool6", 
                        "DOI not fixed in Tool6.json")

if __name__ == "__main__":
    unittest.main()