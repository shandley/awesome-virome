#!/usr/bin/env python3
"""
Tests for the DOI validator.
"""

import unittest
from pathlib import Path
import sys

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent))

from scripts.citation_system.validators.doi_validator import DOIValidator


class DOIValidatorTest(unittest.TestCase):
    """Test cases for the DOI Validator."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.validator = DOIValidator()
    
    def test_valid_doi_formats(self):
        """Test validation of DOIs with valid formats."""
        valid_dois = [
            "10.1093/nar/gkac273",
            "10.1038/s41586-021-03819-2",
            "10.1093/bioinformatics/btac776",
            "10.1093/bioinformatics/btz265",
            "10.1073/pnas.2122636119"
        ]
        
        for doi in valid_dois:
            self.assertTrue(
                self.validator.is_valid_format(doi), 
                f"DOI should be valid: {doi}"
            )
    
    def test_invalid_doi_formats(self):
        """Test validation of DOIs with invalid formats."""
        invalid_dois = [
            "",
            "not a doi",
            "11.1093/nar/gkac273",  # Invalid prefix
            "10.1093/nar/gkac273)",  # Trailing parenthesis
            "10.1093/nar/gkac273.",  # Trailing period
            "10.1038/s41586-021-03819-2]"  # Trailing bracket
        ]
        
        for doi in invalid_dois:
            self.assertFalse(
                self.validator.is_valid_format(doi), 
                f"DOI should be invalid: {doi}"
            )
    
    def test_doi_cleaning(self):
        """Test cleaning of DOIs with common issues."""
        test_cases = [
            # (input, expected output)
            ("https://doi.org/10.1093/nar/gkac273", "10.1093/nar/gkac273"),
            ("doi.org/10.1093/nar/gkac273", "10.1093/nar/gkac273"),
            ("DOI: 10.1093/nar/gkac273", "10.1093/nar/gkac273"),
            ("10.1093/nar/gkac273)", "10.1093/nar/gkac273"),
            ("10.1093/nar/gkac273.", "10.1093/nar/gkac273"),
            ("10.1093/gigascience/giae020](https://doi.org/10.1093/gigascience/giae020)",
             "10.1093/gigascience/giae020")
        ]
        
        for input_doi, expected_output in test_cases:
            cleaned_doi = self.validator._clean_doi(input_doi)
            self.assertEqual(
                cleaned_doi, 
                expected_output,
                f"DOI cleaning failed: {input_doi} -> {cleaned_doi} (expected: {expected_output})"
            )
    
    def test_doi_normalization(self):
        """Test normalization of DOIs."""
        test_cases = [
            # (input, expected output)
            ("https://doi.org/10.1093/nar/gkac273", "10.1093/nar/gkac273"),
            ("10.1093/nar/gkac273)", "10.1093/nar/gkac273"),
            ("invalid doi", "")  # Invalid DOI returns empty string
        ]
        
        for input_doi, expected_output in test_cases:
            normalized_doi = self.validator.normalize_doi(input_doi)
            self.assertEqual(
                normalized_doi, 
                expected_output,
                f"DOI normalization failed: {input_doi} -> {normalized_doi}"
            )


if __name__ == "__main__":
    unittest.main()