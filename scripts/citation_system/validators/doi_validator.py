#!/usr/bin/env python3
"""
DOI validation and normalization utilities.
"""

import logging
import re
import time
from typing import Dict, List, Tuple

import requests

from ..config import DOI_PATTERN, DOI_RESOLVER_URL, MAX_RETRIES, RETRY_DELAY
from ..utils.logging_utils import setup_logging

logger = logging.getLogger(__name__)


class DOIValidator:
    """Validates and normalizes DOIs."""
    
    def __init__(self):
        """Initialize the DOI validator."""
        self.pattern = re.compile(DOI_PATTERN, re.IGNORECASE)
        self.session = requests.Session()
        self.headers = {
            "Accept": "application/json",
            "User-Agent": "AwesomeVirome/1.0 (https://github.com/shandley/awesome-virome; mailto:your.email@example.com)"
        }
    
    def is_valid_format(self, doi: str) -> bool:
        """
        Check if a DOI has a valid format.
        
        Args:
            doi: DOI to validate
        
        Returns:
            True if the DOI format is valid, False otherwise
        """
        if not doi:
            return False
        
        # Clean up DOI - remove any URL prefix, trailing characters, etc.
        cleaned_doi = self._clean_doi(doi)
        
        # Check against regex pattern - only care about the cleaned version
        return bool(self.pattern.match(cleaned_doi))
    
    def _clean_doi(self, doi: str) -> str:
        """
        Clean a DOI by removing URL prefixes and other common issues.
        
        Args:
            doi: DOI to clean
        
        Returns:
            Cleaned DOI
        """
        if not doi:
            return ""
        
        # Remove leading/trailing whitespace
        doi = doi.strip()
        
        # Remove DOI URL prefix
        if doi.startswith('https://doi.org/'):
            doi = doi[len('https://doi.org/'):]
        elif doi.startswith('http://doi.org/'):
            doi = doi[len('http://doi.org/'):]
        elif doi.startswith('doi.org/'):
            doi = doi[len('doi.org/'):]
        elif doi.startswith('DOI:'):
            doi = doi[len('DOI:'):].strip()
        
        # Remove trailing punctuation that's not part of the DOI
        # Handle multiple trailing punctuation (e.g., ")).")
        while doi and doi[-1] in '.),;':
            doi = doi[:-1]
        
        # Handle special cases like multiple closing parentheses
        if doi and doi.endswith(')))'):
            doi = doi[:-3]
        elif doi and doi.endswith('))'):
            doi = doi[:-2]
        elif doi and doi.endswith(')'):
            doi = doi[:-1]
        
        # Remove markdown link formatting
        if '](' in doi:
            doi = doi.split('](')[0]
        
        return doi
    
    def normalize_doi(self, doi: str) -> str:
        """
        Normalize a DOI to a standard format.
        
        Args:
            doi: DOI to normalize
        
        Returns:
            Normalized DOI or empty string if invalid
        """
        cleaned_doi = self._clean_doi(doi)
        
        if not self.is_valid_format(cleaned_doi):
            logger.warning(f"Invalid DOI format: {doi} (cleaned: {cleaned_doi})")
            return ""
        
        return cleaned_doi
    
    def is_resolvable(self, doi: str, max_retries: int = MAX_RETRIES) -> bool:
        """
        Check if a DOI is resolvable (exists in the DOI system).
        
        Args:
            doi: DOI to check
            max_retries: Maximum number of retries for failed requests
        
        Returns:
            True if the DOI is resolvable, False otherwise
        """
        if not doi:
            return False
        
        # Clean and validate the DOI format first
        cleaned_doi = self._clean_doi(doi)
        if not self.is_valid_format(cleaned_doi):
            return False
        
        # Build the URL to check
        url = f"{DOI_RESOLVER_URL}api/handles/{cleaned_doi}"
        
        # Try to resolve the DOI
        retries = 0
        while retries <= max_retries:
            try:
                response = self.session.get(
                    url, 
                    headers=self.headers,
                    timeout=10
                )
                
                if response.status_code == 200:
                    data = response.json()
                    return data.get('responseCode') == 1
                
                # Handle specific error codes
                if response.status_code in [404, 400]:
                    # DOI not found or bad request
                    return False
                
                # Server errors or rate limiting, retry after delay
                if response.status_code in [429, 500, 502, 503, 504]:
                    retries += 1
                    if retries <= max_retries:
                        time.sleep(RETRY_DELAY * retries)
                        continue
                    return False
                
                # Other errors
                return False
                
            except requests.RequestException:
                retries += 1
                if retries <= max_retries:
                    time.sleep(RETRY_DELAY * retries)
                    continue
                return False
        
        return False
    
    def full_validation(self, doi: str) -> Dict[str, bool]:
        """
        Perform full DOI validation (format and resolvability).
        
        Args:
            doi: DOI to validate
        
        Returns:
            Dictionary with validation results
        """
        cleaned_doi = self._clean_doi(doi)
        
        format_valid = self.is_valid_format(cleaned_doi)
        
        # Only check resolvability if format is valid
        resolvable = self.is_resolvable(cleaned_doi) if format_valid else False
        
        return {
            "original": doi,
            "cleaned": cleaned_doi,
            "format_valid": format_valid,
            "resolvable": resolvable,
            "valid": format_valid and resolvable
        }
    
    def validate_multiple(
        self, 
        dois: List[str], 
        check_resolvable: bool = True
    ) -> List[Dict[str, bool]]:
        """
        Validate multiple DOIs.
        
        Args:
            dois: List of DOIs to validate
            check_resolvable: Whether to check if DOIs are resolvable
        
        Returns:
            List of validation results for each DOI
        """
        results = []
        
        for doi in dois:
            cleaned_doi = self._clean_doi(doi)
            format_valid = self.is_valid_format(cleaned_doi)
            
            result = {
                "original": doi,
                "cleaned": cleaned_doi,
                "format_valid": format_valid
            }
            
            if check_resolvable and format_valid:
                result["resolvable"] = self.is_resolvable(cleaned_doi)
                result["valid"] = format_valid and result["resolvable"]
            else:
                result["resolvable"] = None
                result["valid"] = format_valid
            
            results.append(result)
        
        return results
    
    def find_doi_issues(self, dois: List[str]) -> Tuple[int, List[Dict[str, str]]]:
        """
        Find issues with a list of DOIs.
        
        Args:
            dois: List of DOIs to check
        
        Returns:
            Tuple containing count of invalid DOIs and list of issues
        """
        issues = []
        invalid_count = 0
        
        for doi in dois:
            if not doi:
                continue
                
            cleaned_doi = self._clean_doi(doi)
            
            if doi != cleaned_doi:
                issues.append({
                    "original": doi,
                    "cleaned": cleaned_doi,
                    "issue": "formatting",
                    "message": f"DOI has formatting issues: {doi} -> {cleaned_doi}"
                })
            
            if not self.is_valid_format(cleaned_doi):
                invalid_count += 1
                issues.append({
                    "original": doi,
                    "cleaned": cleaned_doi,
                    "issue": "invalid_format",
                    "message": f"Invalid DOI format: {doi}"
                })
        
        return invalid_count, issues


def main():
    """Command-line interface for DOI validation."""
    import argparse
    import json
    import sys
    
    # Set up logging
    logger = setup_logging("doi_validator")
    
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Validate DOIs")
    parser.add_argument(
        "dois", 
        nargs="*", 
        help="DOIs to validate"
    )
    parser.add_argument(
        "-f", "--file", 
        help="File containing DOIs (one per line)"
    )
    parser.add_argument(
        "-r", "--resolve", 
        action="store_true",
        help="Check if DOIs are resolvable"
    )
    parser.add_argument(
        "-o", "--output", 
        help="Output file for validation results (JSON format)"
    )
    
    args = parser.parse_args()
    
    # Collect DOIs from arguments and/or file
    dois = args.dois
    
    if args.file:
        try:
            with open(args.file, 'r') as f:
                file_dois = [line.strip() for line in f if line.strip()]
                dois.extend(file_dois)
        except Exception as e:
            logger.error(f"Error reading DOI file: {e}")
            return 1
    
    if not dois:
        logger.error("No DOIs provided")
        return 1
    
    # Validate DOIs
    validator = DOIValidator()
    results = validator.validate_multiple(dois, args.resolve)
    
    # Summarize results
    valid_count = sum(1 for r in results if r["valid"])
    logger.info(f"DOIs validated: {len(results)}")
    logger.info(f"Valid DOIs: {valid_count}")
    logger.info(f"Invalid DOIs: {len(results) - valid_count}")
    
    # Output results
    if args.output:
        try:
            with open(args.output, 'w') as f:
                json.dump(results, f, indent=2)
            logger.info(f"Results saved to {args.output}")
        except Exception as e:
            logger.error(f"Error writing output file: {e}")
            return 1
    else:
        # Print results to console
        for result in results:
            status = "Valid" if result["valid"] else "Invalid"
            print(f"{status}: {result['original']} -> {result['cleaned']}")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())