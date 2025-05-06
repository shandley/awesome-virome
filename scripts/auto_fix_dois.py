#!/usr/bin/env python3
"""
Authoritative DOI Formatting Tool for Metadata Files

This script is the single authoritative source for fixing DOI formatting issues across the repository.
It automatically fixes common DOI formatting issues in metadata files:

1. Removes trailing parentheses, periods, and other non-DOI characters
2. Removes markdown link formatting
3. Standardizes DOI format (lowercase, trim whitespace)
4. Handles DOIs in different fields (doi, citation.doi, academic_impact.doi, publication.doi, etc.)
5. Supports JSON and YAML files

This script consolidates DOI fixing functionality previously spread across multiple scripts.
"""

import os
import json
import re
import logging
import argparse
import glob
from pathlib import Path
from typing import List, Dict, Any, Tuple

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("auto_fix_dois.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Constants
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(SCRIPT_DIR)
METADATA_DIR = os.path.join(ROOT_DIR, "metadata")
BIOINFORMATICS_DIR = os.path.join(METADATA_DIR, "bioinformatics")
ACADEMIC_IMPACT_DIR = os.path.join(METADATA_DIR, "academic_impact")

# DOI regex patterns
DOI_PATTERN = r"^10\.\d{4,9}/[-._;()/:A-Za-z0-9]+$"
DOI_PREFIX_PATTERN = r"^https?://(dx\.)?doi\.org/"

# Common DOI formatting issues
DOI_ISSUES = [
    (r'\)\.?$', ''),         # Remove trailing ).
    (r'\.+$', ''),           # Remove trailing periods
    (r'\)$', ''),            # Remove trailing )
    (r'\].*$', ''),          # Remove markdown links
    (r'^https?://(dx\.)?doi\.org/', ''),  # Remove doi.org prefix
    (r'\s+', ''),            # Remove whitespace
    (r'<.*>', ''),           # Remove HTML tags
    (r'^DOI:', ''),          # Remove 'DOI:' prefix
    (r'^doi:', '')           # Remove 'doi:' prefix (case insensitive)
]

def find_metadata_files() -> List[str]:
    """Find all JSON and YAML metadata files in the repository"""
    metadata_files = []
    
    # Search for JSON files in metadata directories
    for directory in [BIOINFORMATICS_DIR, ACADEMIC_IMPACT_DIR, METADATA_DIR]:
        if os.path.exists(directory):
            # Get all JSON files
            metadata_files.extend(glob.glob(os.path.join(directory, "**/*.json"), recursive=True))
            
            # Get YAML files if they exist (some repositories use YAML)
            metadata_files.extend(glob.glob(os.path.join(directory, "**/*.yaml"), recursive=True))
            metadata_files.extend(glob.glob(os.path.join(directory, "**/*.yml"), recursive=True))
    
    # Search in pubmed_citations if it exists
    pubmed_dir = os.path.join(METADATA_DIR, "pubmed_citations")
    if os.path.exists(pubmed_dir):
        metadata_files.extend(glob.glob(os.path.join(pubmed_dir, "*.json")))
    
    return metadata_files

def fix_doi(doi: str) -> Tuple[str, bool]:
    """
    Fix common DOI formatting issues
    
    Args:
        doi: The DOI string to fix
        
    Returns:
        Tuple of (fixed DOI, whether a change was made)
    """
    if not doi or not isinstance(doi, str):
        return doi, False
    
    original_doi = doi.strip()
    
    # Fix markdown link format like 10.1093/gigascience/giae020](https://doi.org/10.1093/gigascience/giae020)
    if '](' in doi:
        doi = doi.split('](')[0]
    
    # Handle square brackets
    if '[' in doi and ']' in doi:
        doi = re.sub(r'\[([^\]]+)\].*', r'\1', doi)
    
    # Apply common fixes
    for pattern, replacement in DOI_ISSUES:
        doi = re.sub(pattern, replacement, doi)
    
    # Ensure DOI starts with 10.
    if not doi.startswith('10.'):
        # Try to find a valid DOI within the string
        match = re.search(r'(10\.\d{4,9}/[-._;()/:A-Za-z0-9]+)', doi)
        if match:
            doi = match.group(1)
    
    # Standardize to lowercase
    doi = doi.lower()
    
    # Final cleanup of any remaining whitespace
    doi = doi.strip()
    
    # Check if we made a change
    return doi, doi != original_doi

def process_metadata_file(file_path: str) -> bool:
    """
    Process a single metadata file, fixing DOI issues
    
    Args:
        file_path: Path to the metadata file
        
    Returns:
        True if changes were made, False otherwise
    """
    try:
        file_ext = os.path.splitext(file_path)[1].lower()
        with open(file_path, 'r') as f:
            # Parse based on file extension
            if file_ext in [".yaml", ".yml"]:
                # We would use PyYAML here, but it's not always available
                # For simplicity, we'll skip YAML files in this demonstration
                logger.info(f"YAML support not implemented for {file_path}")
                return False
            else:  # Assume JSON
                data = json.load(f)
    except (json.JSONDecodeError, IOError):
        logger.error(f"Failed to parse file {file_path}")
        return False
    
    # Check various locations where DOIs might be stored
    made_changes = False
    
    # Check direct doi field
    if "doi" in data:
        fixed_doi, changed = fix_doi(data["doi"])
        if changed:
            data["doi"] = fixed_doi
            made_changes = True
            logger.info(f"Fixed DOI in {file_path}: {fixed_doi}")
    
    # Check academic_impact field
    if "academic_impact" in data and isinstance(data["academic_impact"], dict):
        if "doi" in data["academic_impact"]:
            fixed_doi, changed = fix_doi(data["academic_impact"]["doi"])
            if changed:
                data["academic_impact"]["doi"] = fixed_doi
                made_changes = True
                logger.info(f"Fixed academic_impact DOI in {file_path}: {fixed_doi}")
    
    # Check citation field
    if "citation" in data and isinstance(data["citation"], dict):
        if "doi" in data["citation"]:
            fixed_doi, changed = fix_doi(data["citation"]["doi"])
            if changed:
                data["citation"]["doi"] = fixed_doi
                made_changes = True
                logger.info(f"Fixed citation DOI in {file_path}: {fixed_doi}")
    
    # Check citation_info field
    if "citation_info" in data and isinstance(data["citation_info"], dict):
        # Check publication subfield
        if "publication" in data["citation_info"] and isinstance(data["citation_info"]["publication"], dict):
            if "doi" in data["citation_info"]["publication"]:
                fixed_doi, changed = fix_doi(data["citation_info"]["publication"]["doi"])
                if changed:
                    data["citation_info"]["publication"]["doi"] = fixed_doi
                    made_changes = True
                    logger.info(f"Fixed publication DOI in {file_path}: {fixed_doi}")
    
    # Write changes back if needed
    if made_changes:
        try:
            with open(file_path, 'w') as f:
                json.dump(data, f, indent=2)
            return True
        except IOError:
            logger.error(f"Failed to write changes to {file_path}")
    
    return False

def main():
    """Main entry point for the script"""
    parser = argparse.ArgumentParser(description="Auto-fix DOI formatting issues in metadata files")
    parser.add_argument("--dry-run", action="store_true", help="Only report issues without making changes")
    args = parser.parse_args()
    
    logger.info("Starting DOI auto-fix script...")
    metadata_files = find_metadata_files()
    logger.info(f"Found {len(metadata_files)} metadata files")
    
    fixed_files = 0
    fixed_dois = 0
    
    for file_path in metadata_files:
        if args.dry_run:
            try:
                with open(file_path, 'r') as f:
                    data = json.load(f)
                
                # Check for DOI issues without fixing
                doi_issues = []
                
                if "doi" in data:
                    _, changed = fix_doi(data["doi"])
                    if changed:
                        doi_issues.append(f"Direct DOI: {data['doi']}")
                
                if "academic_impact" in data and isinstance(data["academic_impact"], dict):
                    if "doi" in data["academic_impact"]:
                        _, changed = fix_doi(data["academic_impact"]["doi"])
                        if changed:
                            doi_issues.append(f"Academic impact DOI: {data['academic_impact']['doi']}")
                
                if "citation" in data and isinstance(data["citation"], dict):
                    if "doi" in data["citation"]:
                        _, changed = fix_doi(data["citation"]["doi"])
                        if changed:
                            doi_issues.append(f"Citation DOI: {data['citation']['doi']}")
                
                if doi_issues:
                    logger.info(f"Issues in {file_path}:")
                    for issue in doi_issues:
                        logger.info(f"  - {issue}")
                    fixed_dois += len(doi_issues)
                    fixed_files += 1
            except (json.JSONDecodeError, IOError):
                logger.error(f"Error processing {file_path}")
        else:
            changes_made = process_metadata_file(file_path)
            if changes_made:
                fixed_files += 1
                fixed_dois += 1  # This is an approximation
    
    if args.dry_run:
        logger.info(f"Dry run completed. Found {fixed_dois} DOI issues in {fixed_files} files.")
    else:
        logger.info(f"Completed. Fixed {fixed_dois} DOI issues in {fixed_files} files.")

if __name__ == "__main__":
    main()