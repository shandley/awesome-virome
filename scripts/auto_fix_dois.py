#!/usr/bin/env python3
"""
Auto-fix DOI Formatting in Metadata Files

This script automatically fixes common DOI formatting issues in metadata files:
1. Removes trailing parentheses, periods, and other non-DOI characters
2. Removes markdown link formatting
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

# DOI regex pattern
DOI_PATTERN = r"^10\.\d{4,9}/[-._;()/:A-Z0-9]+$"

# Common DOI formatting issues
DOI_ISSUES = [
    (r'\)\.?$', ''),  # Remove trailing ).
    (r'\.+$', ''),    # Remove trailing periods
    (r'\)$', ''),     # Remove trailing )
    (r'\].*$', ''),   # Remove markdown links
]

def find_metadata_files() -> List[str]:
    """Find all JSON metadata files in the repository"""
    metadata_files = []
    
    # Search in metadata/bioinformatics
    if os.path.exists(BIOINFORMATICS_DIR):
        metadata_files.extend(glob.glob(os.path.join(BIOINFORMATICS_DIR, "*.json")))
    
    # Search in metadata/academic_impact
    if os.path.exists(ACADEMIC_IMPACT_DIR):
        metadata_files.extend(glob.glob(os.path.join(ACADEMIC_IMPACT_DIR, "*.json")))
    
    # Search in metadata root
    metadata_files.extend(glob.glob(os.path.join(METADATA_DIR, "*.json")))
    
    return metadata_files

def fix_doi(doi: str) -> Tuple[str, bool]:
    """
    Fix common DOI formatting issues
    
    Args:
        doi: The DOI string to fix
        
    Returns:
        Tuple of (fixed DOI, whether a change was made)
    """
    if not doi:
        return doi, False
    
    original_doi = doi
    
    # Fix markdown link format like 10.1093/gigascience/giae020](https://doi.org/10.1093/gigascience/giae020)
    if '](' in doi:
        doi = doi.split('](')[0]
    
    # Apply common fixes
    for pattern, replacement in DOI_ISSUES:
        doi = re.sub(pattern, replacement, doi)
    
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
        with open(file_path, 'r') as f:
            data = json.load(f)
    except json.JSONDecodeError:
        logger.error(f"Failed to parse JSON in {file_path}")
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