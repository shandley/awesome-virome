#!/usr/bin/env python3
"""
Script to check that version information in README.md is properly formatted.
"""

import re
import sys
import logging
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

def load_readme(readme_path):
    """Load the README.md file."""
    try:
        with open(readme_path, 'r') as f:
            readme_content = f.read()
        return readme_content
    except FileNotFoundError as e:
        logger.error(f"Error loading README.md: {e}")
        sys.exit(1)

def check_version_format(readme_content):
    """Check for properly formatted version information in README.md."""
    # Define regular expressions for version patterns
    version_patterns = [
        (r'\[v\d+\.\d+\.\d+, \d{4}\]', 'Semantic version with year [v1.2.3, 2023]'),
        (r'\[v\d+\.\d+, \d{4}\]', 'Major.minor version with year [v1.2, 2023]'),
        (r'\[v\d+, \d{4}\]', 'Major version with year [v1, 2023]')
    ]
    
    # Find all version strings in the README
    version_matches = []
    for line in readme_content.splitlines():
        for pattern, desc in version_patterns:
            matches = re.findall(pattern, line)
            for match in matches:
                version_matches.append((match, desc, line.strip()))
    
    # Find potential version strings that don't match the expected format
    incorrect_versions = []
    potential_version_pattern = r'(\[v?\d+\.\d+(\.\d+)?[^,\]]*\])'
    potential_matches = re.findall(potential_version_pattern, readme_content)
    
    for match, _ in potential_matches:
        # Check if this is a valid version format
        is_valid = False
        for pattern, _ in version_patterns:
            if re.fullmatch(pattern, match):
                is_valid = True
                break
        
        if not is_valid:
            # Find the line containing this match
            for line in readme_content.splitlines():
                if match in line:
                    incorrect_versions.append((match, line.strip()))
                    break
    
    return version_matches, incorrect_versions

def check_update_dates(readme_content):
    """Check for properly formatted update dates in README.md."""
    # Define the pattern for update dates: [Updated: MM/YYYY]
    update_date_pattern = r'\[Updated: (0[1-9]|1[0-2])/20\d{2}\]'
    
    # Find all update dates in the README
    update_dates = []
    incorrect_dates = []
    
    for line in readme_content.splitlines():
        # Check for correct format
        correct_matches = re.findall(update_date_pattern, line)
        if correct_matches:
            for _ in correct_matches:
                update_dates.append(line.strip())
        
        # Check for potentially incorrect formats
        potential_pattern = r'\[Updated:([^\]]+)\]'
        potential_matches = re.findall(potential_pattern, line)
        
        for match in potential_matches:
            if not re.match(r' (0[1-9]|1[0-2])/20\d{2}$', match):
                incorrect_dates.append((f"[Updated:{match}]", line.strip()))
    
    return update_dates, incorrect_dates

def main():
    """Main function to check version information in README.md."""
    repo_root = Path(__file__).parent.parent
    readme_path = repo_root / 'README.md'
    
    # Load README content
    readme_content = load_readme(readme_path)
    
    # Check version formats
    version_matches, incorrect_versions = check_version_format(readme_content)
    logger.info(f"Found {len(version_matches)} correctly formatted version strings")
    
    # Check update dates
    update_dates, incorrect_dates = check_update_dates(readme_content)
    logger.info(f"Found {len(update_dates)} correctly formatted update dates")
    
    # Report issues
    has_issues = False
    
    if incorrect_versions:
        has_issues = True
        logger.error(f"Found {len(incorrect_versions)} incorrectly formatted version strings:")
        for version, line in incorrect_versions:
            logger.error(f"  - {version} in line: {line}")
    
    if incorrect_dates:
        has_issues = True
        logger.error(f"Found {len(incorrect_dates)} incorrectly formatted update dates:")
        for date, line in incorrect_dates:
            logger.error(f"  - {date} in line: {line}")
    
    if has_issues:
        logger.error("Version information check failed. Please fix the issues above.")
        sys.exit(1)
    else:
        logger.info("All version information is properly formatted.")
        sys.exit(0)

if __name__ == "__main__":
    main()