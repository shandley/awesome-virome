#!/usr/bin/env python3
"""
Script to fix corrupted JSON files in the cache directory
"""

import os
import json
import glob
import logging
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)

def fix_json_files(directory):
    """Find and fix corrupted JSON files in the directory"""
    fixed_count = 0
    examined_count = 0
    
    cache_dir = Path(directory)
    if not cache_dir.exists() or not cache_dir.is_dir():
        logger.error(f"Directory {directory} does not exist or is not a directory")
        return
    
    for filepath in glob.glob(os.path.join(directory, "*.json")):
        examined_count += 1
        try:
            # Try to load the JSON file
            with open(filepath, 'r') as f:
                data = json.load(f)
            # If it loads successfully, it's not corrupted
        except json.JSONDecodeError as e:
            logger.info(f"Found corrupted file: {filepath} - {str(e)}")
            # File is corrupted, let's fix it
            with open(filepath, 'r') as f:
                content = f.read()
            
            # Extract cache_date if it exists
            import re
            date_match = re.search(r'"cache_date":\s*"([^"]+)"', content)
            cache_date = None
            if date_match:
                cache_date = date_match.group(1)
            else:
                # Use current time if no cache date found
                from datetime import datetime
                cache_date = datetime.now().isoformat()
            
            # Create a proper empty response
            fixed_content = f'{{\n  "cache_date": "{cache_date}",\n  "data": []\n}}'
            
            with open(filepath, 'w') as f:
                f.write(fixed_content)
            fixed_count += 1
            logger.info(f"Fixed corrupted file: {filepath}")
                
    logger.info(f"Examined {examined_count} files and fixed {fixed_count} corrupted JSON files")
    return fixed_count

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Fix corrupted JSON files in the cache directory")
    parser.add_argument("--directory", default="metadata/cache", help="Directory containing JSON cache files")
    args = parser.parse_args()
    
    # Get absolute path if relative was provided
    directory = os.path.abspath(args.directory)
    logger.info(f"Scanning directory: {directory}")
    
    fixed = fix_json_files(directory)
    if fixed > 0:
        logger.info(f"Successfully fixed {fixed} corrupted JSON files")
    else:
        logger.info("No corrupted JSON files found or fixed")