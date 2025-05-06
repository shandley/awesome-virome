#!/usr/bin/env python3
"""
Fix Citation Counts Script

This script fixes the citation counts for tools with DOIs by directly querying
the CrossRef API and updating the json files in the pubmed_citations directory.
This is a one-time fix to ensure all citation counts are correctly updated.

Usage:
    python fix_citation_counts.py
"""

import os
import sys
import json
import time
import logging
import glob
import urllib.request
import urllib.error
from pathlib import Path
from datetime import datetime
from typing import Dict, Any, List
from concurrent.futures import ThreadPoolExecutor, as_completed

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("fix_citation_counts.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Directory constants
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
METADATA_DIR = os.path.join(BASE_DIR, "metadata")
PUBMED_CITATIONS_DIR = os.path.join(METADATA_DIR, "pubmed_citations")
CROSSREF_API_URL = "https://api.crossref.org/works/"

def get_citation_count_from_crossref(doi: str) -> int:
    """
    Get citation count directly from CrossRef API.
    
    Args:
        doi: The DOI to look up
        
    Returns:
        int: The citation count or 0 if not found
    """
    if not doi:
        return 0
    
    url = f"{CROSSREF_API_URL}{doi}"
    headers = {
        "User-Agent": "AwesomeVirome/1.0 (https://github.com/scottx611x/awesome-virome; mailto:scott@airy.ai)"
    }
    
    # Add retries to handle transient failures
    max_retries = 3
    retry_delay = 5  # seconds
    
    for attempt in range(max_retries):
        try:
            logger.info(f"Querying CrossRef API for DOI: {doi} (attempt {attempt+1}/{max_retries})")
            
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=30) as response:
                data = json.loads(response.read().decode('utf-8'))
                message = data.get('message', {})
                citation_count = message.get('is-referenced-by-count', 0)
                logger.info(f"Found {citation_count} citations for DOI: {doi}")
                return citation_count
                
        except urllib.error.HTTPError as e:
            if e.code == 429:  # Too Many Requests
                retry_after = e.headers.get('Retry-After', retry_delay)
                try:
                    retry_after = int(retry_after)
                except (ValueError, TypeError):
                    retry_after = retry_delay
                
                logger.warning(f"Rate limited by CrossRef API. Waiting {retry_after}s before retry.")
                time.sleep(retry_after)
            else:
                logger.warning(f"Failed to get citation count for DOI {doi}: HTTP error {e.code}")
                if attempt < max_retries - 1:
                    logger.info(f"Retrying in {retry_delay}s...")
                    time.sleep(retry_delay)
                else:
                    return 0
                    
        except Exception as e:
            logger.error(f"Error getting citation count for DOI {doi}: {e}")
            if attempt < max_retries - 1:
                logger.info(f"Retrying in {retry_delay}s...")
                time.sleep(retry_delay)
            else:
                return 0
    
    return 0

def update_tool_citation(file_path: str) -> Dict[str, Any]:
    """
    Update citation count for a single tool.
    
    Args:
        file_path: Path to the tool's JSON file
        
    Returns:
        dict: The updated tool data
    """
    try:
        with open(file_path, 'r') as f:
            tool_data = json.load(f)
        
        name = tool_data.get('name', os.path.basename(file_path).replace('.json', ''))
        logger.info(f"Processing {name}")
        
        # Get citation info
        citation_info = tool_data.get('citation_info', {})
        publication = citation_info.get('publication', {})
        
        if not publication:
            logger.info(f"No publication info for {name}")
            return tool_data
        
        # Check if we have a DOI to work with
        doi = publication.get('doi')
        if not doi:
            logger.info(f"No DOI for {name}")
            return tool_data
        
        # Get current citation count
        current_count = publication.get('citation_count', 0)
        logger.info(f"{name} current citation count: {current_count}")
        
        # Get updated citation count from CrossRef
        updated_count = get_citation_count_from_crossref(doi)
        
        if updated_count > 0:
            # Update the citation count
            publication['citation_count'] = updated_count
            logger.info(f"Updated {name} citation count: {updated_count}")
        else:
            logger.info(f"No citation count found for {name} ({doi})")
        
        # Update and save the tool data
        tool_data['citation_info']['publication'] = publication
        tool_data['last_updated'] = datetime.now().isoformat()
        
        with open(file_path, 'w') as f:
            json.dump(tool_data, f, indent=2)
        
        # Add a small delay to avoid rate limiting
        time.sleep(2)
        
        return tool_data
        
    except Exception as e:
        logger.error(f"Error updating citation count for {file_path}: {e}")
        return {}

def fix_all_citation_counts():
    """
    Fix citation counts for all tools with DOIs.
    """
    # Get all tool files
    tool_files = glob.glob(os.path.join(PUBMED_CITATIONS_DIR, "*.json"))
    tool_files = [f for f in tool_files if "pubmed_citations.json" not in f and "summary.json" not in f and "citation_counts_summary.json" not in f]
    
    logger.info(f"Found {len(tool_files)} tool files to process")
    
    # Statistics
    stats = {
        "tools_processed": 0,
        "tools_with_doi": 0,
        "tools_with_citation_counts": 0,
        "tools_updated": 0,
        "total_citations": 0,
        "start_time": datetime.now()
    }
    
    # Create a list of high-priority tools to process first (these are known to have citation issues)
    priority_tools = ["VirSorter.json", "VirFinder.json", "Virosaurus.json"]
    priority_files = [os.path.join(PUBMED_CITATIONS_DIR, tool) for tool in priority_tools]
    
    # Filter out files that don't exist
    priority_files = [f for f in priority_files if os.path.exists(f)]
    
    # Remove priority files from the main list to avoid duplication
    remaining_files = [f for f in tool_files if f not in priority_files]
    
    # Process priority tools first
    logger.info(f"Processing {len(priority_files)} priority tools first")
    for file_path in priority_files:
        try:
            tool_name = os.path.basename(file_path).replace('.json', '')
            logger.info(f"Processing priority tool: {tool_name}")
            
            tool_data = update_tool_citation(file_path)
            
            # Update statistics
            if tool_data:
                stats["tools_processed"] += 1
                
                citation_info = tool_data.get('citation_info', {})
                publication = citation_info.get('publication', {})
                
                if publication and publication.get('doi'):
                    stats["tools_with_doi"] += 1
                    
                    citation_count = publication.get('citation_count', 0)
                    if citation_count > 0:
                        stats["tools_with_citation_counts"] += 1
                        stats["total_citations"] += citation_count
                        stats["tools_updated"] += 1
                        
        except Exception as e:
            logger.error(f"Error processing priority tool {file_path}: {e}")
    
    # Process remaining tools
    logger.info(f"Processing {len(remaining_files)} remaining tools")
    for file_path in remaining_files:
        try:
            tool_data = update_tool_citation(file_path)
            
            # Update statistics
            if tool_data:
                stats["tools_processed"] += 1
                
                citation_info = tool_data.get('citation_info', {})
                publication = citation_info.get('publication', {})
                
                if publication and publication.get('doi'):
                    stats["tools_with_doi"] += 1
                    
                    citation_count = publication.get('citation_count', 0)
                    if citation_count > 0:
                        stats["tools_with_citation_counts"] += 1
                        stats["total_citations"] += citation_count
                        stats["tools_updated"] += 1
                        
        except Exception as e:
            logger.error(f"Error processing {file_path}: {e}")
    
    # Calculate elapsed time
    elapsed_time = (datetime.now() - stats["start_time"]).total_seconds()
    
    # Report statistics
    logger.info(f"Citation count fix completed in {elapsed_time:.2f}s")
    logger.info(f"Processed {stats['tools_processed']} tools")
    logger.info(f"Found {stats['tools_with_doi']} tools with DOIs")
    logger.info(f"Updated {stats['tools_with_citation_counts']} tools with citation counts")
    logger.info(f"Total citations found: {stats['total_citations']}")
    
    # Update summary file
    summary_data = {
        "last_updated": datetime.now().isoformat(),
        "tools_processed": stats["tools_processed"],
        "tools_with_doi": stats["tools_with_doi"],
        "tools_with_citation_counts": stats["tools_with_citation_counts"],
        "tools_updated": stats["tools_updated"],
        "total_citations": stats["total_citations"],
        "processing_time": elapsed_time
    }
    
    with open(os.path.join(PUBMED_CITATIONS_DIR, "citation_counts_summary.json"), 'w') as f:
        json.dump(summary_data, f, indent=2)
    
    logger.info("Summary file updated")

if __name__ == "__main__":
    logger.info("Starting fix_citation_counts.py")
    fix_all_citation_counts()
    logger.info("Citation count fix completed")