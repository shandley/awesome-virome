#!/usr/bin/env python3
"""
Incremental Metadata Update Script

This script performs targeted updates to metadata for recently modified repositories
rather than updating all repositories. This approach reduces API load and improves efficiency.
"""

import os
import sys
import json
import logging
import argparse
import datetime
from pathlib import Path
from typing import Dict, List, Set, Any, Optional

# Try to import from the correct directory
script_dir = os.path.dirname(os.path.abspath(__file__))
root_dir = os.path.dirname(script_dir)
sys.path.append(script_dir)

# Import required modules
try:
    from apis.citations_api import cache_manager, GitHubAPI, SemanticScholarAPI, CrossRefAPI
    from enhance_metadata import extract_tool_metadata, get_repo_metadata, METADATA_DIR
except ImportError as e:
    print(f"Error importing required modules: {e}")
    print("Make sure you're running this script from the repository root.")
    sys.exit(1)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("metadata_enhancement.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Constants
METADATA_DIR = Path(os.path.join(root_dir, "metadata"))
DATA_JSON_PATH = os.path.join(root_dir, "data.json")
RECENT_DAYS_DEFAULT = 14  # Consider repos updated in last 14 days as "recent"

def load_data_json() -> Dict[str, Any]:
    """Load the current data.json file"""
    try:
        with open(DATA_JSON_PATH, 'r') as f:
            return json.load(f)
    except (json.JSONDecodeError, FileNotFoundError) as e:
        logger.error(f"Error loading data.json: {e}")
        return {"tools": []}

def get_recent_repos(data: Dict[str, Any], days: int = RECENT_DAYS_DEFAULT) -> List[Dict[str, Any]]:
    """
    Extract repositories that have been updated recently.
    
    Args:
        data: The current data.json content
        days: Number of days to consider as "recent"
        
    Returns:
        List of repository entries that were recently updated
    """
    cutoff_date = datetime.datetime.now() - datetime.timedelta(days=days)
    recent_repos = []
    
    for tool in data.get("tools", []):
        # Skip tools without repository info
        if not tool.get("repository"):
            continue
            
        # Check if the tool has update date information
        last_updated = tool.get("last_updated")
        if last_updated:
            try:
                update_date = datetime.datetime.fromisoformat(last_updated.replace('Z', '+00:00'))
                if update_date >= cutoff_date:
                    recent_repos.append(tool)
                    continue
            except (ValueError, TypeError):
                # If date parsing fails, fall back to next method
                pass
        
        # Check if tool has been recently added (no metadata file yet)
        tool_id = tool.get("id")
        if tool_id:
            metadata_path = METADATA_DIR / f"{tool_id}.json"
            if not metadata_path.exists():
                recent_repos.append(tool)
                
    logger.info(f"Found {len(recent_repos)} recently updated repositories out of {len(data.get('tools', []))} total")
    return recent_repos

def update_tool_metadata(tool: Dict[str, Any], github_token: Optional[str] = None) -> Dict[str, Any]:
    """
    Update metadata for a specific tool.
    
    Args:
        tool: Tool entry from data.json
        github_token: Optional GitHub API token
        
    Returns:
        Updated tool metadata
    """
    tool_id = tool.get("id")
    repo_url = tool.get("repository")
    
    if not tool_id or not repo_url:
        logger.warning(f"Skipping tool with missing ID or repository URL: {tool.get('name')}")
        return {}
    
    logger.info(f"Updating metadata for {tool_id}: {repo_url}")
    
    try:
        # Initialize API clients
        github_api = GitHubAPI(token=github_token)
        semantic_scholar_api = SemanticScholarAPI()
        crossref_api = CrossRefAPI()
        
        # Get repository metadata
        metadata = get_repo_metadata(repo_url, github_api)
        
        # Extract tool-specific metadata
        tool_metadata = extract_tool_metadata(tool_id, repo_url, tool.get("name", ""), 
                                              github_api, semantic_scholar_api, crossref_api)
        
        # Merge repository and tool metadata
        metadata.update(tool_metadata)
        
        # Add timestamp
        metadata["last_updated"] = datetime.datetime.now().isoformat()
        
        # Save to file
        metadata_path = METADATA_DIR / f"{tool_id}.json"
        with open(metadata_path, 'w') as f:
            json.dump(metadata, f, indent=2)
            
        logger.info(f"Successfully updated metadata for {tool_id}")
        return metadata
        
    except Exception as e:
        logger.error(f"Error updating metadata for {tool_id}: {e}")
        return {}

def main():
    """Main entry point for the script."""
    parser = argparse.ArgumentParser(description="Incremental metadata update for recently modified repositories")
    parser.add_argument("--days", type=int, default=RECENT_DAYS_DEFAULT, 
                        help=f"Number of days to consider as recent (default: {RECENT_DAYS_DEFAULT})")
    parser.add_argument("--recent-only", action="store_true", 
                        help="Only update repositories that have been recently modified")
    parser.add_argument("--limit", type=int, help="Limit the number of repositories to update")
    args = parser.parse_args()
    
    # Create metadata directory if it doesn't exist
    METADATA_DIR.mkdir(exist_ok=True)
    
    # Load current data
    data = load_data_json()
    
    # Get GitHub token from environment
    github_token = os.environ.get("GITHUB_TOKEN")
    if not github_token:
        logger.warning("No GitHub token found. API rate limits will be restricted.")
    
    # Get repositories to update
    if args.recent_only:
        repos_to_update = get_recent_repos(data, args.days)
    else:
        repos_to_update = data.get("tools", [])
    
    # Apply limit if specified
    if args.limit and args.limit > 0:
        repos_to_update = repos_to_update[:args.limit]
    
    logger.info(f"Starting incremental update for {len(repos_to_update)} repositories")
    
    # Update metadata for each repository
    updated_count = 0
    for tool in repos_to_update:
        metadata = update_tool_metadata(tool, github_token)
        if metadata:
            updated_count += 1
    
    logger.info(f"Incremental update completed. Updated {updated_count} repositories.")

if __name__ == "__main__":
    main()