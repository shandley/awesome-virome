#!/usr/bin/env python3
"""
Cache Warming Script

This script proactively warms the cache with important data before major updates,
reducing API load and improving performance during heavy operations.
"""

import os
import sys
import json
import logging
import argparse
import datetime
import random
from pathlib import Path
from typing import Dict, List, Set, Any, Optional

# Try to import from the correct directory
script_dir = os.path.dirname(os.path.abspath(__file__))
root_dir = os.path.dirname(script_dir)
sys.path.append(script_dir)

# Import required modules
try:
    from apis.citations_api import cache_manager, GitHubAPI
except ImportError as e:
    print(f"Error importing required modules: {e}")
    print("Make sure you're running this script from the repository root.")
    sys.exit(1)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("cache_warming.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Constants
DATA_JSON_PATH = os.path.join(root_dir, "data.json")
METADATA_DIR = Path(os.path.join(root_dir, "metadata"))
POPULAR_TOOLS_COUNT = 20  # Number of most popular tools to warm

def load_data_json() -> Dict[str, Any]:
    """Load the current data.json file"""
    try:
        with open(DATA_JSON_PATH, 'r') as f:
            return json.load(f)
    except (json.JSONDecodeError, FileNotFoundError) as e:
        logger.error(f"Error loading data.json: {e}")
        return {"tools": []}

def load_metadata_for_tool(tool_id: str) -> Dict[str, Any]:
    """Load the metadata file for a specific tool"""
    metadata_path = METADATA_DIR / f"{tool_id}.json"
    try:
        if os.path.exists(metadata_path):
            with open(metadata_path, 'r') as f:
                return json.load(f)
        return {}
    except (json.JSONDecodeError, FileNotFoundError) as e:
        logger.warning(f"Error loading metadata for {tool_id}: {e}")
        return {}

def get_popular_tools(data: Dict[str, Any], count: int = POPULAR_TOOLS_COUNT) -> List[Dict[str, Any]]:
    """
    Get the most popular tools based on stars.

    Args:
        data: The data.json content
        count: Number of popular tools to return

    Returns:
        List of popular tool entries
    """
    # Combine data.json info with metadata
    tools_with_metadata = []

    for tool in data.get("tools", []):
        tool_id = tool.get("id")
        if tool_id:
            metadata = load_metadata_for_tool(tool_id)
            combined = {**tool, **metadata}
            tools_with_metadata.append(combined)

    # Sort by stars (descending)
    sorted_tools = sorted(
        tools_with_metadata,
        key=lambda x: x.get("stars", 0),
        reverse=True
    )

    return sorted_tools[:count]

def warm_cache_for_tool(tool: Dict[str, Any], github_api: GitHubAPI) -> bool:
    """
    Warm the cache for a specific tool by fetching its repository contents.

    Args:
        tool: Tool entry with metadata
        github_api: GitHub API client

    Returns:
        bool: True if successful, False otherwise
    """
    tool_id = tool.get("id")
    repo_url = tool.get("repository")

    if not tool_id or not repo_url:
        logger.warning(f"Skipping tool with missing ID or repository URL: {tool.get('name')}")
        return False

    logger.info(f"Warming cache for {tool_id}: {repo_url}")

    try:
        # Extract owner/repo for GitHub repos
        if "github.com" in repo_url:
            parts = repo_url.rstrip("/").split("/")
            if len(parts) >= 2:
                owner = parts[-2]
                repo = parts[-1].replace(".git", "")

                # Fetch repository contents
                github_api.get_repo_contents(owner, repo, "", repo_url)

        logger.info(f"Successfully warmed cache for {tool_id}")
        return True
    except Exception as e:
        logger.error(f"Error warming cache for {tool_id}: {e}")
        return False

def warm_important_caches(github_token: Optional[str] = None) -> int:
    """
    Warm caches for the most popular repositories.

    Args:
        github_token: Optional GitHub API token

    Returns:
        Number of successfully warmed caches
    """
    # Load data.json
    data = load_data_json()

    # Get popular tools to warm
    popular_tools = get_popular_tools(data)
    tools_to_warm = {tool.get("id"): tool for tool in popular_tools if tool.get("id")}

    logger.info(f"Found {len(tools_to_warm)} popular tools to warm cache for")

    # Initialize API client
    github_api = GitHubAPI(token=github_token)

    # Warm caches
    success_count = 0
    for tool_id, tool in tools_to_warm.items():
        success = warm_cache_for_tool(tool, github_api)
        if success:
            success_count += 1

    logger.info(f"Successfully warmed cache for {success_count} out of {len(tools_to_warm)} tools")
    return success_count

def warm_random_sample(count: int = 10, github_token: Optional[str] = None) -> int:
    """
    Warm caches for a random sample of repositories.

    Args:
        count: Number of random repositories to warm
        github_token: Optional GitHub API token

    Returns:
        Number of successfully warmed caches
    """
    # Load data.json
    data = load_data_json()

    # Get all tools with repositories
    tools_with_repos = [tool for tool in data.get("tools", []) if tool.get("repository") and tool.get("id")]

    # Select random sample
    sample_size = min(count, len(tools_with_repos))
    random_sample = random.sample(tools_with_repos, sample_size)

    logger.info(f"Selected {sample_size} random tools to warm cache for")

    # Initialize API client
    github_api = GitHubAPI(token=github_token)

    # Warm caches
    success_count = 0
    for tool in random_sample:
        success = warm_cache_for_tool(tool, github_api)
        if success:
            success_count += 1

    logger.info(f"Successfully warmed cache for {success_count} out of {sample_size} random tools")
    return success_count

def main():
    """Main entry point for the script."""
    parser = argparse.ArgumentParser(description="Cache warming for important repositories")
    parser.add_argument("--warm-important", action="store_true",
                        help="Warm cache for the most popular repositories")
    parser.add_argument("--warm-random", action="store_true",
                        help="Warm cache for a random sample of repositories")
    parser.add_argument("--count", type=int, default=10,
                        help="Number of random repositories to warm (for --warm-random)")
    args = parser.parse_args()

    # Get GitHub token from environment
    github_token = os.environ.get("GITHUB_TOKEN")
    if not github_token:
        logger.warning("No GitHub token found. API rate limits will be restricted.")

    start_time = datetime.datetime.now()
    logger.info(f"Starting cache warming at {start_time.isoformat()}")

    # Perform cache warming based on args
    if args.warm_random:
        success_count = warm_random_sample(args.count, github_token)
    else:
        # Default to warming popular caches
        success_count = warm_important_caches(github_token)

    end_time = datetime.datetime.now()
    duration = (end_time - start_time).total_seconds()

    logger.info(f"Cache warming completed in {duration:.1f} seconds")
    logger.info(f"Successfully warmed {success_count} caches")

    # Get updated cache metrics
    metrics = cache_manager.get_metrics()
    logger.info(f"Current cache metrics: {metrics['cache_files']} entries, {metrics['hit_rate']*100:.1f}% hit rate")

    return 0

if __name__ == "__main__":
    sys.exit(main())
