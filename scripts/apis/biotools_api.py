#!/usr/bin/env python3
"""
API module for interacting with Bio.tools registry.
https://bio.tools/ is a comprehensive registry of bioinformatics tools and services.
"""

import requests
import json
import time
import logging
import os
from pathlib import Path
from datetime import datetime, timedelta

# Configure logging
logger = logging.getLogger(__name__)

# API Base URL
BIOTOOLS_API_BASE = "https://bio.tools/api"

# Cache directory for storing API responses
CACHE_DIR = Path(os.path.dirname(os.path.dirname(__file__))) / "cache" / "biotools"
CACHE_DIR.mkdir(exist_ok=True, parents=True)

# Cache expiration time (in days)
CACHE_EXPIRATION_DAYS = 30

class RateLimiter:
    """Rate limiter for Bio.tools API calls to avoid overloading the service."""
    def __init__(self, requests_per_minute=20):
        self.delay = 60.0 / requests_per_minute
        self.last_request_time = 0
        
    def wait(self):
        """Wait if needed to respect rate limits."""
        current_time = time.time()
        time_since_last_request = current_time - self.last_request_time
        
        if time_since_last_request < self.delay:
            sleep_time = self.delay - time_since_last_request
            time.sleep(sleep_time)
            
        self.last_request_time = time.time()

# Create a single rate limiter instance
rate_limiter = RateLimiter(requests_per_minute=20)

def get_cached_data(cache_key, expiration_days=CACHE_EXPIRATION_DAYS):
    """Get data from cache if available and not expired."""
    cache_file = CACHE_DIR / f"{cache_key}.json"
    
    if not cache_file.exists():
        return None
    
    try:
        file_modified_time = datetime.fromtimestamp(cache_file.stat().st_mtime)
        if datetime.now() - file_modified_time > timedelta(days=expiration_days):
            logger.debug(f"Cache for {cache_key} has expired")
            return None
            
        with open(cache_file, 'r', encoding='utf-8') as f:
            return json.load(f)
    except Exception as e:
        logger.warning(f"Error reading cache for {cache_key}: {e}")
        return None

def save_to_cache(cache_key, data):
    """Save data to cache."""
    cache_file = CACHE_DIR / f"{cache_key}.json"
    
    try:
        with open(cache_file, 'w', encoding='utf-8') as f:
            json.dump(data, f, indent=2)
        return True
    except Exception as e:
        logger.warning(f"Error saving cache for {cache_key}: {e}")
        return False

def search_tool(query, page=1, page_size=25, use_cache=True):
    """
    Search for tools in Bio.tools registry by name or keywords.
    
    Args:
        query (str): The search query
        page (int): Page number for pagination
        page_size (int): Number of results per page
        use_cache (bool): Whether to use cached results if available
        
    Returns:
        dict: Search results or None if the request failed
    """
    cache_key = f"search_{query.lower().replace(' ', '_')}_{page}_{page_size}"
    
    if use_cache:
        cached_data = get_cached_data(cache_key)
        if cached_data:
            logger.debug(f"Using cached data for Bio.tools search: {query}")
            return cached_data
    
    # Wait to respect rate limits
    rate_limiter.wait()
    
    params = {
        "q": query,
        "page": page,
        "format": "json",
        "page_size": page_size
    }
    
    try:
        logger.info(f"Searching Bio.tools for: {query}")
        response = requests.get(f"{BIOTOOLS_API_BASE}/tool", params=params)
        
        if response.status_code != 200:
            logger.warning(f"Bio.tools API returned status code {response.status_code} for query: {query}")
            return None
            
        data = response.json()
        
        # Cache the results
        save_to_cache(cache_key, data)
        
        return data
    except Exception as e:
        logger.error(f"Error searching Bio.tools API: {e}")
        return None

def get_tool_details(tool_id, use_cache=True):
    """
    Get detailed information about a specific tool from Bio.tools.
    
    Args:
        tool_id (str): The Bio.tools ID of the tool
        use_cache (bool): Whether to use cached results if available
        
    Returns:
        dict: Tool details or None if the request failed
    """
    cache_key = f"tool_{tool_id.lower()}"
    
    if use_cache:
        cached_data = get_cached_data(cache_key)
        if cached_data:
            logger.debug(f"Using cached data for Bio.tools tool: {tool_id}")
            return cached_data
    
    # Wait to respect rate limits
    rate_limiter.wait()
    
    try:
        logger.info(f"Fetching details from Bio.tools for tool: {tool_id}")
        response = requests.get(f"{BIOTOOLS_API_BASE}/tool/{tool_id}")
        
        if response.status_code != 200:
            logger.warning(f"Bio.tools API returned status code {response.status_code} for tool: {tool_id}")
            return None
            
        data = response.json()
        
        # Cache the results
        save_to_cache(cache_key, data)
        
        return data
    except Exception as e:
        logger.error(f"Error fetching tool details from Bio.tools API: {e}")
        return None

def get_tools_by_category(category, page=1, page_size=25, use_cache=True):
    """
    Get tools from Bio.tools by category (EDAM Operation or Topic).
    
    Args:
        category (str): The EDAM Operation or Topic term
        page (int): Page number for pagination
        page_size (int): Number of results per page
        use_cache (bool): Whether to use cached results if available
        
    Returns:
        dict: Search results or None if the request failed
    """
    cache_key = f"category_{category.lower().replace(' ', '_').replace(':', '_')}_{page}_{page_size}"
    
    if use_cache:
        cached_data = get_cached_data(cache_key)
        if cached_data:
            logger.debug(f"Using cached data for Bio.tools category: {category}")
            return cached_data
    
    # Wait to respect rate limits
    rate_limiter.wait()
    
    params = {
        "topic": category if category.startswith("topic") else None,
        "operation": category if category.startswith("operation") else None,
        "page": page,
        "format": "json",
        "page_size": page_size
    }
    
    # Remove None values
    params = {k: v for k, v in params.items() if v is not None}
    
    try:
        logger.info(f"Searching Bio.tools for category: {category}")
        response = requests.get(f"{BIOTOOLS_API_BASE}/tool", params=params)
        
        if response.status_code != 200:
            logger.warning(f"Bio.tools API returned status code {response.status_code} for category: {category}")
            return None
            
        data = response.json()
        
        # Cache the results
        save_to_cache(cache_key, data)
        
        return data
    except Exception as e:
        logger.error(f"Error searching Bio.tools API by category: {e}")
        return None

def extract_tool_metadata(tool_data):
    """
    Extract and organize relevant metadata from Bio.tools tool data.
    
    Args:
        tool_data (dict): The raw tool data from Bio.tools API
        
    Returns:
        dict: Structured metadata with relevant fields
    """
    if not tool_data:
        return None
        
    metadata = {
        "name": tool_data.get("name"),
        "biotoolsID": tool_data.get("biotoolsID"),
        "description": tool_data.get("description"),
        "homepage": tool_data.get("homepage"),
        "version": tool_data.get("version", [{}])[0].get("version") if tool_data.get("version") else None,
        "license": tool_data.get("license"),
        "publication_count": len(tool_data.get("publication", [])),
        "last_updated": tool_data.get("editPermission", {}).get("timeOfLastUpdate"),
        "topics": [topic.get("term") for topic in tool_data.get("topic", [])],
        "operations": [op.get("term") for op in tool_data.get("function", [{}])[0].get("operation", []) if tool_data.get("function")],
        "input_formats": [],
        "output_formats": [],
        "documentation": {url.get("url") for url in tool_data.get("documentation", []) if url.get("type") == "General"},
        "installation": {
            "conda": False,
            "bioconda": False,
            "pip": False,
            "container": False
        }
    }
    
    # Extract input/output formats
    for function in tool_data.get("function", []):
        for input_data in function.get("input", []):
            for format_data in input_data.get("format", []):
                metadata["input_formats"].append(format_data.get("term"))
                
        for output_data in function.get("output", []):
            for format_data in output_data.get("format", []):
                metadata["output_formats"].append(format_data.get("term"))
    
    # Check for installation methods
    for download in tool_data.get("download", []):
        if download.get("type") == "Container":
            metadata["installation"]["container"] = True
            
    # Check for Bioconda/Conda in links
    for link in tool_data.get("link", []):
        if "conda" in link.get("url", "").lower() or "anaconda" in link.get("url", "").lower():
            metadata["installation"]["conda"] = True
            if "bioconda" in link.get("url", "").lower():
                metadata["installation"]["bioconda"] = True
                
        if "pypi" in link.get("url", "").lower() or "pip" in link.get("url", "").lower():
            metadata["installation"]["pip"] = True
    
    # Remove duplicates from lists
    metadata["input_formats"] = list(set(metadata["input_formats"]))
    metadata["output_formats"] = list(set(metadata["output_formats"]))
    
    return metadata