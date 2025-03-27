#!/usr/bin/env python3
"""
API module for interacting with Bioconda packages.
Bioconda is a channel for the conda package manager specializing in bioinformatics software.
"""

import requests
import json
import time
import logging
import os
import re
from pathlib import Path
from datetime import datetime, timedelta

# Configure logging
logger = logging.getLogger(__name__)

# API Base URLs
ANACONDA_API_BASE = "https://api.anaconda.org"
BIOCONDA_GITHUB_API = "https://api.github.com/repos/bioconda/bioconda-recipes"

# Cache directory for storing API responses
CACHE_DIR = Path(os.path.dirname(os.path.dirname(__file__))) / "cache" / "bioconda"
CACHE_DIR.mkdir(exist_ok=True, parents=True)

# Cache expiration time (in days)
CACHE_EXPIRATION_DAYS = 30

class RateLimiter:
    """Rate limiter for Bioconda API calls to avoid overloading the service."""
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

# Create rate limiter instances
anaconda_limiter = RateLimiter(requests_per_minute=20)
github_limiter = RateLimiter(requests_per_minute=10)  # GitHub has stricter rate limits

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

def search_package(package_name, use_cache=True):
    """
    Search for a package in the Bioconda channel.
    
    Args:
        package_name (str): The name of the package to search for
        use_cache (bool): Whether to use cached results if available
        
    Returns:
        dict: Package information or None if the request failed
    """
    # Sanitize package name for cache key
    sanitized_name = re.sub(r'[^\w\-\.]', '_', package_name.lower())
    cache_key = f"package_{sanitized_name}"
    
    if use_cache:
        cached_data = get_cached_data(cache_key)
        if cached_data:
            logger.debug(f"Using cached data for Bioconda package: {package_name}")
            return cached_data
    
    # Wait to respect rate limits
    anaconda_limiter.wait()
    
    try:
        logger.info(f"Searching Bioconda for package: {package_name}")
        response = requests.get(f"{ANACONDA_API_BASE}/package/bioconda/{package_name}")
        
        if response.status_code != 200:
            if response.status_code == 404:
                logger.warning(f"Package not found in Bioconda: {package_name}")
            else:
                logger.warning(f"Anaconda API returned status code {response.status_code} for package: {package_name}")
            return None
            
        data = response.json()
        
        # Cache the results
        save_to_cache(cache_key, data)
        
        return data
    except Exception as e:
        logger.error(f"Error searching Bioconda package: {e}")
        return None

def get_package_files(package_name, version=None, use_cache=True):
    """
    Get files for a specific version of a Bioconda package.
    
    Args:
        package_name (str): The name of the package
        version (str): The version of the package, or None for latest
        use_cache (bool): Whether to use cached results if available
        
    Returns:
        dict: Package files or None if the request failed
    """
    # Generate a cache key
    sanitized_name = re.sub(r'[^\w\-\.]', '_', package_name.lower())
    version_str = f"_{version.replace('.', '_')}" if version else ""
    cache_key = f"files_{sanitized_name}{version_str}"
    
    if use_cache:
        cached_data = get_cached_data(cache_key)
        if cached_data:
            logger.debug(f"Using cached data for Bioconda package files: {package_name} {version or 'latest'}")
            return cached_data
    
    # Wait to respect rate limits
    anaconda_limiter.wait()
    
    try:
        logger.info(f"Fetching files for Bioconda package: {package_name} {version or 'latest'}")
        url = f"{ANACONDA_API_BASE}/package/bioconda/{package_name}/files"
        if version:
            url += f"?version={version}"
            
        response = requests.get(url)
        
        if response.status_code != 200:
            logger.warning(f"Anaconda API returned status code {response.status_code} for files of package: {package_name}")
            return None
            
        data = response.json()
        
        # Cache the results
        save_to_cache(cache_key, data)
        
        return data
    except Exception as e:
        logger.error(f"Error fetching Bioconda package files: {e}")
        return None

def get_package_recipe(package_name, use_cache=True):
    """
    Try to get the recipe (meta.yaml) for a Bioconda package from GitHub.
    
    Args:
        package_name (str): The name of the package
        use_cache (bool): Whether to use cached results if available
        
    Returns:
        dict: Recipe information or None if not found
    """
    sanitized_name = re.sub(r'[^\w\-\.]', '_', package_name.lower())
    cache_key = f"recipe_{sanitized_name}"
    
    if use_cache:
        cached_data = get_cached_data(cache_key)
        if cached_data:
            logger.debug(f"Using cached data for Bioconda recipe: {package_name}")
            return cached_data
    
    # Wait to respect rate limits
    github_limiter.wait()
    
    recipe_path = f"recipes/{package_name}/meta.yaml"
    github_token = os.environ.get("GITHUB_TOKEN", "")
    headers = {}
    if github_token:
        headers["Authorization"] = f"token {github_token}"
    
    try:
        logger.info(f"Fetching Bioconda recipe for: {package_name}")
        # First, get the recipe content
        content_url = f"{BIOCONDA_GITHUB_API}/contents/{recipe_path}"
        response = requests.get(content_url, headers=headers)
        
        if response.status_code != 200:
            logger.warning(f"Could not find recipe for {package_name}, status code: {response.status_code}")
            return None
            
        content_data = response.json()
        download_url = content_data.get("download_url")
        
        if not download_url:
            logger.warning(f"No download URL for recipe: {package_name}")
            return None
            
        # Get the actual recipe content
        github_limiter.wait()
        recipe_response = requests.get(download_url, headers=headers)
        
        if recipe_response.status_code != 200:
            logger.warning(f"Failed to download recipe for {package_name}")
            return None
            
        recipe_content = recipe_response.text
        
        # Parse the recipe to extract dependencies
        dependencies = extract_dependencies_from_recipe(recipe_content)
        
        result = {
            "name": package_name,
            "recipe_path": recipe_path,
            "recipe_url": content_data.get("html_url"),
            "dependencies": dependencies,
            "raw_content": recipe_content
        }
        
        # Cache the results
        save_to_cache(cache_key, result)
        
        return result
    except Exception as e:
        logger.error(f"Error fetching Bioconda recipe: {e}")
        return None

def extract_dependencies_from_recipe(recipe_content):
    """
    Extract dependencies from a Bioconda recipe (meta.yaml).
    
    Args:
        recipe_content (str): The content of the meta.yaml file
        
    Returns:
        dict: Dictionary with different types of dependencies
    """
    if not recipe_content:
        return {}
        
    dependencies = {
        "build": [],
        "run": [],
        "host": []
    }
    
    # Simple regex extraction for dependencies sections
    # This is a basic implementation - a proper YAML parser would be more robust
    build_deps = re.search(r'requirements:.*?build:.*?\n(.*?)(?:\n\n|\n[a-z]+:)', recipe_content, re.DOTALL)
    run_deps = re.search(r'requirements:.*?run:.*?\n(.*?)(?:\n\n|\n[a-z]+:)', recipe_content, re.DOTALL)
    host_deps = re.search(r'requirements:.*?host:.*?\n(.*?)(?:\n\n|\n[a-z]+:)', recipe_content, re.DOTALL)
    
    # Extract build dependencies
    if build_deps:
        build_lines = build_deps.group(1).strip().split('\n')
        for line in build_lines:
            dep = line.strip().strip('-').strip()
            if dep and not dep.startswith('#'):
                # Remove version specifications
                dep_name = re.sub(r'\s*[<>=!].+$', '', dep)
                dependencies["build"].append(dep_name.strip())
    
    # Extract run dependencies
    if run_deps:
        run_lines = run_deps.group(1).strip().split('\n')
        for line in run_lines:
            dep = line.strip().strip('-').strip()
            if dep and not dep.startswith('#'):
                # Remove version specifications
                dep_name = re.sub(r'\s*[<>=!].+$', '', dep)
                dependencies["run"].append(dep_name.strip())
    
    # Extract host dependencies
    if host_deps:
        host_lines = host_deps.group(1).strip().split('\n')
        for line in host_lines:
            dep = line.strip().strip('-').strip()
            if dep and not dep.startswith('#'):
                # Remove version specifications
                dep_name = re.sub(r'\s*[<>=!].+$', '', dep)
                dependencies["host"].append(dep_name.strip())
    
    return dependencies

def extract_package_metadata(package_data, files_data=None, recipe_data=None):
    """
    Extract and organize relevant metadata from Bioconda package data.
    
    Args:
        package_data (dict): The raw package data from Anaconda API
        files_data (dict): Optional file data for the package
        recipe_data (dict): Optional recipe data for the package
        
    Returns:
        dict: Structured metadata with relevant fields
    """
    if not package_data:
        return None
        
    # Basic metadata
    metadata = {
        "name": package_data.get("name"),
        "summary": package_data.get("summary"),
        "description": package_data.get("description"),
        "license": package_data.get("license"),
        "latest_version": package_data.get("latest_version"),
        "versions": package_data.get("versions", []),
        "dev_url": package_data.get("dev_url"),
        "home": package_data.get("home"),
        "publication": package_data.get("doc_url"),
        "conda_platforms": package_data.get("conda_platforms", []),
        "installation": {
            "conda": True,
            "bioconda": True
        },
        "dependencies": {}
    }
    
    # Add files metadata if available
    if files_data:
        platforms = set()
        python_versions = set()
        
        for file_data in files_data:
            if "attrs" in file_data:
                attrs = file_data["attrs"]
                if "subdir" in attrs:
                    platforms.add(attrs["subdir"])
                if "depends" in attrs:
                    for dep in attrs["depends"]:
                        if dep.startswith("python "):
                            python_version = dep.split(" ")[1].split(".")[0]
                            python_versions.add(python_version)
        
        metadata["supported_platforms"] = list(platforms)
        metadata["python_versions"] = list(python_versions)
    
    # Add recipe dependencies if available
    if recipe_data and "dependencies" in recipe_data:
        metadata["dependencies"] = recipe_data["dependencies"]
    
    return metadata