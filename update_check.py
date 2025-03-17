#!/usr/bin/env python3
"""
Script to automatically check when GitHub/GitLab/Bitbucket repositories in the README were last updated.
Also identifies unavailable repositories (404 errors).
"""

import re
import os
import sys
import argparse
import time
import json
import logging
import requests
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Tuple, Optional, Dict, List, Any

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("update_check.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Define API base URLs
GITHUB_API = "https://api.github.com/repos/"
GITLAB_API = "https://gitlab.com/api/v4/projects/"
BITBUCKET_API = "https://api.bitbucket.org/2.0/repositories/"

# Maximum number of concurrent requests
MAX_WORKERS = 3

# Rate limiting delay (seconds) - increased to avoid GitHub rate limiting
RATE_LIMIT_DELAY = 2.0

# GitHub API token (set as environment variable GITHUB_TOKEN)
GITHUB_TOKEN = os.environ.get("GITHUB_TOKEN", "")

# Import cache manager if available
HAS_CACHE_MANAGER = False
try:
    # Add scripts directory to path to find the APIs module
    scripts_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "scripts")
    sys.path.append(scripts_dir)
    
    from apis.citations_api import cache_manager, RateLimiter as SmartRateLimiter
    HAS_CACHE_MANAGER = True
    logger.info("Smart cache system loaded. Repository information will be cached.")
    
    # Use the imported rate limiter
    USE_SMART_RATE_LIMITER = True
except ImportError as e:
    logger.warning(f"Smart cache system not available: {e}. Will use simple rate limiting.")
    USE_SMART_RATE_LIMITER = False


class RateLimiter:
    """Simple rate limiter for API requests."""
    def __init__(self, calls_per_minute=30, max_retries=3):
        self.delay = 60.0 / calls_per_minute
        self.last_request_time = 0
        self.max_retries = max_retries
        
    def wait(self):
        """Wait if needed to respect rate limits."""
        current_time = time.time()
        time_since_last_request = current_time - self.last_request_time
        
        if time_since_last_request < self.delay:
            sleep_time = self.delay - time_since_last_request
            time.sleep(sleep_time)
            
        self.last_request_time = time.time()
    
    def update_rate_limit(self, headers):
        """Placeholder for compatibility with smart rate limiter."""
        pass
    
    def handle_rate_limit_response(self, response):
        """Simple handling of rate limit errors."""
        if response.status_code == 429 or (response.status_code == 403 and 'rate limit' in response.text.lower()):
            logger.warning("Rate limit exceeded. Waiting 60 seconds before retry.")
            time.sleep(60)
            return True
        return False
        
    def reset_rate_limit(self):
        """Placeholder for compatibility with smart rate limiter."""
        pass


# Create rate limiters based on whether we have the smart implementation
if USE_SMART_RATE_LIMITER:
    github_limiter = SmartRateLimiter(calls_per_minute=30)  # GitHub allows more with token
    gitlab_limiter = SmartRateLimiter(calls_per_minute=30)
    bitbucket_limiter = SmartRateLimiter(calls_per_minute=30)
    logger.info("Using smart rate limiter with automatic rate limit detection")
else:
    github_limiter = RateLimiter(calls_per_minute=30) 
    gitlab_limiter = RateLimiter(calls_per_minute=30)
    bitbucket_limiter = RateLimiter(calls_per_minute=30)
    logger.info("Using simple rate limiter with fixed delays")

def extract_repos_from_readme(readme_path):
    """Extract repository URLs from the README file."""
    with open(readme_path, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Regular expression to match markdown links
    link_pattern = r'\[([^\]]+)\]\(([^)]+)\)'
    links = re.findall(link_pattern, content)
    
    repos = []
    for name, url in links:
        # Filter for GitHub, GitLab, and Bitbucket repositories
        if ('github.com/' in url or 'gitlab.com/' in url or 'bitbucket.org/' in url) and not url.endswith('.md'):
            repos.append((name, url.strip()))
    
    return repos

def get_github_repo_info(repo_url):
    """Get the last updated time, star count, and version info for a GitHub repository."""
    try:
        # Extract owner/repo from URL
        match = re.search(r'github\.com/([^/]+/[^/]+)', repo_url)
        if not match:
            return None, None, None, None
        
        repo_path = match.group(1).rstrip('/')
        # Handle branch or subdirectory in URL
        if '/tree/' in repo_path:
            repo_path = repo_path.split('/tree/')[0]
            
        # Check cache if available
        if HAS_CACHE_MANAGER:
            cache_key = f"repo_info_{repo_path}"
            cached_data = cache_manager.get(cache_key)
            if cached_data is not None:
                logger.info(f"Using cached GitHub repo data for {repo_path}")
                return cached_data
        
        # Use rate limiter
        github_limiter.wait()
        
        # Create API URL
        api_url = f"{GITHUB_API}{repo_path}"
        
        headers = {"Accept": "application/vnd.github.v3+json"}
        if GITHUB_TOKEN:
            headers["Authorization"] = f"token {GITHUB_TOKEN}"
        
        logger.info(f"Fetching GitHub repo: {repo_path}")
        response = requests.get(api_url, headers=headers)
        
        # Update rate limit information
        github_limiter.update_rate_limit(response.headers)
        
        # Check for rate limiting
        if response.status_code == 403 and 'X-RateLimit-Remaining' in response.headers:
            remaining = int(response.headers['X-RateLimit-Remaining'])
            reset_time = int(response.headers['X-RateLimit-Reset'])
            reset_datetime = datetime.fromtimestamp(reset_time)
            now = datetime.now()
            wait_time = (reset_datetime - now).total_seconds()
            
            if remaining == 0 and wait_time > 0:
                logger.warning(f"GitHub API rate limit exceeded. Reset at {reset_datetime}. Waiting {wait_time} seconds.")
                time.sleep(min(wait_time + 1, 3600))  # Wait until reset, but not more than an hour
                return get_github_repo_info(repo_url)  # Retry after waiting
        
        if response.status_code == 404:
            logger.warning(f"Repository not found: {repo_url}")
            result = (None, "not_found", None, None)
            # Cache the not found result
            if HAS_CACHE_MANAGER:
                cache_manager.set(cache_key, result, repo_url)
            return result
        
        if response.status_code == 429 or github_limiter.handle_rate_limit_response(response):
            logger.warning(f"Rate limit hit for {repo_url}, waiting and retrying")
            time.sleep(10)  # Wait 10 seconds before retrying
            return get_github_repo_info(repo_url)  # Retry
            
        response.raise_for_status()
        
        data = response.json()
        updated_at = data.get('updated_at') or data.get('pushed_at')
        stars = data.get('stargazers_count', 0)
        
        # Get version information from releases
        version_info = None
        if GITHUB_TOKEN:  # Only fetch releases if we have a token to avoid rate limits
            try:
                # Check cache for releases if available
                releases_cache_key = f"releases_{repo_path}"
                if HAS_CACHE_MANAGER:
                    cached_releases = cache_manager.get(releases_cache_key)
                    if cached_releases is not None:
                        logger.info(f"Using cached GitHub releases for {repo_path}")
                        releases = cached_releases
                    else:
                        github_limiter.wait()
                        releases_url = f"{GITHUB_API}{repo_path}/releases"
                        releases_response = requests.get(
                            releases_url, 
                            headers=headers, 
                            params={"per_page": 1}  # Only get the latest release
                        )
                        
                        github_limiter.update_rate_limit(releases_response.headers)
                        
                        if releases_response.status_code == 200:
                            releases = releases_response.json()
                            # Cache the releases data
                            if HAS_CACHE_MANAGER:
                                cache_manager.set(releases_cache_key, releases, repo_url)
                        else:
                            releases = []
                else:
                    github_limiter.wait()
                    releases_url = f"{GITHUB_API}{repo_path}/releases"
                    releases_response = requests.get(
                        releases_url, 
                        headers=headers, 
                        params={"per_page": 1}  # Only get the latest release
                    )
                    
                    github_limiter.update_rate_limit(releases_response.headers)
                    
                    if releases_response.status_code == 200:
                        releases = releases_response.json()
                    else:
                        releases = []
                
                if releases:
                    latest_release = releases[0]
                    tag_name = latest_release.get('tag_name', '')
                    # Clean up tag name to get version
                    version = tag_name
                    if version.startswith('v'):
                        version = version  # Keep the v prefix
                    elif version.startswith('release-'):
                        version = 'v' + version[8:]
                    
                    # Get release date
                    release_date = latest_release.get('published_at')
                    if release_date:
                        release_year = datetime.strptime(release_date, '%Y-%m-%dT%H:%M:%SZ').year
                        version_info = (version, release_year)
            except Exception as e:
                logger.warning(f"Error fetching releases for {repo_path}: {e}")
        
        if updated_at:
            result = (datetime.strptime(updated_at, '%Y-%m-%dT%H:%M:%SZ'), None, stars, version_info)
        else:
            result = (None, None, stars, version_info)
            
        # Cache the result if caching is available
        if HAS_CACHE_MANAGER:
            cache_manager.set(cache_key, result, repo_url)
            
        return result
    except requests.exceptions.RequestException as e:
        if hasattr(e, 'response') and e.response is not None:
            if e.response.status_code == 404:
                logger.warning(f"Repository not found: {repo_url}")
                result = (None, "not_found", None, None)
                # Cache the not found result
                if HAS_CACHE_MANAGER:
                    cache_key = f"repo_info_{repo_path if 'repo_path' in locals() else 'unknown'}"
                    cache_manager.set(cache_key, result, repo_url)
                return result
            elif e.response.status_code == 429:
                logger.warning(f"Rate limit hit for {repo_url}, waiting 60 seconds")
                time.sleep(60)  # Wait 60 seconds before retrying
                return get_github_repo_info(repo_url)  # Retry
        logger.error(f"Error fetching GitHub repo {repo_url}: {e}")
        return None, "error", None, None
    except Exception as e:
        logger.error(f"Unexpected error for {repo_url}: {e}")
        return None, "error", None, None

def get_gitlab_repo_info(repo_url):
    """Get the last updated time, star count, and version info for a GitLab repository."""
    try:
        # Extract project ID or path from URL
        match = re.search(r'gitlab\.com/([^/]+/[^/]+)', repo_url)
        if not match:
            return None, None, None, None
        
        project_path = match.group(1).rstrip('/')
        # Handle branch or subdirectory in URL
        if '/tree/' in project_path:
            project_path = project_path.split('/tree/')[0]
        
        # Check cache if available
        if HAS_CACHE_MANAGER:
            cache_key = f"gitlab_repo_info_{project_path}"
            cached_data = cache_manager.get(cache_key)
            if cached_data is not None:
                logger.info(f"Using cached GitLab repo data for {project_path}")
                return cached_data
            
        # Use rate limiter
        gitlab_limiter.wait()
        
        # URL encode the project path
        encoded_path = requests.utils.quote(project_path, safe='')
        
        # Create API URL for project lookup
        api_url = f"{GITLAB_API}?search={encoded_path}"
        
        logger.info(f"Fetching GitLab repo: {project_path}")
        response = requests.get(api_url)
        
        # Update rate limit headers
        gitlab_limiter.update_rate_limit(response.headers)
        
        if response.status_code == 404:
            logger.warning(f"Repository not found: {repo_url}")
            result = (None, "not_found", None, None)
            # Cache the not found result
            if HAS_CACHE_MANAGER:
                cache_manager.set(cache_key, result, repo_url)
            return result
        
        if response.status_code == 429 or gitlab_limiter.handle_rate_limit_response(response):
            logger.warning(f"Rate limit hit for {repo_url}, waiting 60 seconds")
            time.sleep(60)  # Wait 60 seconds before retrying
            return get_gitlab_repo_info(repo_url)  # Retry
            
        response.raise_for_status()
        
        projects = response.json()
        for project in projects:
            if project['path_with_namespace'].lower() == project_path.lower():
                project_id = project['id']
                
                # Check cache for project details if available
                project_details_cache_key = f"gitlab_project_{project_id}"
                project_data = None
                
                if HAS_CACHE_MANAGER:
                    project_data = cache_manager.get(project_details_cache_key)
                    if project_data is not None:
                        logger.info(f"Using cached GitLab project details for {project_path}")
                
                if project_data is None:
                    # Get detailed project info from API
                    project_url = f"{GITLAB_API}{project_id}"
                    gitlab_limiter.wait()
                    proj_response = requests.get(project_url)
                    gitlab_limiter.update_rate_limit(proj_response.headers)
                    proj_response.raise_for_status()
                    
                    project_data = proj_response.json()
                    
                    # Cache the project details
                    if HAS_CACHE_MANAGER:
                        cache_manager.set(project_details_cache_key, project_data, repo_url)
                
                updated_at = project_data.get('last_activity_at')
                stars = project_data.get('star_count', 0)
                
                # Get version information from tags
                version_info = None
                try:
                    # Check cache for tags if available
                    tags_cache_key = f"gitlab_tags_{project_id}"
                    tags = None
                    
                    if HAS_CACHE_MANAGER:
                        tags = cache_manager.get(tags_cache_key)
                        if tags is not None:
                            logger.info(f"Using cached GitLab tags for {project_path}")
                    
                    if tags is None:
                        gitlab_limiter.wait()
                        tags_url = f"{GITLAB_API}{project_id}/repository/tags"
                        tags_response = requests.get(
                            tags_url,
                            params={"per_page": 1, "order_by": "updated", "sort": "desc"}
                        )
                        gitlab_limiter.update_rate_limit(tags_response.headers)
                        
                        if tags_response.status_code == 200:
                            tags = tags_response.json()
                            # Cache the tags
                            if HAS_CACHE_MANAGER:
                                cache_manager.set(tags_cache_key, tags, repo_url)
                        else:
                            tags = []
                    
                    if tags:
                        latest_tag = tags[0]
                        tag_name = latest_tag.get('name', '')
                        # Clean up tag name to get version
                        version = tag_name
                        if version.startswith('v'):
                            version = version  # Keep the v prefix
                        elif version.startswith('release-'):
                            version = 'v' + version[8:]
                        
                        # Get release date (use commit date as GitLab doesn't have release dates)
                        if 'commit' in latest_tag and 'committed_date' in latest_tag['commit']:
                            commit_date = latest_tag['commit']['committed_date']
                            if commit_date:
                                commit_year = datetime.strptime(commit_date, '%Y-%m-%dT%H:%M:%S.%fZ').year
                                version_info = (version, commit_year)
                except Exception as e:
                    logger.warning(f"Error fetching tags for {project_path}: {e}")
                
                if updated_at:
                    result = (datetime.strptime(updated_at, '%Y-%m-%dT%H:%M:%S.%fZ'), None, stars, version_info)
                else:
                    result = (None, None, stars, version_info)
                
                # Cache the final result
                if HAS_CACHE_MANAGER:
                    cache_manager.set(cache_key, result, repo_url)
                
                return result
                
        # No matching project found
        result = (None, None, None, None)
        if HAS_CACHE_MANAGER:
            cache_manager.set(cache_key, result, repo_url)
        return result
    except requests.exceptions.RequestException as e:
        if hasattr(e, 'response') and e.response is not None:
            if e.response.status_code == 404:
                logger.warning(f"Repository not found: {repo_url}")
                result = (None, "not_found", None, None)
                # Cache the not found result if we have extracted the project path
                if HAS_CACHE_MANAGER and 'project_path' in locals():
                    cache_key = f"gitlab_repo_info_{project_path}"
                    cache_manager.set(cache_key, result, repo_url)
                return result
            elif e.response.status_code == 429:
                logger.warning(f"Rate limit hit for {repo_url}, waiting 60 seconds")
                time.sleep(60)  # Wait 60 seconds before retrying
                return get_gitlab_repo_info(repo_url)  # Retry
        logger.error(f"Error fetching GitLab repo {repo_url}: {e}")
        return None, "error", None, None
    except Exception as e:
        logger.error(f"Unexpected error for {repo_url}: {e}")
        return None, "error", None, None

def get_bitbucket_repo_info(repo_url):
    """Get the last updated time, watchers count, and version info for a Bitbucket repository."""
    try:
        # Extract workspace/repo from URL
        match = re.search(r'bitbucket\.org/([^/]+/[^/]+)', repo_url)
        if not match:
            return None, None, None, None
        
        repo_path = match.group(1).rstrip('/')
        # Handle branch or subdirectory in URL
        if '/src/' in repo_path:
            repo_path = repo_path.split('/src/')[0]
        
        # Check cache if available
        if HAS_CACHE_MANAGER:
            cache_key = f"bitbucket_repo_info_{repo_path}"
            cached_data = cache_manager.get(cache_key)
            if cached_data is not None:
                logger.info(f"Using cached Bitbucket repo data for {repo_path}")
                return cached_data
            
        # Use rate limiter
        bitbucket_limiter.wait()
        
        # Create API URL
        api_url = f"{BITBUCKET_API}{repo_path}"
        
        logger.info(f"Fetching Bitbucket repo: {repo_path}")
        response = requests.get(api_url)
        
        # Update rate limit information
        bitbucket_limiter.update_rate_limit(response.headers)
        
        if response.status_code == 404:
            logger.warning(f"Repository not found: {repo_url}")
            result = (None, "not_found", None, None)
            # Cache the not found result
            if HAS_CACHE_MANAGER:
                cache_manager.set(cache_key, result, repo_url)
            return result
        
        if response.status_code == 429 or bitbucket_limiter.handle_rate_limit_response(response):
            logger.warning(f"Rate limit hit for {repo_url}, waiting 60 seconds")
            time.sleep(60)  # Wait 60 seconds before retrying
            return get_bitbucket_repo_info(repo_url)  # Retry
            
        response.raise_for_status()
        
        data = response.json()
        updated_on = data.get('updated_on')
        # Bitbucket doesn't have stars but has watchers
        watchers = data.get('watchers_count', 0)
        
        # Get version information from tags
        version_info = None
        try:
            # Check cache for tags if available
            tags_cache_key = f"bitbucket_tags_{repo_path}"
            tags_data = None
            
            if HAS_CACHE_MANAGER:
                tags_data = cache_manager.get(tags_cache_key)
                if tags_data is not None:
                    logger.info(f"Using cached Bitbucket tags for {repo_path}")
            
            if tags_data is None:
                bitbucket_limiter.wait()
                tags_url = f"{BITBUCKET_API}{repo_path}/refs/tags?sort=-target.date&pagelen=1"
                tags_response = requests.get(tags_url)
                bitbucket_limiter.update_rate_limit(tags_response.headers)
                
                if tags_response.status_code == 200:
                    tags_data = tags_response.json()
                    # Cache the tags data
                    if HAS_CACHE_MANAGER:
                        cache_manager.set(tags_cache_key, tags_data, repo_url)
                else:
                    tags_data = {"values": []}
            
            if 'values' in tags_data and tags_data['values']:
                latest_tag = tags_data['values'][0]
                tag_name = latest_tag.get('name', '')
                # Clean up tag name to get version
                version = tag_name
                if version.startswith('v'):
                    version = version  # Keep the v prefix
                elif version.startswith('release-'):
                    version = 'v' + version[8:]
                
                # Get tag date
                if 'target' in latest_tag and 'date' in latest_tag['target']:
                    tag_date = latest_tag['target']['date']
                    if tag_date:
                        # Bitbucket date format might vary
                        try:
                            tag_year = datetime.strptime(tag_date, '%Y-%m-%dT%H:%M:%S.%f%z').year
                            version_info = (version, tag_year)
                        except ValueError:
                            try:
                                tag_year = datetime.strptime(tag_date, '%Y-%m-%dT%H:%M:%S%z').year
                                version_info = (version, tag_year)
                            except ValueError:
                                logger.warning(f"Couldn't parse date format for {repo_path} tag: {tag_date}")
        except Exception as e:
            logger.warning(f"Error fetching tags for {repo_path}: {e}")
        
        if updated_on:
            result = (datetime.strptime(updated_on, '%Y-%m-%dT%H:%M:%S.%f%z'), None, watchers, version_info)
        else:
            result = (None, None, watchers, version_info)
            
        # Cache the final result
        if HAS_CACHE_MANAGER:
            cache_manager.set(cache_key, result, repo_url)
            
        return result
    except requests.exceptions.RequestException as e:
        if hasattr(e, 'response') and e.response is not None:
            if e.response.status_code == 404:
                logger.warning(f"Repository not found: {repo_url}")
                result = (None, "not_found", None, None)
                # Cache the not found result if we have the repository path
                if HAS_CACHE_MANAGER and 'repo_path' in locals():
                    cache_key = f"bitbucket_repo_info_{repo_path}"
                    cache_manager.set(cache_key, result, repo_url)
                return result
            elif e.response.status_code == 429:
                logger.warning(f"Rate limit hit for {repo_url}, waiting 60 seconds")
                time.sleep(60)  # Wait 60 seconds before retrying
                return get_bitbucket_repo_info(repo_url)  # Retry
        logger.error(f"Error fetching Bitbucket repo {repo_url}: {e}")
        return None, "error", None, None
    except Exception as e:
        logger.error(f"Unexpected error for {repo_url}: {e}")
        return None, "error", None, None

def get_repo_last_updated(repo_name, repo_url):
    """Get the last updated time, stars, and version for a repository based on its URL."""
    # Handle different repository hosts
    if 'github.com' in repo_url:
        last_updated, status, stars, version_info = get_github_repo_info(repo_url)
        return repo_name, repo_url, last_updated, status, stars, version_info
    elif 'gitlab.com' in repo_url:
        last_updated, status, stars, version_info = get_gitlab_repo_info(repo_url)
        return repo_name, repo_url, last_updated, status, stars, version_info
    elif 'bitbucket.org' in repo_url:
        last_updated, status, stars, version_info = get_bitbucket_repo_info(repo_url)
        return repo_name, repo_url, last_updated, status, stars, version_info
    return repo_name, repo_url, None, None, None, None

def batch_process_repos(repos, batch_size=5, batch_delay=10, max_retries=3):
    """Process repositories in batches to avoid rate limiting."""
    results = []
    total_repos = len(repos)
    
    # Validate input
    if not isinstance(repos, list):
        logger.error(f"Expected repos to be a list, got {type(repos)}")
        return []
        
    if not all(isinstance(repo, tuple) and len(repo) == 2 for repo in repos):
        logger.warning("Some repos have invalid format. Expected a list of (name, url) tuples.")
        # Try to fix by filtering
        repos = [repo for repo in repos if isinstance(repo, tuple) and len(repo) == 2]
        if not repos:
            logger.error("No valid repositories to process")
            return []
    
    # Retry mechanism for batches that completely fail
    failed_batches = []
    retry_attempts = 0
    
    for i in range(0, total_repos, batch_size):
        batch = repos[i:i + batch_size]
        logger.info(f"Processing batch {i//batch_size + 1}/{(total_repos + batch_size - 1)//batch_size} ({len(batch)} repos)")
        
        batch_results = []
        error_count = 0
        
        with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
            futures = {executor.submit(get_repo_last_updated, name, url): (name, url) for name, url in batch}
            
            for future in as_completed(futures):
                name, url = futures[future]
                try:
                    result = future.result()
                    # Validate result
                    if not isinstance(result, tuple):
                        logger.error(f"Invalid result format for {name}: {result}")
                        result = (name, url, None, "error", None, None)
                        error_count += 1
                    elif len(result) not in (5, 6):
                        logger.error(f"Invalid result tuple length for {name}: {result}")
                        result = (name, url, None, "error", None, None)
                        error_count += 1
                        
                    batch_results.append(result)
                    logger.info(f"Processed {name} ({url})")
                except Exception as e:
                    logger.error(f"Error processing {name} ({url}): {e}")
                    batch_results.append((name, url, None, "error", None, None))
                    error_count += 1
        
        # Check if this batch completely failed
        if error_count == len(batch) and retry_attempts < max_retries:
            logger.warning(f"Batch {i//batch_size + 1} completely failed. Adding to retry queue.")
            failed_batches.append(batch)
        else:
            results.extend(batch_results)
        
        # If there are more batches to process, wait to avoid rate limiting
        if i + batch_size < total_repos:
            logger.info(f"Waiting {batch_delay} seconds before processing next batch...")
            time.sleep(batch_delay)
    
    # Retry failed batches
    if failed_batches and retry_attempts < max_retries:
        retry_attempts += 1
        logger.warning(f"Retrying {len(failed_batches)} failed batches (attempt {retry_attempts}/{max_retries})...")
        
        # Increase delay for retries to avoid rate limits
        retry_delay = batch_delay * 2
        
        # Flatten the list of batches
        retry_repos = [repo for batch in failed_batches for repo in batch]
        
        # Give some time before retrying
        time.sleep(retry_delay)
        
        # Recursively process the failed batches with a smaller batch size
        retry_results = batch_process_repos(
            retry_repos, 
            batch_size=max(1, batch_size // 2),  # Reduce batch size for retries
            batch_delay=retry_delay,
            max_retries=max_retries - 1
        )
        
        # Add retry results
        results.extend(retry_results)
    
    return results

def extract_repo_categories(readme_path):
    """Extract repositories and their categories from the README file."""
    with open(readme_path, 'r', encoding='utf-8') as f:
        content = f.read()

    # Define the main categories we want to track
    main_categories = {
        "Virus and Phage Identification": [
            "## Virus and Phage Identification",
            "### Metagenome Analysis",
            "### Integrated Viruses",
            "### RNA Virus Identification"
        ],
        "Host Prediction": ["## Host Prediction"],
        "Genome Analysis": [
            "## Genome Analysis",
            "### Genome Annotation",
            "### Genome Assembly",
            "### Genome Completeness",
            "### Genome Comparison",
            "### Gene Finding"
        ],
        "Taxonomy": ["## Taxonomy"]
    }
    
    # Map each section header to its main category
    section_to_category = {}
    for main_cat, sections in main_categories.items():
        for section in sections:
            section_to_category[section] = main_cat
    
    # Initialize repo categories
    repo_categories = {}
    current_section = None
    current_category = None
    
    # Regular expression to match markdown links
    link_pattern = r'\[([^\]]+)\]\(([^)]+)\)'
    
    # Process the README line by line
    for line in content.split('\n'):
        # Check if this line is a section header
        if line.startswith('#'):
            for section, category in section_to_category.items():
                if line.strip() == section:
                    current_section = section
                    current_category = category
                    break
        
        # If we're in a tracked section, look for repository links
        if current_category:
            matches = re.findall(link_pattern, line)
            for name, url in matches:
                # Consider only GitHub, GitLab, and Bitbucket repositories
                if ('github.com/' in url or 'gitlab.com/' in url or 'bitbucket.org/' in url) and not url.endswith('.md'):
                    if name not in repo_categories:
                        repo_categories[name] = current_category
    
    return repo_categories

def validate_repo_data(repo_data):
    """Validate repository data and filter out invalid entries."""
    valid_repo_data = []
    for data in repo_data:
        if not isinstance(data, tuple):
            logger.error(f"Invalid repo data format: {data} is not a tuple, skipping")
            continue
            
        if len(data) not in (5, 6):
            logger.error(f"Invalid repo data tuple length: {data} has {len(data)} elements (expected 5 or 6), skipping")
            continue
            
        # Verify that the first two elements are strings (name and URL)
        if len(data) >= 2:
            if not isinstance(data[0], str) or not isinstance(data[1], str):
                logger.error(f"Invalid repo data: name and URL must be strings: {data}")
                continue
                
            # Make sure strings are not empty
            if not data[0].strip() or not data[1].strip():
                logger.error(f"Invalid repo data: name or URL is empty: {data}")
                continue
        
        valid_repo_data.append(data)
    
    if len(valid_repo_data) < len(repo_data):
        logger.warning(f"Filtered out {len(repo_data) - len(valid_repo_data)} invalid repo data entries")
        
    return valid_repo_data

def safe_write_file(file_path, content, mode="w"):
    """Safely write content to a file with proper error handling."""
    try:
        # Create directory if it doesn't exist
        os.makedirs(os.path.dirname(file_path), exist_ok=True)
        
        with open(file_path, mode, encoding='utf-8') as f:
            f.write(content)
        return True
    except Exception as e:
        logger.error(f"Error writing to file {file_path}: {e}")
        return False

def update_readme_with_dates_status_and_stars(readme_path, repo_data):
    """Update the README with the last updated information, availability status, and create a popular packages section."""
    try:
        with open(readme_path, 'r', encoding='utf-8') as f:
            content = f.read()
    except Exception as e:
        logger.error(f"Error reading README file {readme_path}: {e}")
        return None, [], []
    
    updated_content = content
    unavailable_repos = []
    github_repos_with_stars = []
    
    # Extract repository categories
    repo_categories = extract_repo_categories(readme_path)
    
    # Validate repo_data before processing
    repo_data = validate_repo_data(repo_data)
    
    # Sort repositories by name for consistent processing
    try:
        sorted_repos = sorted(repo_data, key=lambda x: x[0])
    except Exception as e:
        logger.error(f"Error sorting repositories: {e}, proceeding with unsorted data")
        sorted_repos = repo_data
    
    for data in sorted_repos:
        try:
            # Handle both old and new format of data tuple
            if len(data) == 6:
                repo_name, repo_url, last_updated, status, stars, version_info = data
            else:
                repo_name, repo_url, last_updated, status, stars = data
                version_info = None
                
            # Validate essential data
            if not repo_name or not repo_url:
                logger.warning(f"Skipping entry with empty name or URL: {data}")
                continue
                
            # Collect GitHub repos with stars for popular packages section
            if 'github.com' in repo_url and stars is not None and stars > 0:
                category = repo_categories.get(repo_name, "Other")
                github_repos_with_stars.append((repo_name, repo_url, stars, category))
                
            # Check if repository is unavailable
            if status == "not_found":
                unavailable_repos.append((repo_name, repo_url))
                try:
                    # Check if the repo is already marked as unavailable
                    repo_pattern = re.escape(f"[{repo_name}]({repo_url})") + r".*?(\[unavailable\])?"
                    unavailable_match = re.search(repo_pattern, updated_content)
                    
                    if unavailable_match and unavailable_match.group(1):
                        # Already marked as unavailable, no need to change
                        pass
                    else:
                        # Mark as unavailable
                        # Replace all existing tags with the new tag
                        # First remove all existing tags
                        temp_content = re.sub(
                            re.escape(f"[{repo_name}]({repo_url})") + r"(?:\s+\[.*?\])*",
                            f"[{repo_name}]({repo_url})",
                            updated_content
                        )
                        # Then add the new tag
                        updated_content = temp_content.replace(
                            f"[{repo_name}]({repo_url})",
                            f"[{repo_name}]({repo_url}) [unavailable]"
                        )
                except Exception as e:
                    logger.error(f"Error marking unavailable repo {repo_name}: {e}")
            else:
                try:
                    # Get existing version tag from README if it exists
                    version_pattern = re.escape(f"[{repo_name}]({repo_url})") + r".*?(\[v[\d\.]+,\s*\d{4}\])"
                    existing_version_match = re.search(version_pattern, updated_content)
                    
                    # Keep track of tags to add
                    tags = []
                    
                    # Add update date if available
                    if last_updated:
                        try:
                            year = last_updated.year
                            month = last_updated.month
                            tags.append(f"[Updated: {month:02d}/{year}]")
                        except Exception as e:
                            logger.error(f"Error formatting update date for {repo_name}: {e}")
                    
                    # Add version info if available from API and not already in README
                    if version_info and not existing_version_match:
                        try:
                            version, year = version_info
                            tags.append(f"[{version}, {year}]")
                        except Exception as e:
                            logger.error(f"Error formatting version info for {repo_name}: {e}")
                    elif existing_version_match:
                        # Keep existing version tag
                        tags.append(existing_version_match.group(1))
                    
                    if tags:
                        try:
                            # First remove all existing tags
                            tag_removal_pattern = re.escape(f"[{repo_name}]({repo_url})") + r"(?:\s+\[.*?\])*"
                            entry_match = re.search(tag_removal_pattern, updated_content)
                            
                            if entry_match:
                                entry_without_tags = entry_match.group(0)
                                base_entry = f"[{repo_name}]({repo_url})"
                                
                                # Build new entry with all tags
                                new_entry = base_entry + " " + " ".join(tags)
                                
                                # Replace the entire entry
                                updated_content = updated_content.replace(entry_without_tags, new_entry)
                            else:
                                logger.warning(f"Could not find entry for {repo_name} in README when trying to add tags")
                        except Exception as e:
                            logger.error(f"Error updating README entry for {repo_name}: {e}")
                except Exception as e:
                    logger.error(f"Error processing version information for {repo_name}: {e}")
        except Exception as e:
            logger.error(f"Error processing repo data: {e} for {data}")
            continue
    
    # Create overall popular packages section
    if github_repos_with_stars:
        top_repos = sorted(github_repos_with_stars, key=lambda x: x[2], reverse=True)[:20]  # Top 20 repos by stars
        
        popular_packages_section = "## Popular Packages\n\nRanked by GitHub stars:\n\n"
        for i, (name, url, stars, _) in enumerate(top_repos, 1):
            popular_packages_section += f"{i}. [{name}]({url}) - ⭐ {stars} stars\n"
        popular_packages_section += "\n"
        
        # Check if the section already exists
        if "## Popular Packages" in updated_content:
            # Replace existing section
            pattern = r"## Popular Packages.*?(?=\n## |\Z)"
            updated_content = re.sub(pattern, popular_packages_section, updated_content, flags=re.DOTALL)
        else:
            # Add after the introduction
            updated_content = updated_content.replace(
                "## Getting Started", 
                f"{popular_packages_section}## Getting Started"
            )
    
    # Create "Top Packages by Category" section
    if github_repos_with_stars:
        # Group repositories by category
        repos_by_category = {}
        for name, url, stars, category in github_repos_with_stars:
            if category not in repos_by_category:
                repos_by_category[category] = []
            repos_by_category[category].append((name, url, stars))
        
        # Sort repositories within each category by stars
        for category in repos_by_category:
            repos_by_category[category] = sorted(repos_by_category[category], key=lambda x: x[2], reverse=True)
        
        # Create the section content
        top_packages_section = "## Top Packages by Category\n\nHere are the most starred packages in key categories:\n\n"
        
        # Include only specific categories of interest
        categories_of_interest = [
            "Virus and Phage Identification",
            "Host Prediction",
            "Genome Analysis",
            "Taxonomy"
        ]
        
        for category in categories_of_interest:
            if category in repos_by_category and repos_by_category[category]:
                top_packages_section += f"### {category}\n"
                # Take top 3 from each category
                for i, (name, url, stars) in enumerate(repos_by_category[category][:3], 1):
                    top_packages_section += f"{i}. [{name}]({url}) - ⭐ {stars} stars\n"
                top_packages_section += "\n"
        
        # Check if the section already exists
        if "## Top Packages by Category" in updated_content:
            # Replace existing section
            pattern = r"## Top Packages by Category.*?(?=\n## |\Z)"
            updated_content = re.sub(pattern, top_packages_section, updated_content, flags=re.DOTALL)
        else:
            # Add after Popular Packages section
            if "## Popular Packages" in updated_content:
                updated_content = updated_content.replace(
                    "## Getting Started", 
                    f"{top_packages_section}## Getting Started"
                )
            else:
                # Add before Getting Started if Popular Packages doesn't exist
                updated_content = updated_content.replace(
                    "## Getting Started", 
                    f"{top_packages_section}## Getting Started"
                )
    
    # Write the updated content back to the README
    if not safe_write_file(readme_path, updated_content):
        logger.error("Failed to update README file")
    
    # Create a file with the list of unavailable repositories
    unavailable_content = "# Unavailable Repositories\n\n"
    unavailable_content += "The following repositories were not found (404 errors):\n\n"
    for name, url in unavailable_repos:
        unavailable_content += f"- [{name}]({url})\n"
    
    if not safe_write_file(Path(__file__).parent / "unavailable_repos.md", unavailable_content):
        logger.error("Failed to write unavailable repositories file")
    
    # Create a file with the list of starred repositories
    starred_content = "# GitHub Repositories by Stars\n\n"
    
    if github_repos_with_stars:
        try:
            sorted_by_stars = sorted(github_repos_with_stars, key=lambda x: x[2], reverse=True)
            starred_content += "| Repository | Category | Stars |\n"
            starred_content += "|------------|----------|-------|\n"
            for name, url, stars, category in sorted_by_stars:
                starred_content += f"| [{name}]({url}) | {category} | {stars} |\n"
        except Exception as e:
            logger.error(f"Error sorting starred repositories: {e}")
            starred_content += "Error occurred while sorting repositories.\n"
    else:
        starred_content += "No GitHub repositories with stars found.\n"
    
    if not safe_write_file(Path(__file__).parent / "starred_repos.md", starred_content):
        logger.error("Failed to write starred repositories file")
    
    return updated_content, unavailable_repos, github_repos_with_stars

def invalidate_repo_cache(repo_url):
    """Invalidate all cached data for a specific repository."""
    if not HAS_CACHE_MANAGER:
        return 0
        
    try:
        # Use the cache manager to invalidate all caches related to this repo
        count = cache_manager.invalidate_repo_caches(repo_url)
        if count > 0:
            logger.info(f"Invalidated {count} cached items for {repo_url}")
        return count
    except Exception as e:
        logger.error(f"Error invalidating cache for {repo_url}: {e}")
        return 0

def main(readme):
    """Main function to check and update repository information."""
    try:
        start_time = datetime.now()
        logger.info(f"Starting repository update check at {start_time}")
        
        # Check if GITHUB_TOKEN is set
        if not GITHUB_TOKEN:
            logger.warning("GITHUB_TOKEN environment variable not set. Rate limiting will be more restrictive.")
        
        # Log cache status if available
        if HAS_CACHE_MANAGER:
            try:
                metrics = cache_manager.get_metrics()
                logger.info(f"Cache statistics: {metrics['cache_files']} cached files, " +
                           f"{metrics['hit_rate']*100:.1f}% hit rate, " +
                           f"{metrics['hits']} hits, {metrics['misses']} misses")
            except Exception as e:
                logger.warning(f"Error getting cache metrics: {e}")
                
        readme_path = os.path.join(Path(__file__).parent, readme)
        
        print(f"Checking and updating repository information in {readme_path}")
        if not os.path.exists(readme_path):
            logger.error(f"README file not found at {readme_path}")
            sys.exit(1)
        
        # Step 1: Extract repositories from README
        try:
            logger.info(f"Extracting repositories from {readme_path}...")
            repos = extract_repos_from_readme(readme_path)
            if not repos:
                logger.error("No repositories found in README. Aborting.")
                return False
            logger.info(f"Found {len(repos)} repository URLs in the README")
        except Exception as e:
            logger.error(f"Error extracting repositories from README: {e}")
            return False
        
        # Step 2: Extract repository categories
        try:
            repo_categories = extract_repo_categories(readme_path)
            logger.info(f"Categorized {len(repo_categories)} repositories into their respective sections")
        except Exception as e:
            logger.error(f"Error extracting repository categories: {e}")
            repo_categories = {}  # Continue with empty categories
        
        # Step 3: Process repositories in batches
        try:
            results = batch_process_repos(repos, batch_size=10, batch_delay=15)
            if not results:
                logger.error("No results returned from batch processing. Aborting.")
                return False
        except Exception as e:
            logger.error(f"Error processing repositories in batches: {e}")
            return False
        
        # Step 4: Validate results before processing
        results = validate_repo_data(results)
        if not results:
            logger.error("No valid results after validation. Aborting.")
            return False
        
        # Step 5: Calculate statistics
        try:
            updated_repos = [r for r in results if r[2] is not None]
            unavailable_repos = [r for r in results if r[3] == "not_found"]
            error_repos = [r for r in results if r[3] == "error"]
            github_repos_with_stars = [r for r in results if 'github.com' in r[1] and r[4] is not None and r[4] > 0]
            
            logger.info(f"Found update times for {len(updated_repos)}/{len(repos)} repositories")
            logger.info(f"Found {len(unavailable_repos)} unavailable repositories")
            logger.info(f"Encountered errors with {len(error_repos)} repositories")
            logger.info(f"Found {len(github_repos_with_stars)} GitHub repositories with stars")
            
            # If we're using the cache manager, invalidate caches for repositories that were updated recently
            if HAS_CACHE_MANAGER:
                # Load the previous scan results if available
                prev_results = {}
                if os.path.exists(Path(__file__).parent / "repo_updates.json"):
                    try:
                        with open(Path(__file__).parent / "repo_updates.json", 'r') as f:
                            prev_data = json.load(f)
                            for item in prev_data:
                                if "url" in item and "last_updated" in item:
                                    prev_results[item["url"]] = item["last_updated"]
                    except (json.JSONDecodeError, IOError) as e:
                        logger.warning(f"Error loading previous repository data: {e}")
                
                # Invalidate caches for repositories that have been updated
                invalidated_count = 0
                for name, url, last_updated, status, stars, *rest in results:
                    if last_updated and url in prev_results and prev_results[url]:
                        prev_updated = datetime.fromisoformat(prev_results[url])
                        # If the repo has been updated since our last check, invalidate its caches
                        if last_updated > prev_updated:
                            logger.info(f"Repository {name} updated since last check. Invalidating caches.")
                            invalidated_count += invalidate_repo_cache(url)
                
                if invalidated_count > 0:
                    logger.info(f"Invalidated {invalidated_count} cached items for updated repositories")
        except Exception as e:
            logger.error(f"Error calculating statistics: {e}")
            # Continue despite statistics error
        
        # Step 6: Group repositories by category
        try:
            repos_by_category = {}
            for data in results:
                try:
                    # Handle both old and new format of data tuple
                    if len(data) == 6:
                        name, url, _, _, stars, _ = data
                    else:
                        name, url, _, _, stars = data
                        
                    if not name or not url:
                        continue
                        
                    if name in repo_categories and 'github.com' in url and stars is not None and stars > 0:
                        category = repo_categories[name]
                        if category not in repos_by_category:
                            repos_by_category[category] = []
                        repos_by_category[category].append((name, url, stars))
                except Exception as e:
                    logger.error(f"Error processing data tuple {data}: {e}")
                    continue
        except Exception as e:
            logger.error(f"Error grouping repositories by category: {e}")
            repos_by_category = {}  # Continue with empty categories
        
        # Step 7: Log information about the top repositories by category
        try:
            logger.info("Top repositories by category:")
            categories_of_interest = [
                "Virus and Phage Identification",
                "Host Prediction",
                "Genome Analysis",
                "Taxonomy"
            ]
            
            for category in categories_of_interest:
                if category in repos_by_category and repos_by_category[category]:
                    try:
                        top_repos = sorted(repos_by_category[category], key=lambda x: x[2], reverse=True)[:3]
                        logger.info(f"  {category}:")
                        for name, _, stars in top_repos:
                            logger.info(f"    - {name}: {stars} stars")
                    except Exception as e:
                        logger.error(f"Error sorting repositories for category {category}: {e}")
        except Exception as e:
            logger.error(f"Error logging top repositories: {e}")
        
        # Step 8: Update the README with the dates, availability status, and stars
        try:
            updated_content, unavailable, repos_with_stars = update_readme_with_dates_status_and_stars(readme_path, results)
            if updated_content:
                logger.info(f"README updated with repository information")
                logger.info(f"Updated both Popular Packages and Top Packages by Category sections")
                logger.info(f"Created unavailable_repos.md with {len(unavailable)} unavailable repositories")
                logger.info(f"Created starred_repos.md with {len(repos_with_stars)} starred repositories")
            else:
                logger.error("Failed to update README with repository information")
        except Exception as e:
            logger.error(f"Error updating README: {e}")
            # Continue despite README update error
        
        # Step 9: Save results to JSON for later reference
        try:
            json_results = []
            for data in results:
                try:
                    # Handle both old and new format of data tuple
                    if len(data) == 6:
                        name, url, date, status, stars, version_info = data
                    else:
                        name, url, date, status, stars = data
                        version_info = None
                        
                    category = repo_categories.get(name, "Other")
                    result_dict = {
                        "name": name,
                        "url": url,
                        "category": category,
                        "last_updated": date.isoformat() if date else None,
                        "status": status,
                        "stars": stars
                    }
                    
                    # Add version info if available
                    if version_info:
                        try:
                            version, year = version_info
                            result_dict["version"] = version
                            result_dict["version_year"] = year
                        except Exception as e:
                            logger.error(f"Error adding version info for {name}: {e}")
                            
                    json_results.append(result_dict)
                except Exception as e:
                    logger.error(f"Error creating JSON for {data}: {e}")
                    continue
                    
            # Write the JSON to file
            json_content = json.dumps(json_results, indent=2)
            if safe_write_file(Path(__file__).parent / "repo_updates.json", json_content):
                logger.info("Results saved to repo_updates.json")
            else:
                logger.error("Failed to save results to repo_updates.json")
        except Exception as e:
            logger.error(f"Error saving results to JSON: {e}")
        
        # Log final cache statistics if available
        if HAS_CACHE_MANAGER:
            try:
                metrics = cache_manager.get_metrics()
                logger.info(f"Final cache statistics: {metrics['cache_files']} cached files, " +
                           f"{metrics['hit_rate']*100:.1f}% hit rate, " +
                           f"{metrics['hits']} hits, {metrics['misses']} misses, " +
                           f"{metrics['sets']} sets, {metrics['invalidations']} invalidations")
            except Exception as e:
                logger.warning(f"Error getting final cache metrics: {e}")
        
        end_time = datetime.now()
        logger.info(f"Repository update check completed at {end_time} (duration: {end_time - start_time})")
        return True
        
    except Exception as e:
        logger.error(f"Unhandled exception in main function: {e}")
        return False

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Check and update repository information in README")
    parser.add_argument("-r", "--readme", type=str, help="Path to the README file", default="README.md")
    parser.add_argument("--exit-on-error", action="store_true", help="Exit with non-zero status code if an error occurs")
    args = parser.parse_args()

    result = main(args.readme)
    
    if not result and args.exit_on_error:
        logger.error("Exiting with error due to --exit-on-error flag")
        sys.exit(1)
    elif not result:
        logger.warning("Script completed with errors, but continuing due to no --exit-on-error flag")
    else:
        logger.info("Script completed successfully")
