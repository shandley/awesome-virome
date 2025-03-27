#!/usr/bin/env python3
"""
Module for integrating with various citation and academic impact APIs.
"""

import os
import json
import time
import requests
import logging
import base64
import re
from datetime import datetime, timedelta
from typing import Dict, List, Optional, Any, Tuple, Union, Set
from urllib.parse import quote, urlencode
from pathlib import Path
import hashlib

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class RateLimiter:
    """
    Advanced rate limiter with adaptive behavior based on remaining API limits.
    Implements exponential backoff on rate limit errors and automatically
    adjusts wait times based on rate limits.
    """
    
    def __init__(self, calls_per_minute: int = 10, max_retries: int = 3):
        self.default_rate = calls_per_minute
        self.calls_per_minute = calls_per_minute
        self.interval = 60 / calls_per_minute
        self.last_call = 0
        self.max_retries = max_retries
        
        # Track rate limit information
        self.rate_limit_remaining = None
        self.rate_limit_reset = None
        self.rate_limit_total = None
        self.last_rate_limit_check = 0
        self.rate_limit_adjusted = False
    
    def wait(self):
        """Wait if necessary to respect the rate limit."""
        current_time = time.time()
        elapsed = current_time - self.last_call
        
        # If we have rate limit information, adjust dynamically
        if self.rate_limit_remaining is not None and self.rate_limit_reset is not None:
            # If we're close to the limit, slow down dramatically
            if self.rate_limit_remaining < 10:
                time_to_reset = max(0, self.rate_limit_reset - current_time)
                if time_to_reset > 0:
                    # Calculate a conservative interval to make the remaining requests last
                    if self.rate_limit_remaining > 0:
                        new_interval = time_to_reset / self.rate_limit_remaining
                        # Don't wait more than 60 seconds
                        if new_interval < 60:
                            logger.info(f"Rate limit low ({self.rate_limit_remaining}/{self.rate_limit_total}). " +
                                       f"Adjusting interval to {new_interval:.2f}s until reset " +
                                       f"in {time_to_reset:.1f}s")
                            self.interval = new_interval
                            self.rate_limit_adjusted = True
                    else:
                        # No requests remaining, wait until reset (up to a minute)
                        wait_time = min(time_to_reset, 60)
                        logger.warning(f"No API requests remaining. Waiting {wait_time:.1f}s for reset.")
                        time.sleep(wait_time)
                        # Reset to default after waiting
                        self.reset_rate_limit()
                        return
        
        # Standard rate limiting
        if elapsed < self.interval:
            time.sleep(self.interval - elapsed)
        
        self.last_call = time.time()
    
    def update_rate_limit(self, headers: Dict[str, str]):
        """
        Update rate limit information from API response headers.
        
        Args:
            headers: Response headers containing rate limit information
        """
        # GitHub style rate limit headers
        if 'X-RateLimit-Remaining' in headers and 'X-RateLimit-Reset' in headers:
            try:
                self.rate_limit_remaining = int(headers['X-RateLimit-Remaining'])
                self.rate_limit_reset = int(headers['X-RateLimit-Reset'])
                if 'X-RateLimit-Limit' in headers:
                    self.rate_limit_total = int(headers['X-RateLimit-Limit'])
                self.last_rate_limit_check = time.time()
                
                # Log if we're getting low on requests
                if self.rate_limit_remaining < 20 and self.rate_limit_total:
                    logger.warning(f"API rate limit low: {self.rate_limit_remaining}/{self.rate_limit_total} " +
                                 f"requests remaining. Reset in {self.rate_limit_reset - time.time():.1f}s")
            except (ValueError, TypeError):
                pass
        
        # Zenodo style rate limit headers (Retry-After)
        elif 'Retry-After' in headers:
            try:
                retry_after = int(headers['Retry-After'])
                if retry_after > 0:
                    logger.warning(f"API rate limited. Retry after {retry_after}s")
                    self.rate_limit_remaining = 0
                    self.rate_limit_reset = time.time() + retry_after
                    self.last_rate_limit_check = time.time()
            except (ValueError, TypeError):
                pass
    
    def handle_rate_limit_response(self, response) -> bool:
        """
        Handle a response that indicates rate limiting.
        Returns True if the request should be retried after waiting.
        
        Args:
            response: The HTTP response object
        """
        self.update_rate_limit(response.headers)
        
        if response.status_code == 429 or (response.status_code == 403 and 'rate limit' in response.text.lower()):
            # Get retry time
            retry_after = None
            if 'Retry-After' in response.headers:
                try:
                    retry_after = int(response.headers['Retry-After'])
                except (ValueError, TypeError):
                    pass
            
            if retry_after is None and self.rate_limit_reset:
                retry_after = max(1, int(self.rate_limit_reset - time.time()))
            
            # Default to reasonable wait time
            if retry_after is None or retry_after <= 0:
                retry_after = 60
            
            # Cap to reasonable max
            retry_after = min(retry_after, 300)
            
            logger.warning(f"Rate limit exceeded. Waiting {retry_after}s before retry.")
            time.sleep(retry_after)
            
            # Reset our tracking
            self.reset_rate_limit()
            return True
            
        return False
    
    def reset_rate_limit(self):
        """Reset rate limiting to default values."""
        if self.rate_limit_adjusted:
            logger.info("Resetting rate limit parameters to default values")
            self.interval = 60 / self.default_rate
            self.rate_limit_adjusted = False
        
        # Clear rate limit information
        self.rate_limit_remaining = None
        self.rate_limit_reset = None


class CacheManager:
    """
    Advanced cache manager with smart invalidation capabilities.
    Manages dependencies between cached items and provides selective invalidation.
    Includes built-in instrumentation for tracking cache performance.
    """
    
    def __init__(self, cache_dir: str, expiry_days: int = 30):
        """
        Initialize the cache manager.
        
        Args:
            cache_dir: Base directory for cache storage
            expiry_days: Default expiration time in days
        """
        self.base_cache_dir = Path(cache_dir)
        self.base_cache_dir.mkdir(exist_ok=True, parents=True)
        
        # Create a directory for dependency tracking
        self.deps_dir = self.base_cache_dir / "_dependencies"
        self.deps_dir.mkdir(exist_ok=True)
        
        # Create a directory for performance metrics
        self.metrics_dir = self.base_cache_dir / "_metrics"
        self.metrics_dir.mkdir(exist_ok=True)
        
        self.default_expiry = timedelta(days=expiry_days)
        
        # Map of repo URLs to their dependency caches
        self._dependency_map: Dict[str, Set[str]] = {}
        self._load_dependency_maps()
        
        # Performance metrics
        self._metrics = {
            "hits": 0,
            "misses": 0,
            "invalidations": 0,
            "sets": 0,
            "start_time": datetime.now().isoformat()
        }
        self._load_metrics()
    
    def _load_dependency_maps(self):
        """Load existing dependency maps from disk."""
        for deps_file in self.deps_dir.glob("*.json"):
            try:
                with open(deps_file, 'r') as f:
                    repo_key = deps_file.stem
                    self._dependency_map[repo_key] = set(json.load(f))
            except (json.JSONDecodeError, IOError) as e:
                logger.warning(f"Failed to load dependency map {deps_file}: {e}")
    
    def _load_metrics(self):
        """Load existing performance metrics from disk."""
        metrics_file = self.metrics_dir / "cache_metrics.json"
        if metrics_file.exists():
            try:
                with open(metrics_file, 'r') as f:
                    stored_metrics = json.load(f)
                    # Update only the metrics we track, ignore others
                    for key in ['hits', 'misses', 'invalidations', 'sets']:
                        if key in stored_metrics:
                            self._metrics[key] = stored_metrics[key]
            except (json.JSONDecodeError, IOError) as e:
                logger.warning(f"Failed to load cache metrics: {e}")
    
    def _save_metrics(self):
        """Save performance metrics to disk."""
        metrics_file = self.metrics_dir / "cache_metrics.json"
        try:
            # Add calculated metrics
            current_metrics = self._metrics.copy()
            total_requests = current_metrics['hits'] + current_metrics['misses']
            current_metrics['hit_rate'] = (
                current_metrics['hits'] / total_requests if total_requests > 0 else 0
            )
            current_metrics['last_updated'] = datetime.now().isoformat()
            
            with open(metrics_file, 'w') as f:
                json.dump(current_metrics, f, indent=2)
        except IOError as e:
            logger.warning(f"Failed to save cache metrics: {e}")
    
    def _save_dependency_map(self, repo_id: str):
        """Save the dependency map for a repository to disk."""
        if repo_id not in self._dependency_map:
            return
            
        deps_file = self.deps_dir / f"{repo_id}.json"
        try:
            with open(deps_file, 'w') as f:
                json.dump(list(self._dependency_map[repo_id]), f)
        except IOError as e:
            logger.warning(f"Failed to save dependency map for {repo_id}: {e}")
    
    def _hash_key(self, key: str) -> str:
        """Create a filesystem-safe hash of a cache key."""
        return hashlib.md5(key.encode()).hexdigest()
    
    def _get_cache_path(self, cache_key: str) -> Path:
        """Get the path for a cached item."""
        hashed_key = self._hash_key(cache_key)
        return self.base_cache_dir / f"{hashed_key}.json"
    
    def register_dependency(self, repo_url: str, cache_key: str):
        """
        Register a dependency between a repository and a cache item.
        
        Args:
            repo_url: Repository URL that the cache depends on
            cache_key: The cache key to associate with this repository
        """
        # Create a stable identifier for the repository
        repo_id = self._hash_key(repo_url)
        
        # Add the cache key to this repo's dependencies
        if repo_id not in self._dependency_map:
            self._dependency_map[repo_id] = set()
        
        self._dependency_map[repo_id].add(cache_key)
        self._save_dependency_map(repo_id)
    
    def is_valid(self, cache_key: str, expiry: Optional[timedelta] = None) -> bool:
        """
        Check if a cache item is valid (exists and not expired).
        
        Args:
            cache_key: The cache key to check
            expiry: Optional custom expiration time
            
        Returns:
            bool: True if the cache is valid, False otherwise
        """
        cache_path = self._get_cache_path(cache_key)
        
        if not cache_path.exists():
            return False
        
        try:
            # Check if expired
            expiry = expiry or self.default_expiry
            file_modified_time = datetime.fromtimestamp(cache_path.stat().st_mtime)
            if datetime.now() - file_modified_time > expiry:
                logger.debug(f"Cache for {cache_key} has expired")
                return False
                
            # The cache exists and is not expired
            with open(cache_path, 'r') as f:
                cache_data = json.load(f)
                
            # Extra validation: check if the data has the expected structure
            if not isinstance(cache_data, dict) or 'cache_date' not in cache_data:
                logger.warning(f"Invalid cache structure for {cache_key}")
                return False
                
            return True
            
        except (json.JSONDecodeError, ValueError, KeyError, IOError) as e:
            logger.warning(f"Error validating cache for {cache_key}: {e}")
            return False
    
    def get(self, cache_key: str, default=None):
        """
        Get data from cache.
        
        Args:
            cache_key: The cache key to retrieve
            default: Default value to return if cache miss
            
        Returns:
            The cached data or the default value
        """
        cache_path = self._get_cache_path(cache_key)
        
        if not self.is_valid(cache_key):
            # Track cache miss
            self._metrics['misses'] += 1
            self._save_metrics()
            return default
        
        try:
            with open(cache_path, 'r') as f:
                # Track cache hit
                self._metrics['hits'] += 1
                self._save_metrics()
                return json.load(f)['data']
        except (json.JSONDecodeError, KeyError, IOError) as e:
            # Track cache miss
            self._metrics['misses'] += 1
            self._save_metrics()
            logger.warning(f"Error reading cache for {cache_key}: {e}")
            return default
    
    def set(self, cache_key: str, data: Any, repo_url: Optional[str] = None):
        """
        Save data to cache.
        
        Args:
            cache_key: The cache key to store
            data: The data to cache
            repo_url: Optional repository URL to register as a dependency
            
        Returns:
            bool: True if successful, False otherwise
        """
        cache_path = self._get_cache_path(cache_key)
        
        try:
            # Track cache set
            self._metrics['sets'] += 1
            
            # Add timing information to the cache
            cache_data = {
                'cache_date': datetime.now().isoformat(),
                'data': data
            }
            
            with open(cache_path, 'w') as f:
                json.dump(cache_data, f, indent=2)
            
            # Register dependency if a repo URL was provided
            if repo_url:
                self.register_dependency(repo_url, cache_key)
            
            # Save metrics
            self._save_metrics()
                
            return True
            
        except IOError as e:
            logger.warning(f"Error saving cache for {cache_key}: {e}")
            return False
    
    def invalidate_repo_caches(self, repo_url: str):
        """
        Invalidate all caches associated with a specific repository.
        
        Args:
            repo_url: Repository URL whose caches should be invalidated
            
        Returns:
            int: Number of cache items invalidated
        """
        repo_id = self._hash_key(repo_url)
        
        if repo_id not in self._dependency_map:
            logger.debug(f"No cached dependencies found for repository: {repo_url}")
            return 0
        
        invalidated_count = 0
        for cache_key in list(self._dependency_map[repo_id]):  # Create a copy of the set to safely iterate
            cache_path = self._get_cache_path(cache_key)
            if cache_path.exists():
                try:
                    cache_path.unlink()
                    invalidated_count += 1
                    logger.debug(f"Invalidated cache: {cache_key}")
                except IOError as e:
                    logger.warning(f"Failed to invalidate cache {cache_key}: {e}")
        
        # Track invalidations
        if invalidated_count > 0:
            self._metrics['invalidations'] += invalidated_count
            self._save_metrics()
            
        # Clear the dependency list after invalidation
        self._dependency_map[repo_id] = set()
        self._save_dependency_map(repo_id)
        
        logger.info(f"Invalidated {invalidated_count} cache items for repository: {repo_url}")
        return invalidated_count
    
    def invalidate(self, cache_key: str):
        """
        Invalidate a specific cache item.
        
        Args:
            cache_key: The cache key to invalidate
            
        Returns:
            bool: True if successful, False otherwise
        """
        cache_path = self._get_cache_path(cache_key)
        
        if not cache_path.exists():
            return True
            
        try:
            cache_path.unlink()
            
            # Track invalidation
            self._metrics['invalidations'] += 1
            self._save_metrics()
            
            # Remove this cache key from all dependency maps
            for repo_id, cache_keys in self._dependency_map.items():
                if cache_key in cache_keys:
                    cache_keys.remove(cache_key)
                    self._save_dependency_map(repo_id)
                    
            return True
            
        except IOError as e:
            logger.warning(f"Failed to invalidate cache {cache_key}: {e}")
            return False
    
    def clear_all(self):
        """
        Clear all cache data.
        
        Returns:
            int: Number of cache items cleared
        """
        count = 0
        for cache_file in self.base_cache_dir.glob("*.json"):
            # Skip special directories
            if cache_file.is_dir() or cache_file.name.startswith("_"):
                continue
                
            try:
                cache_file.unlink()
                count += 1
            except IOError as e:
                logger.warning(f"Failed to clear cache file {cache_file}: {e}")
        
        # Clear dependency maps
        for deps_file in self.deps_dir.glob("*.json"):
            try:
                deps_file.unlink()
            except IOError:
                pass
                
        self._dependency_map = {}
        
        # Track invalidations but preserve metrics history
        if count > 0:
            self._metrics['invalidations'] += count
            self._save_metrics()
        
        logger.info(f"Cleared {count} cache items")
        return count
        
    def get_metrics(self) -> Dict[str, Any]:
        """
        Get cache performance metrics.
        
        Returns:
            Dict containing cache performance statistics
        """
        # Calculate derived metrics
        metrics = self._metrics.copy()
        total_requests = metrics['hits'] + metrics['misses']
        
        # Add calculated metrics
        metrics['total_requests'] = total_requests
        metrics['hit_rate'] = (
            metrics['hits'] / total_requests if total_requests > 0 else 0
        )
        metrics['cache_efficiency'] = (
            (metrics['hits'] - metrics['invalidations']) / metrics['sets'] 
            if metrics['sets'] > 0 else 0
        )
        
        # Add file count information
        metrics['cache_files'] = len(list(self.base_cache_dir.glob("*.json")))
        metrics['dependency_maps'] = len(self._dependency_map)
        
        return metrics

# Create a global cache manager
cache_manager = CacheManager(os.path.join("metadata", "cache"))

class ZenodoAPI:
    """Integration with Zenodo for DOI information."""
    
    BASE_URL = "https://zenodo.org/api"
    CACHE_EXPIRY = 30  # days
    
    def __init__(self, token: Optional[str] = None):
        self.token = token
        self.rate_limiter = RateLimiter(calls_per_minute=20)
    
    def _get_cache_key(self, query: str) -> str:
        """Generate a standardized cache key for the query."""
        return f"zenodo_{base64.urlsafe_b64encode(query.encode()).decode()}"
    
    def search_records(self, query: str, repo_url: Optional[str] = None) -> List[Dict[str, Any]]:
        """
        Search Zenodo records for the given query.
        
        Args:
            query: The search query string
            repo_url: Optional repository URL to associate this cache with
        """
        cache_key = self._get_cache_key(f"search_{query}")
        
        # Check if we have a valid cached result
        cached_data = cache_manager.get(cache_key)
        if cached_data is not None:
            logger.info(f"Using cached Zenodo results for '{query}'")
            return cached_data
        
        self.rate_limiter.wait()
        url = f"{self.BASE_URL}/records"
        params = {
            'q': query,
            'sort': 'mostrecent',
            'size': 10
        }
        
        headers = {}
        if self.token:
            headers['Authorization'] = f"Bearer {self.token}"
        
        try:
            response = requests.get(url, params=params, headers=headers)
            
            # Update rate limit information from headers
            self.rate_limiter.update_rate_limit(response.headers)
            
            response.raise_for_status()
            results = response.json().get('hits', {}).get('hits', [])
            
            # Cache the results and associate with repo_url if provided
            cache_manager.set(cache_key, results, repo_url)
            
            return results
        except requests.RequestException as e:
            logger.error(f"Error searching Zenodo: {e}")
            return []
    
    def get_doi_metadata(self, doi: str, repo_url: Optional[str] = None) -> Dict[str, Any]:
        """
        Get metadata for a specific DOI.
        
        Args:
            doi: The DOI to look up
            repo_url: Optional repository URL to associate this cache with
        """
        cache_key = self._get_cache_key(f"doi_{doi}")
        
        # Check if we have a valid cached result
        cached_data = cache_manager.get(cache_key)
        if cached_data is not None:
            logger.info(f"Using cached Zenodo DOI metadata for '{doi}'")
            return cached_data
        
        self.rate_limiter.wait()
        url = f"{self.BASE_URL}/records"
        params = {'q': f"doi:\"{doi}\""} 
        
        headers = {}
        if self.token:
            headers['Authorization'] = f"Bearer {self.token}"
        
        try:
            response = requests.get(url, params=params, headers=headers)
            
            # Update rate limit information from headers
            self.rate_limiter.update_rate_limit(response.headers)
            
            response.raise_for_status()
            results = response.json().get('hits', {}).get('hits', [])
            
            if results:
                # Cache the results and associate with repo_url if provided
                cache_manager.set(cache_key, results[0], repo_url)
                return results[0]
            return {}
        except requests.RequestException as e:
            logger.error(f"Error getting DOI metadata from Zenodo: {e}")
            return {}
    
    def find_tool_doi(self, tool_name: str, repo_url: str) -> Optional[str]:
        """Find a DOI for a tool based on name and repository URL."""
        search_terms = [
            tool_name,
            tool_name.replace('-', ' '),
            f"{tool_name} software",
            f"{tool_name} bioinformatics"
        ]
        
        for term in search_terms:
            results = self.search_records(term, repo_url)
            for record in results:
                metadata = record.get('metadata', {})
                # Check if this record matches our tool
                description = metadata.get('description', '').lower()
                related_ids = metadata.get('related_identifiers', [])
                
                # Look for GitHub URL in related identifiers
                for related in related_ids:
                    if related.get('relation', '') == 'isSupplementTo':
                        identifier = related.get('identifier', '')
                        if repo_url.lower() in identifier.lower():
                            return metadata.get('doi')
                
                # Check if tool name is in title
                title = metadata.get('title', '').lower()
                if tool_name.lower() in title and (
                    'software' in title or 
                    'tool' in title or 
                    'bioinformatics' in title or
                    'pipeline' in title or
                    'workflow' in title
                ):
                    return metadata.get('doi')
        
        return None


class SemanticScholarAPI:
    """Integration with Semantic Scholar API for citation metrics."""
    
    BASE_URL = "https://api.semanticscholar.org/graph/v1"
    CACHE_EXPIRY = 30  # days
    
    def __init__(self, api_key: Optional[str] = None):
        self.api_key = api_key
        self.rate_limiter = RateLimiter(calls_per_minute=10)
    
    def _get_cache_key(self, query: str) -> str:
        """Generate a standardized cache key for the query."""
        return f"semanticscholar_{base64.urlsafe_b64encode(query.encode()).decode()}"
    
    def search_paper(self, query: str, repo_url: Optional[str] = None) -> List[Dict[str, Any]]:
        """
        Search for papers matching a query.
        
        Args:
            query: The search query
            repo_url: Optional repository URL to associate this cache with
        """
        cache_key = self._get_cache_key(f"search_{query}")
        
        # Check if we have a valid cached result
        cached_data = cache_manager.get(cache_key)
        if cached_data is not None:
            logger.info(f"Using cached Semantic Scholar results for '{query}'")
            return cached_data
        
        self.rate_limiter.wait()
        url = f"{self.BASE_URL}/paper/search"
        params = {
            'query': query,
            'limit': 10,
            'fields': 'title,authors,year,venue,citationCount,externalIds,url,abstract'
        }
        
        headers = {}
        if self.api_key:
            headers['x-api-key'] = self.api_key
        
        try:
            response = requests.get(url, params=params, headers=headers)
            
            # Update rate limit information from headers
            self.rate_limiter.update_rate_limit(response.headers)
            
            response.raise_for_status()
            results = response.json().get('data', [])
            
            # Cache the results and associate with repo_url if provided
            cache_manager.set(cache_key, results, repo_url)
            
            return results
        except requests.RequestException as e:
            logger.error(f"Error searching Semantic Scholar: {e}")
            return []
    
    def get_paper_by_doi(self, doi: str, repo_url: Optional[str] = None) -> Dict[str, Any]:
        """
        Get paper information using DOI.
        
        Args:
            doi: The DOI to look up
            repo_url: Optional repository URL to associate this cache with
        """
        cache_key = self._get_cache_key(f"doi_{doi}")
        
        # Check if we have a valid cached result
        cached_data = cache_manager.get(cache_key)
        if cached_data is not None:
            logger.info(f"Using cached Semantic Scholar DOI data for '{doi}'")
            return cached_data
        
        self.rate_limiter.wait()
        url = f"{self.BASE_URL}/paper/DOI:{doi}"
        params = {
            'fields': 'title,authors,year,venue,citationCount,referenceCount,influentialCitationCount,externalIds,url,abstract,references,citations'
        }
        
        headers = {}
        if self.api_key:
            headers['x-api-key'] = self.api_key
        
        try:
            response = requests.get(url, params=params, headers=headers)
            response.raise_for_status()
            data = response.json()
            
            # Cache the results and associate with repo_url if provided
            cache_manager.set(cache_key, data, repo_url)
            
            return data
        except requests.RequestException as e:
            logger.error(f"Error getting paper by DOI from Semantic Scholar: {e}")
            return {}
    
    def get_paper_by_title_author(self, title: str, authors: List[str], repo_url: Optional[str] = None) -> Dict[str, Any]:
        """
        Find a paper by its title and authors.
        
        Args:
            title: The paper title
            authors: List of author names
            repo_url: Optional repository URL to associate this cache with
        """
        # Generate a search query with title and first author
        first_author = authors[0] if authors else ""
        query = f"{title} {first_author}"
        
        cache_key = self._get_cache_key(f"title_author_{base64.urlsafe_b64encode(query.encode()).decode()}")
        
        # Check if we have a valid cached result
        cached_data = cache_manager.get(cache_key)
        if cached_data is not None:
            logger.info(f"Using cached Semantic Scholar title/author data for '{title}'")
            return cached_data
        
        results = self.search_paper(query, repo_url)
        
        # Try to find an exact title match
        for paper in results:
            if paper.get('title', '').lower() == title.lower():
                # Further validate with author matching
                paper_authors = [author.get('name', '') for author in paper.get('authors', [])]
                # Check if at least one author matches
                for author in authors:
                    if any(author.lower() in paper_author.lower() for paper_author in paper_authors):
                        # Cache the results and associate with repo_url if provided
                        cache_manager.set(cache_key, paper, repo_url)
                        return paper
        
        # If no exact match, return the top result if it seems relevant
        if results and title.lower() in results[0].get('title', '').lower():
            # Cache the results and associate with repo_url if provided
            cache_manager.set(cache_key, results[0], repo_url)
            return results[0]
        
        # Cache the empty result to avoid repeated searches
        cache_manager.set(cache_key, {}, repo_url)
        return {}
    
    def get_citation_metrics(self, paper_id: str, repo_url: Optional[str] = None) -> Dict[str, Any]:
        """
        Get detailed citation metrics for a paper.
        
        Args:
            paper_id: Semantic Scholar paper ID
            repo_url: Optional repository URL to associate this cache with
        """
        cache_key = self._get_cache_key(f"metrics_{paper_id}")
        
        # Check if we have a valid cached result
        cached_data = cache_manager.get(cache_key)
        if cached_data is not None:
            logger.info(f"Using cached Semantic Scholar metrics for '{paper_id}'")
            return cached_data
        
        self.rate_limiter.wait()
        url = f"{self.BASE_URL}/paper/{paper_id}"
        params = {
            'fields': 'citationCount,influentialCitationCount,citations.year'
        }
        
        headers = {}
        if self.api_key:
            headers['x-api-key'] = self.api_key
        
        try:
            response = requests.get(url, params=params, headers=headers)
            response.raise_for_status()
            data = response.json()
            
            # Process citation data to get citations by year
            citations_by_year = {}
            for citation in data.get('citations', []):
                year = citation.get('year')
                if year:
                    citations_by_year[str(year)] = citations_by_year.get(str(year), 0) + 1
            
            metrics = {
                'total_citations': data.get('citationCount', 0),
                'influential_citations': data.get('influentialCitationCount', 0),
                'citations_by_year': citations_by_year
            }
            
            # Cache the results and associate with repo_url if provided
            cache_manager.set(cache_key, metrics, repo_url)
            
            return metrics
        except requests.RequestException as e:
            logger.error(f"Error getting citation metrics from Semantic Scholar: {e}")
            default_metrics = {
                'total_citations': 0,
                'influential_citations': 0,
                'citations_by_year': {}
            }
            # Cache the default metrics to avoid repeated failures
            cache_manager.set(cache_key, default_metrics, repo_url)
            return default_metrics
    
    def find_related_papers(self, paper_id: str, field: str = "virology", repo_url: Optional[str] = None) -> List[Dict[str, Any]]:
        """
        Find related papers in a specific field.
        
        Args:
            paper_id: Semantic Scholar paper ID
            field: Field to filter related papers by
            repo_url: Optional repository URL to associate this cache with
        """
        cache_key = self._get_cache_key(f"related_{paper_id}_{field}")
        
        # Check if we have a valid cached result
        cached_data = cache_manager.get(cache_key)
        if cached_data is not None:
            logger.info(f"Using cached related papers for '{paper_id}' in '{field}'")
            return cached_data
        
        self.rate_limiter.wait()
        url = f"{self.BASE_URL}/paper/{paper_id}/related"
        params = {
            'fields': 'title,authors,year,venue,citationCount,url,abstract',
            'limit': 20
        }
        
        headers = {}
        if self.api_key:
            headers['x-api-key'] = self.api_key
        
        try:
            response = requests.get(url, params=params, headers=headers)
            response.raise_for_status()
            all_related = response.json().get('data', [])
            
            # Filter for papers related to the specified field
            field_related = []
            for paper in all_related:
                abstract = paper.get('abstract', '').lower()
                title = paper.get('title', '').lower()
                venue = paper.get('venue', '').lower()
                
                if (field.lower() in abstract or
                    field.lower() in title or
                    'virus' in abstract or
                    'viral' in abstract or
                    'pathogen' in abstract):
                    field_related.append(paper)
            
            # Cache the results and associate with repo_url if provided
            result = field_related[:10]  # Store top 10
            cache_manager.set(cache_key, result, repo_url)
            
            return result
        except requests.RequestException as e:
            logger.error(f"Error finding related papers from Semantic Scholar: {e}")
            return []


class CrossRefAPI:
    """Integration with CrossRef API for publication information."""
    
    BASE_URL = "https://api.crossref.org"
    CACHE_EXPIRY = 30  # days
    
    def __init__(self, email: Optional[str] = None):
        """Initialize CrossRef API client with optional email for Polite Pool."""
        self.email = email
        self.rate_limiter = RateLimiter(calls_per_minute=10)
    
    def _get_cache_key(self, query: str) -> str:
        """Generate a standardized cache key for the query."""
        return f"crossref_{base64.urlsafe_b64encode(query.encode()).decode()}"
    
    def search_works(self, query: str, repo_url: Optional[str] = None) -> List[Dict[str, Any]]:
        """
        Search CrossRef works with the given query.
        
        Args:
            query: The search query
            repo_url: Optional repository URL to associate this cache with
        """
        cache_key = self._get_cache_key(f"search_{query}")
        
        # Check if we have a valid cached result
        cached_data = cache_manager.get(cache_key)
        if cached_data is not None:
            logger.info(f"Using cached CrossRef results for '{query}'")
            return cached_data
        
        self.rate_limiter.wait()
        url = f"{self.BASE_URL}/works"
        params = {
            'query': query,
            'rows': 10,
            'sort': 'relevance'
        }
        
        headers = {}
        if self.email:
            headers['User-Agent'] = f"AwesomeVirome/1.0 (mailto:{self.email})"
        
        try:
            response = requests.get(url, params=params, headers=headers)
            
            # Update rate limit information from headers
            self.rate_limiter.update_rate_limit(response.headers)
            
            response.raise_for_status()
            results = response.json().get('message', {}).get('items', [])
            
            # Cache the results and associate with repo_url if provided
            cache_manager.set(cache_key, results, repo_url)
            
            return results
        except requests.RequestException as e:
            logger.error(f"Error searching CrossRef: {e}")
            return []
    
    def get_work_by_doi(self, doi: str, repo_url: Optional[str] = None) -> Dict[str, Any]:
        """
        Get detailed work information by DOI.
        
        Args:
            doi: DOI to look up
            repo_url: Optional repository URL to associate this cache with
        """
        cache_key = self._get_cache_key(f"doi_{doi}")
        
        # Check if we have a valid cached result
        cached_data = cache_manager.get(cache_key)
        if cached_data is not None:
            logger.info(f"Using cached CrossRef DOI data for '{doi}'")
            return cached_data
        
        self.rate_limiter.wait()
        url = f"{self.BASE_URL}/works/{quote(doi)}"
        
        headers = {}
        if self.email:
            headers['User-Agent'] = f"AwesomeVirome/1.0 (mailto:{self.email})"
        
        try:
            response = requests.get(url, headers=headers)
            response.raise_for_status()
            data = response.json().get('message', {})
            
            # Cache the results and associate with repo_url if provided
            cache_manager.set(cache_key, data, repo_url)
            
            return data
        except requests.RequestException as e:
            logger.error(f"Error getting work by DOI from CrossRef: {e}")
            return {}
    
    def format_citation(self, work: Dict[str, Any], style: str = "apa") -> str:
        """Format a citation in the specified style."""
        if not work:
            return ""
        
        try:
            # Extract publication data
            title = work.get('title', [''])[0]
            authors = []
            for author in work.get('author', []):
                given = author.get('given', '')
                family = author.get('family', '')
                if given and family:
                    authors.append(f"{family}, {given[0]}.")
                elif family:
                    authors.append(family)
            
            year = work.get('published', {}).get('date-parts', [['']])[0][0]
            journal = work.get('container-title', [''])[0]
            volume = work.get('volume', '')
            issue = work.get('issue', '')
            pages = work.get('page', '')
            doi = work.get('DOI', '')
            
            if style.lower() == "apa":
                # Format authors
                if len(authors) == 0:
                    author_text = ""
                elif len(authors) == 1:
                    author_text = f"{authors[0]}"
                elif len(authors) < 8:
                    author_text = ", ".join(authors[:-1]) + f", & {authors[-1]}"
                else:
                    author_text = ", ".join(authors[:6]) + f", ... {authors[-1]}"
                
                # Format APA style citation
                citation = f"{author_text} ({year}). {title}."
                if journal:
                    citation += f" {journal}"
                    if volume:
                        citation += f", {volume}"
                        if issue:
                            citation += f"({issue})"
                    if pages:
                        citation += f", {pages}"
                citation += f". https://doi.org/{doi}" if doi else ""
                
                return citation
            
            elif style.lower() == "bibtex":
                # Create a BibTeX key from first author's last name and year
                first_author = work.get('author', [{}])[0].get('family', 'Unknown').lower()
                bibtex_key = f"{first_author}{year}"
                
                # Format BibTeX entry
                bibtex = f"@article{{{bibtex_key},\n"
                bibtex += f"  title = {{{title}}},\n"
                
                # Authors
                if authors:
                    bibtex += f"  author = {{{' and '.join(authors)}}},\n"
                
                # Other fields
                if year:
                    bibtex += f"  year = {{{year}}},\n"
                if journal:
                    bibtex += f"  journal = {{{journal}}},\n"
                if volume:
                    bibtex += f"  volume = {{{volume}}},\n"
                if issue:
                    bibtex += f"  number = {{{issue}}},\n"
                if pages:
                    bibtex += f"  pages = {{{pages}}},\n"
                if doi:
                    bibtex += f"  doi = {{{doi}}},\n"
                
                bibtex += "}"
                return bibtex
            
            # Default to a simple citation format
            return f"{', '.join(authors[:3])} {'et al. ' if len(authors) > 3 else ''}({year}). {title}. {journal}. DOI: {doi}"
        
        except Exception as e:
            logger.error(f"Error formatting citation: {e}")
            return ""


class PubMedAPI:
    """
    Integration with PubMed API (E-utilities) for citation information.
    
    This class handles interaction with NCBI's E-utilities API to search
    PubMed, retrieve publication metadata, and get citation information.
    """
    
    ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    ESUMMARY_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    ELINK_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi"
    CACHE_EXPIRY = 30  # days
    
    def __init__(self, api_key: Optional[str] = None, email: Optional[str] = None):
        """
        Initialize the PubMed API client.
        
        Args:
            api_key: NCBI API key for higher rate limits
            email: Contact email for NCBI API usage
        """
        self.api_key = api_key
        self.email = email
        
        # NCBI limits to 3 requests per second without API key, 10 with API key
        self.rate_limiter = RateLimiter(calls_per_minute=180 if api_key else 60)
        
        # Log warning if no API key provided
        if not api_key:
            logger.warning("No NCBI API key provided. Rate limits will be restricted (3 requests/second).")
    
    def _get_cache_key(self, query: str) -> str:
        """Generate a standardized cache key for the query."""
        return f"pubmed_{base64.urlsafe_b64encode(query.encode()).decode()}"
    
    def _build_params(self, base_params: Dict[str, Any]) -> Dict[str, Any]:
        """Build API parameters with API key and tool information."""
        params = base_params.copy()
        
        # Add API key if available
        if self.api_key:
            params['api_key'] = self.api_key
            
        # Add tool and email information (good API citizenship)
        params['tool'] = 'awesome-virome'
        if self.email:
            params['email'] = self.email
            
        return params
    
    def search_publications(self, query: str, max_results: int = 10, repo_url: Optional[str] = None) -> List[Dict[str, Any]]:
        """
        Search for publications in PubMed using the E-utilities API.
        
        Args:
            query: Search query string
            max_results: Maximum number of results to return
            repo_url: Optional repository URL to associate this cache with
        """
        cache_key = self._get_cache_key(f"search_{query}_{max_results}")
        
        # Check if we have a valid cached result
        cached_data = cache_manager.get(cache_key)
        if cached_data is not None:
            logger.info(f"Using cached PubMed search results for '{query}'")
            return cached_data
        
        self.rate_limiter.wait()
        
        # Step 1: Search for PMIDs matching the query
        search_params = self._build_params({
            'db': 'pubmed',
            'term': query,
            'retmode': 'json',
            'retmax': max_results,
            'sort': 'relevance'
        })
        
        try:
            response = requests.get(self.ESEARCH_URL, params=search_params)
            response.raise_for_status()
            
            search_results = response.json()
            pmids = search_results.get('esearchresult', {}).get('idlist', [])
            
            if not pmids:
                logger.info(f"No PubMed results found for query: {query}")
                cache_manager.set(cache_key, [], repo_url)
                return []
            
            # Step 2: Get detailed information for the found PMIDs
            self.rate_limiter.wait()
            pmid_string = ",".join(pmids)
            summary_params = self._build_params({
                'db': 'pubmed',
                'id': pmid_string,
                'retmode': 'json'
            })
            
            summary_response = requests.get(self.ESUMMARY_URL, params=summary_params)
            summary_response.raise_for_status()
            
            summary_results = summary_response.json()
            result_list = []
            
            # Process the results
            result_dict = summary_results.get('result', {})
            for pmid in pmids:
                if pmid in result_dict:
                    pub_data = result_dict[pmid]
                    
                    # Extract the DOI if available
                    doi = None
                    for id_type in pub_data.get('articleids', []):
                        if id_type.get('idtype') == 'doi':
                            doi = id_type.get('value')
                    
                    # Format authors
                    authors = []
                    for author in pub_data.get('authors', []):
                        authors.append(author.get('name', ''))
                    
                    # Extract publication data
                    publication = {
                        'pmid': pmid,
                        'doi': doi,
                        'title': pub_data.get('title', ''),
                        'authors': authors,
                        'journal': pub_data.get('fulljournalname', ''),
                        'publication_date': pub_data.get('pubdate', ''),
                        'year': pub_data.get('pubdate', '')[:4] if pub_data.get('pubdate', '') else '',
                        'volume': pub_data.get('volume', ''),
                        'issue': pub_data.get('issue', ''),
                        'pages': pub_data.get('pages', ''),
                        'pub_types': pub_data.get('pubtypes', []),
                        'is_software': 'Software' in pub_data.get('pubtypes', []) or 'software' in pub_data.get('title', '').lower(),
                        'is_tool': 'tool' in pub_data.get('title', '').lower() or 'pipeline' in pub_data.get('title', '').lower()
                    }
                    
                    result_list.append(publication)
            
            # Cache the results and associate with repo_url if provided
            cache_manager.set(cache_key, result_list, repo_url)
            
            return result_list
            
        except requests.RequestException as e:
            logger.error(f"Error searching PubMed: {e}")
            return []
    
    def get_publication_by_pmid(self, pmid: str, repo_url: Optional[str] = None) -> Dict[str, Any]:
        """
        Get detailed publication information by PMID.
        
        Args:
            pmid: PubMed ID
            repo_url: Optional repository URL to associate this cache with
        """
        cache_key = self._get_cache_key(f"pmid_{pmid}")
        
        # Check if we have a valid cached result
        cached_data = cache_manager.get(cache_key)
        if cached_data is not None:
            logger.info(f"Using cached PubMed publication data for PMID: {pmid}")
            return cached_data
        
        self.rate_limiter.wait()
        summary_params = self._build_params({
            'db': 'pubmed',
            'id': pmid,
            'retmode': 'json'
        })
        
        try:
            response = requests.get(self.ESUMMARY_URL, params=summary_params)
            response.raise_for_status()
            
            result = response.json().get('result', {})
            if pmid not in result:
                logger.warning(f"PMID {pmid} not found in PubMed")
                cache_manager.set(cache_key, {}, repo_url)
                return {}
            
            pub_data = result[pmid]
            
            # Extract the DOI if available
            doi = None
            for id_type in pub_data.get('articleids', []):
                if id_type.get('idtype') == 'doi':
                    doi = id_type.get('value')
            
            # Format authors
            authors = []
            for author in pub_data.get('authors', []):
                authors.append(author.get('name', ''))
            
            # Extract publication data
            publication = {
                'pmid': pmid,
                'doi': doi,
                'title': pub_data.get('title', ''),
                'authors': authors,
                'journal': pub_data.get('fulljournalname', ''),
                'publication_date': pub_data.get('pubdate', ''),
                'year': pub_data.get('pubdate', '')[:4] if pub_data.get('pubdate', '') else '',
                'volume': pub_data.get('volume', ''),
                'issue': pub_data.get('issue', ''),
                'pages': pub_data.get('pages', ''),
                'pub_types': pub_data.get('pubtypes', []),
                'abstract': self.get_abstract(pmid, repo_url)
            }
            
            # Cache the results and associate with repo_url if provided
            cache_manager.set(cache_key, publication, repo_url)
            
            return publication
            
        except requests.RequestException as e:
            logger.error(f"Error getting publication from PubMed: {e}")
            return {}
    
    def get_abstract(self, pmid: str, repo_url: Optional[str] = None) -> str:
        """
        Get the abstract for a publication by PMID.
        
        Args:
            pmid: PubMed ID
            repo_url: Optional repository URL to associate this cache with
        """
        cache_key = self._get_cache_key(f"abstract_{pmid}")
        
        # Check if we have a valid cached result
        cached_data = cache_manager.get(cache_key)
        if cached_data is not None:
            logger.info(f"Using cached PubMed abstract for PMID: {pmid}")
            return cached_data
        
        self.rate_limiter.wait()
        fetch_params = self._build_params({
            'db': 'pubmed',
            'id': pmid,
            'rettype': 'abstract',
            'retmode': 'text'
        })
        
        try:
            response = requests.get(self.EFETCH_URL, params=fetch_params)
            response.raise_for_status()
            
            abstract = response.text
            
            # Cache the results and associate with repo_url if provided
            cache_manager.set(cache_key, abstract, repo_url)
            
            return abstract
            
        except requests.RequestException as e:
            logger.error(f"Error getting abstract from PubMed: {e}")
            return ""
    
    def get_citation_count(self, pmid: str, repo_url: Optional[str] = None) -> int:
        """
        Get citation count for a publication using PubMed Central.
        
        Args:
            pmid: PubMed ID
            repo_url: Optional repository URL to associate this cache with
        """
        cache_key = self._get_cache_key(f"citations_{pmid}")
        
        # Check if we have a valid cached result
        cached_data = cache_manager.get(cache_key)
        if cached_data is not None:
            logger.info(f"Using cached PubMed citation count for PMID: {pmid}")
            return cached_data
        
        self.rate_limiter.wait()
        link_params = self._build_params({
            'dbfrom': 'pubmed',
            'db': 'pubmed',
            'id': pmid,
            'cmd': 'citation',
            'retmode': 'json'
        })
        
        try:
            response = requests.get(self.ELINK_URL, params=link_params)
            response.raise_for_status()
            
            data = response.json()
            
            # Extract the citation count
            citations = 0
            link_sets = data.get('linksets', [])
            for link_set in link_sets:
                if 'linksetdbs' in link_set:
                    for link_db in link_set['linksetdbs']:
                        if link_db.get('linkname') == 'pubmed_pubmed_citedin':
                            citations = len(link_db.get('links', []))
            
            # Cache the result and associate with repo_url if provided
            cache_manager.set(cache_key, citations, repo_url)
            
            return citations
            
        except requests.RequestException as e:
            logger.error(f"Error getting citation count from PubMed: {e}")
            return 0
    
    def format_citation(self, publication: Dict[str, Any], style: str = "apa") -> str:
        """
        Format a citation in the specified style.
        
        Args:
            publication: Publication data
            style: Citation style to use ('apa', 'bibtex', 'mla', 'chicago')
        """
        if not publication:
            return ""
        
        try:
            # Extract publication data
            title = publication.get('title', '').rstrip('.')
            authors = publication.get('authors', [])
            year = publication.get('year', '')
            journal = publication.get('journal', '')
            volume = publication.get('volume', '')
            issue = publication.get('issue', '')
            pages = publication.get('pages', '')
            doi = publication.get('doi', '')
            pmid = publication.get('pmid', '')
            
            if style.lower() == "apa":
                # Format authors for APA style
                if len(authors) == 0:
                    author_text = ""
                elif len(authors) == 1:
                    # Get last name and initials
                    parts = authors[0].split()
                    if len(parts) > 1:
                        last_name = parts[-1]
                        initials = ''.join([p[0] + '.' for p in parts[:-1]])
                        author_text = f"{last_name}, {initials}"
                    else:
                        author_text = authors[0]
                elif len(authors) < 8:
                    # Format each author properly
                    formatted_authors = []
                    for author in authors:
                        parts = author.split()
                        if len(parts) > 1:
                            last_name = parts[-1]
                            initials = ''.join([p[0] + '.' for p in parts[:-1]])
                            formatted_authors.append(f"{last_name}, {initials}")
                        else:
                            formatted_authors.append(author)
                    author_text = ", ".join(formatted_authors[:-1]) + f", & {formatted_authors[-1]}"
                else:
                    # APA 7th edition: first 6 authors, ..., last author
                    formatted_authors = []
                    for author in authors[:6]:
                        parts = author.split()
                        if len(parts) > 1:
                            last_name = parts[-1]
                            initials = ''.join([p[0] + '.' for p in parts[:-1]])
                            formatted_authors.append(f"{last_name}, {initials}")
                        else:
                            formatted_authors.append(author)
                    last_author = authors[-1].split()
                    if len(last_author) > 1:
                        last_name = last_author[-1]
                        initials = ''.join([p[0] + '.' for p in last_author[:-1]])
                        formatted_last = f"{last_name}, {initials}"
                    else:
                        formatted_last = authors[-1]
                    author_text = ", ".join(formatted_authors) + f", ... {formatted_last}"
                
                # Format APA style citation
                citation = f"{author_text} ({year}). {title}."
                if journal:
                    citation += f" {journal}"
                    if volume:
                        citation += f", {volume}"
                        if issue:
                            citation += f"({issue})"
                    if pages:
                        citation += f", {pages}"
                citation += f". https://doi.org/{doi}" if doi else ""
                
                return citation
            
            elif style.lower() == "bibtex":
                # Create a BibTeX key from first author's last name and year
                if authors and year:
                    first_author = authors[0].split()[-1].lower() if authors[0].split() else "unknown"
                    bibtex_key = f"{first_author}{year}"
                else:
                    bibtex_key = f"pubmed{pmid}"
                
                # Format BibTeX entry
                bibtex = f"@article{{{bibtex_key},\n"
                bibtex += f"  title = {{{title}}},\n"
                
                # Authors
                if authors:
                    # Format each author for BibTeX
                    formatted_authors = []
                    for author in authors:
                        parts = author.split()
                        if len(parts) > 1:
                            last_name = parts[-1]
                            given_names = ' '.join(parts[:-1])
                            formatted_authors.append(f"{last_name}, {given_names}")
                        else:
                            formatted_authors.append(author)
                    bibtex += f"  author = {{{' and '.join(formatted_authors)}}},\n"
                
                # Other fields
                if year:
                    bibtex += f"  year = {{{year}}},\n"
                if journal:
                    bibtex += f"  journal = {{{journal}}},\n"
                if volume:
                    bibtex += f"  volume = {{{volume}}},\n"
                if issue:
                    bibtex += f"  number = {{{issue}}},\n"
                if pages:
                    bibtex += f"  pages = {{{pages}}},\n"
                if doi:
                    bibtex += f"  doi = {{{doi}}},\n"
                if pmid:
                    bibtex += f"  pmid = {{{pmid}}},\n"
                
                bibtex += "}"
                return bibtex
            
            # MLA style
            elif style.lower() == "mla":
                # Format authors for MLA style
                if len(authors) == 0:
                    author_text = ""
                elif len(authors) == 1:
                    parts = authors[0].split()
                    if len(parts) > 1:
                        last_name = parts[-1]
                        first_names = ' '.join(parts[:-1])
                        author_text = f"{last_name}, {first_names}"
                    else:
                        author_text = authors[0]
                elif len(authors) == 2:
                    parts1 = authors[0].split()
                    if len(parts1) > 1:
                        last_name1 = parts1[-1]
                        first_names1 = ' '.join(parts1[:-1])
                        author1 = f"{last_name1}, {first_names1}"
                    else:
                        author1 = authors[0]
                    
                    parts2 = authors[1].split()
                    if len(parts2) > 1:
                        first_names2 = ' '.join(parts2[:-1])
                        last_name2 = parts2[-1]
                        author2 = f"{first_names2} {last_name2}"
                    else:
                        author2 = authors[1]
                    
                    author_text = f"{author1}, and {author2}"
                else:
                    parts = authors[0].split()
                    if len(parts) > 1:
                        last_name = parts[-1]
                        first_names = ' '.join(parts[:-1])
                        author_text = f"{last_name}, {first_names}, et al"
                    else:
                        author_text = f"{authors[0]}, et al"
                
                # Format MLA style citation
                citation = f"{author_text}. \"{title}.\" {journal}"
                if volume or issue:
                    citation += ", vol." if volume else ""
                    citation += f" {volume}" if volume else ""
                    citation += ", no." if issue else ""
                    citation += f" {issue}" if issue else ""
                citation += f", {year}" if year else ""
                citation += f", pp. {pages}" if pages else ""
                citation += f". DOI: {doi}" if doi else ""
                
                return citation
            
            # Default to a simple citation format
            return f"{', '.join(authors[:3])} {'et al. ' if len(authors) > 3 else ''}({year}). {title}. {journal}. DOI: {doi}"
        
        except Exception as e:
            logger.error(f"Error formatting citation: {e}")
            return ""
    
    def find_best_publication_for_tool(self, tool_name: str, repo_url: Optional[str] = None) -> Dict[str, Any]:
        """
        Find the best publication match for a bioinformatics tool.
        
        Args:
            tool_name: Name of the tool
            repo_url: Optional repository URL
        """
        # Try different search queries in order of specificity
        search_queries = [
            f'"{tool_name}"[Title] AND (software[Publication Type] OR tool[Title] OR pipeline[Title])',
            f'"{tool_name}"[Title] AND (bioinformatics OR computational OR algorithm)',
            f'"{tool_name}"[Title] AND (virus OR viral OR virome OR phage)',
            f'"{tool_name}"[Title/Abstract] AND (software[Publication Type] OR tool OR pipeline)',
            f'"{tool_name}"[Title/Abstract] AND (bioinformatics OR computational)',
            f'"{tool_name.replace("-", " ")}"[Title] AND (software OR tool OR pipeline)'
        ]
        
        for query in search_queries:
            results = self.search_publications(query, max_results=5, repo_url=repo_url)
            
            if results:
                # Prioritize software publications
                software_pubs = [pub for pub in results if pub.get('is_software', False)]
                if software_pubs:
                    best_match = software_pubs[0]
                    logger.info(f"Found software publication for {tool_name}: {best_match.get('title')}")
                    
                    # Get citation count
                    pmid = best_match.get('pmid')
                    if pmid:
                        citation_count = self.get_citation_count(pmid, repo_url)
                        best_match['citation_count'] = citation_count
                    
                    return best_match
                
                # If no software publications, return the first result
                best_match = results[0]
                logger.info(f"Found publication for {tool_name}: {best_match.get('title')}")
                
                # Get citation count
                pmid = best_match.get('pmid')
                if pmid:
                    citation_count = self.get_citation_count(pmid, repo_url)
                    best_match['citation_count'] = citation_count
                
                return best_match
        
        logger.info(f"No relevant publications found for {tool_name}")
        return {}


class GitHubAPI:
    """Integration with GitHub API to extract citation information."""
    
    BASE_URL = "https://api.github.com"
    CACHE_EXPIRY = 7  # days
    
    def __init__(self, token: Optional[str] = None):
        self.token = token
        # Use a more conservative rate limit if no token is provided
        self.rate_limiter = RateLimiter(calls_per_minute=5 if not token else 30)
        
        # Log warning if no token provided
        if not token:
            logger.warning("No GitHub token provided. Rate limits will be severely restricted (60 requests/hour).")
    
    def _get_cache_key(self, query: str) -> str:
        """Generate a standardized cache key for the query."""
        return f"github_{base64.urlsafe_b64encode(query.encode()).decode()}"
    
    def get_repo_contents(self, owner: str, repo: str, path: str = "", repo_url: Optional[str] = None) -> List[Dict[str, Any]]:
        """
        Get contents of a repository at the given path.
        
        Args:
            owner: Repository owner
            repo: Repository name
            path: Path within the repository
            repo_url: Optional full repository URL to associate this cache with
        """
        cache_key = self._get_cache_key(f"contents_{owner}_{repo}_{path}")
        repo_url = repo_url or f"https://github.com/{owner}/{repo}"
        
        # Check if we have a valid cached result
        cached_data = cache_manager.get(cache_key)
        if cached_data is not None:
            logger.info(f"Using cached GitHub contents for '{owner}/{repo}/{path}'")
            return cached_data
        
        self.rate_limiter.wait()
        url = f"{self.BASE_URL}/repos/{owner}/{repo}/contents/{path}"
        
        headers = {
            'Accept': 'application/vnd.github.v3+json'
        }
        if self.token:
            headers['Authorization'] = f"token {self.token}"
        
        try:
            response = requests.get(url, headers=headers)
            
            # Check for rate limiting and use the enhanced rate limiter
            if response.status_code == 403 and 'rate limit exceeded' in response.text.lower():
                # Update the rate limiter with the response headers
                self.rate_limiter.update_rate_limit(response.headers)
                
                # Handle rate limit with intelligent retry logic
                if self.rate_limiter.handle_rate_limit_response(response):
                    # Retry the request after waiting
                    return self.get_repo_contents(owner, repo, path, repo_url)
                
                logger.error("GitHub rate limit exceeded and wait time too long. Skipping request.")
                return []
            
            response.raise_for_status()
            contents = response.json()
            
            # Cache the results and associate with repo_url
            result = contents if isinstance(contents, list) else [contents]
            cache_manager.set(cache_key, result, repo_url)
            
            return result
        except requests.RequestException as e:
            if isinstance(e, requests.exceptions.HTTPError) and e.response.status_code == 404:
                # Log 404s as info rather than error - this is expected for many repos
                logger.info(f"File not found on GitHub: {owner}/{repo}/{path}")
            else:
                logger.error(f"Error getting repo contents from GitHub: {e}")
            
            # Cache an empty list for 404s to avoid repeated lookups
            if isinstance(e, requests.exceptions.HTTPError) and e.response.status_code == 404:
                cache_manager.set(cache_key, [], repo_url)
                
            return []
    
    def get_file_content(self, owner: str, repo: str, path: str, repo_url: Optional[str] = None) -> Optional[str]:
        """
        Get the content of a specific file.
        
        Args:
            owner: Repository owner
            repo: Repository name
            path: Path to the file within the repository
            repo_url: Optional full repository URL to associate this cache with
        """
        cache_key = self._get_cache_key(f"file_{owner}_{repo}_{path}")
        repo_url = repo_url or f"https://github.com/{owner}/{repo}"
        
        # Check if we have a valid cached result
        cached_data = cache_manager.get(cache_key)
        if cached_data is not None:
            logger.info(f"Using cached GitHub file for '{owner}/{repo}/{path}'")
            return cached_data
        
        self.rate_limiter.wait()
        url = f"{self.BASE_URL}/repos/{owner}/{repo}/contents/{path}"
        
        headers = {
            'Accept': 'application/vnd.github.v3+json'
        }
        if self.token:
            headers['Authorization'] = f"token {self.token}"
        
        try:
            response = requests.get(url, headers=headers)
            
            # Check for rate limiting and use the enhanced rate limiter
            if response.status_code == 403 and 'rate limit exceeded' in response.text.lower():
                # Update the rate limiter with the response headers
                self.rate_limiter.update_rate_limit(response.headers)
                
                # Handle rate limit with intelligent retry logic
                if self.rate_limiter.handle_rate_limit_response(response):
                    # Retry the request after waiting
                    return self.get_file_content(owner, repo, path, repo_url)
                
                logger.error("GitHub rate limit exceeded and wait time too long. Skipping request.")
                return None
            
            response.raise_for_status()
            data = response.json()
            
            if isinstance(data, dict) and 'content' in data:
                content = base64.b64decode(data['content']).decode('utf-8')
                
                # Cache the content and associate with repo_url
                cache_manager.set(cache_key, content, repo_url)
                
                return content
            return None
        except requests.RequestException as e:
            if isinstance(e, requests.exceptions.HTTPError) and e.response.status_code == 404:
                # Log 404s as info rather than error - this is expected for many repos
                logger.info(f"Error getting file content from GitHub: 404 Not Found for {owner}/{repo}/{path}")
            else:
                logger.error(f"Error getting file content from GitHub: {e}")
            return None
    
    def find_citation_file(self, repo_url: str) -> Optional[Dict[str, Any]]:
        """Find CITATION.cff or other citation files in a repository."""
        # Skip if no repo_url provided
        if not repo_url or not isinstance(repo_url, str):
            logger.warning("Invalid or missing repository URL")
            return None
            
        # Skip non-GitHub repositories
        if 'github.com' not in repo_url:
            logger.info(f"Not a GitHub repository: {repo_url}")
            return None
        
        # Extract owner and repo from URL
        match = re.search(r'github\.com/([^/]+)/([^/]+)', repo_url)
        if not match:
            logger.error(f"Could not extract owner/repo from URL: {repo_url}")
            return None
        
        try:
            owner, repo = match.groups()
            repo = repo.rstrip('.git')
            
            cache_key = self._get_cache_key(f"citation_search_{owner}_{repo}")
            
            # Check if we have a valid cached result
            cached_data = cache_manager.get(cache_key)
            if cached_data is not None:
                logger.info(f"Using cached citation search results for '{owner}/{repo}'")
                return cached_data
            
            # Common citation file paths to check
            citation_paths = [
                "CITATION.cff",
                "CITATION",
                "citation.cff",
                "citation.txt",
                "CITATION.txt",
                "CITATION.md",
                "citation.md",
                "CITATION.bib",
                "citation.bib"
            ]
            
            # Search for citation files in a more controlled manner
            result = None
            
            # First check main citation files (most common)
            for path in citation_paths[:3]:  # Limit to the most common ones first
                try:
                    content = self.get_file_content(owner, repo, path, repo_url)
                    if content:
                        result = {
                            'path': path,
                            'content': content,
                            'format': path.split('.')[-1] if '.' in path else 'txt'
                        }
                        break
                except Exception as e:
                    logger.warning(f"Error checking for citation file {path}: {e}")
            
            # If no result from common files, try .github directory
            if not result:
                try:
                    # Only check the most common location in .github
                    github_path = ".github/CITATION.cff"
                    content = self.get_file_content(owner, repo, github_path, repo_url)
                    if content:
                        result = {
                            'path': github_path,
                            'content': content,
                            'format': 'cff'
                        }
                except Exception as e:
                    logger.warning(f"Error checking for citation in .github directory: {e}")
            
            # As a last resort, check README.md
            if not result:
                try:
                    readme_content = self.get_file_content(owner, repo, "README.md", repo_url)
                    if readme_content:
                        # Look for citation section in README
                        citation_section = None
                        citation_patterns = [
                            r'(?i)#+\s*citation\s*\n+(.*?)(?:\n#+\s|\Z)',
                            r'(?i)#+\s*how\s+to\s+cite\s*\n+(.*?)(?:\n#+\s|\Z)',
                            r'(?i)#+\s*citing\s+this\s+work\s*\n+(.*?)(?:\n#+\s|\Z)'
                        ]
                        
                        for pattern in citation_patterns:
                            match = re.search(pattern, readme_content, re.DOTALL)
                            if match:
                                citation_section = match.group(1).strip()
                                break
                        
                        if citation_section:
                            result = {
                                'path': 'README.md',
                                'content': citation_section,
                                'format': 'md'
                            }
                except Exception as e:
                    logger.warning(f"Error checking for citation in README: {e}")
            
            # Cache the result and associate with repo_url
            cache_manager.set(cache_key, result if result else {}, repo_url)
            
            return result
        
        except Exception as e:
            logger.error(f"Error in find_citation_file for {repo_url}: {e}")
            return None
    
    def parse_citation_cff(self, content: str) -> Dict[str, Any]:
        """Parse CITATION.cff file for citation information."""
        try:
            import yaml
            citation_data = yaml.safe_load(content)
            
            # Extract relevant citation information
            citation = {
                'title': citation_data.get('title', ''),
                'authors': [
                    f"{author.get('family-names', '')} {author.get('given-names', '')}"
                    for author in citation_data.get('authors', [])
                ],
                'year': citation_data.get('date-released', '').split('-')[0] if citation_data.get('date-released') else '',
                'doi': citation_data.get('doi', ''),
                'version': citation_data.get('version', ''),
                'url': citation_data.get('url', ''),
                'repository': citation_data.get('repository-code', '')
            }
            
            # Format according to specified styles
            formatted = {}
            if citation_data.get('preferred-citation'):
                preferred = citation_data.get('preferred-citation', {})
                if 'doi' in preferred:
                    # If the preferred citation has a DOI, we can get formatted citations later
                    citation['preferred_citation_doi'] = preferred.get('doi')
                else:
                    # Otherwise, build from the preferred citation data
                    authors = []
                    for author in preferred.get('authors', []):
                        given = author.get('given-names', '')
                        family = author.get('family-names', '')
                        if given and family:
                            authors.append(f"{family}, {given[0]}.")
                        elif family:
                            authors.append(family)
                    
                    journal = preferred.get('journal', '')
                    volume = preferred.get('volume', '')
                    issue = preferred.get('issue', '')
                    pages = preferred.get('pages', '')
                    year = preferred.get('year', '') or citation['year']
                    
                    # APA style
                    if authors:
                        if len(authors) == 1:
                            apa_authors = authors[0]
                        elif len(authors) < 8:
                            apa_authors = ", ".join(authors[:-1]) + f", & {authors[-1]}"
                        else:
                            apa_authors = ", ".join(authors[:6]) + f", ... {authors[-1]}"
                        
                        apa = f"{apa_authors} ({year}). {preferred.get('title', '')}."
                        if journal:
                            apa += f" {journal}"
                            if volume:
                                apa += f", {volume}"
                                if issue:
                                    apa += f"({issue})"
                            if pages:
                                apa += f", {pages}"
                        formatted['apa'] = apa
            
            citation['formatted'] = formatted
            return citation
        
        except Exception as e:
            logger.error(f"Error parsing CITATION.cff: {e}")
            return {}
    
    def extract_doi_from_text(self, text: str) -> Optional[str]:
        """Extract DOI from text content."""
        # DOI pattern
        doi_pattern = r'(?:doi:|DOI:|https?://doi.org/)(10\.\d+/[^\s\'"]+)'
        match = re.search(doi_pattern, text)
        if match:
            return match.group(1)
        return None