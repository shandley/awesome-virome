#!/usr/bin/env python3
"""
Enhanced cache management system with automatic expiration and LRU eviction.
Extends the existing CacheManager with improved features while maintaining compatibility.
"""

import os
import json
import time
import logging
import hashlib
import heapq
import threading
from datetime import datetime, timedelta
from typing import Dict, List, Any, Optional, Set, Tuple
from pathlib import Path

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class EnhancedCacheManager:
    """
    Enhanced cache manager with automatic expiration, LRU eviction, and size limiting.
    Extends the original CacheManager while maintaining backward compatibility.
    """
    
    def __init__(self, 
                cache_dir: str, 
                expiry_days: int = 30,
                max_size_mb: int = 500,
                lru_check_interval: int = 3600,  # Check LRU every hour
                size_limit_percentage: float = 0.9):  # Start cleanup at 90% of max size
        """
        Initialize the enhanced cache manager.
        
        Args:
            cache_dir: Base directory for cache storage
            expiry_days: Default expiration time in days
            max_size_mb: Maximum cache size in MB (default 500MB)
            lru_check_interval: How often to check for LRU eviction (seconds)
            size_limit_percentage: When to start cleanup (as percentage of max size)
        """
        self.base_cache_dir = Path(cache_dir)
        self.base_cache_dir.mkdir(exist_ok=True, parents=True)
        
        # Create a directory for dependency tracking
        self.deps_dir = self.base_cache_dir / "_dependencies"
        self.deps_dir.mkdir(exist_ok=True)
        
        # Create a directory for performance metrics
        self.metrics_dir = self.base_cache_dir / "_metrics"
        self.metrics_dir.mkdir(exist_ok=True)

        # Create directory for LRU tracking
        self.lru_dir = self.base_cache_dir / "_lru"
        self.lru_dir.mkdir(exist_ok=True)
        
        # Configuration
        self.default_expiry = timedelta(days=expiry_days)
        self.max_size_bytes = max_size_mb * 1024 * 1024
        self.size_limit_percentage = size_limit_percentage
        self.lru_check_interval = lru_check_interval
        
        # LRU and size tracking
        self._lru_data = {}  # Maps cache key hash to last access time
        self._size_data = 0  # Current cache size in bytes
        self._lru_lock = threading.RLock()  # Thread safety for LRU data
        
        # Map of repo URLs to their dependency caches
        self._dependency_map: Dict[str, Set[str]] = {}
        
        # Performance metrics
        self._metrics = {
            "hits": 0,
            "misses": 0,
            "invalidations": 0,
            "sets": 0,
            "expirations": 0,
            "lru_evictions": 0,
            "size_evictions": 0,
            "start_time": datetime.now().isoformat()
        }
        
        # Load existing data
        self._load_dependency_maps()
        self._load_metrics()
        self._load_lru_data()
        self._calculate_current_size()
        
        # Start background maintenance thread
        self._stop_maintenance = False
        self._maintenance_thread = threading.Thread(
            target=self._maintenance_loop,
            daemon=True
        )
        self._maintenance_thread.start()
    
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
                    for key in ['hits', 'misses', 'invalidations', 'sets', 
                              'expirations', 'lru_evictions', 'size_evictions']:
                        if key in stored_metrics:
                            self._metrics[key] = stored_metrics[key]
            except (json.JSONDecodeError, IOError) as e:
                logger.warning(f"Failed to load cache metrics: {e}")
    
    def _load_lru_data(self):
        """Load LRU tracking data."""
        lru_file = self.lru_dir / "lru_data.json"
        if lru_file.exists():
            try:
                with open(lru_file, 'r') as f:
                    self._lru_data = json.load(f)
            except (json.JSONDecodeError, IOError) as e:
                logger.warning(f"Failed to load LRU data: {e}")
                self._lru_data = {}
        
        # Initialize LRU data for existing cache files if needed
        for cache_file in self.base_cache_dir.glob("*.json"):
            if cache_file.is_file() and not cache_file.name.startswith("_"):
                key_hash = cache_file.stem
                if key_hash not in self._lru_data:
                    # Set default timestamp to file modification time
                    self._lru_data[key_hash] = cache_file.stat().st_mtime
    
    def _save_lru_data(self):
        """Save LRU tracking data to disk."""
        lru_file = self.lru_dir / "lru_data.json"
        try:
            with self._lru_lock:
                with open(lru_file, 'w') as f:
                    json.dump(self._lru_data, f)
        except IOError as e:
            logger.warning(f"Failed to save LRU data: {e}")
    
    def _calculate_current_size(self):
        """Calculate the current total size of the cache."""
        total_size = 0
        for cache_file in self.base_cache_dir.glob("*.json"):
            if cache_file.is_file() and not cache_file.name.startswith("_"):
                total_size += cache_file.stat().st_size
        
        self._size_data = total_size
        logger.debug(f"Current cache size: {self._size_data / (1024 * 1024):.2f} MB")
    
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
            current_metrics['cache_size_mb'] = self._size_data / (1024 * 1024)
            
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
    
    def _update_lru(self, cache_key: str):
        """Update the LRU data for a cache key."""
        hashed_key = self._hash_key(cache_key)
        with self._lru_lock:
            self._lru_data[hashed_key] = time.time()
        
        # Periodically save LRU data (don't save on every access for performance)
        if hash(cache_key) % 20 == 0:  # Save roughly every 20 accesses
            self._save_lru_data()
    
    def _check_cache_size(self):
        """Check if cache size exceeds limits and perform eviction if needed."""
        # Skip if size is below threshold
        threshold = self.max_size_bytes * self.size_limit_percentage
        if self._size_data < threshold:
            return
        
        # Calculate how much to remove (target 80% of max)
        target_size = self.max_size_bytes * 0.8
        bytes_to_remove = self._size_data - target_size
        if bytes_to_remove <= 0:
            return
        
        logger.info(f"Cache size ({self._size_data / (1024*1024):.2f} MB) exceeds threshold. "
                   f"Need to remove {bytes_to_remove / (1024*1024):.2f} MB")
        
        bytes_removed = self._evict_by_size(bytes_to_remove)
        logger.info(f"Removed {bytes_removed / (1024*1024):.2f} MB from cache")
    
    def _evict_by_lru(self, count: int) -> int:
        """
        Evict least recently used cache entries.
        
        Args:
            count: Number of entries to evict
            
        Returns:
            int: Number of entries actually evicted
        """
        with self._lru_lock:
            # Create a list of (timestamp, key_hash) tuples
            lru_items = [(timestamp, key_hash) for key_hash, timestamp in self._lru_data.items()]
            
            # Sort by timestamp (oldest first)
            lru_items.sort()
            
            # Take the oldest 'count' items
            items_to_evict = lru_items[:count]
            
        evicted = 0
        total_bytes_removed = 0
        
        for _, key_hash in items_to_evict:
            cache_path = self.base_cache_dir / f"{key_hash}.json"
            if cache_path.exists():
                try:
                    # Track size being removed
                    size = cache_path.stat().st_size
                    
                    # Remove the file
                    cache_path.unlink()
                    
                    # Update tracking
                    with self._lru_lock:
                        if key_hash in self._lru_data:
                            del self._lru_data[key_hash]
                    
                    evicted += 1
                    total_bytes_removed += size
                    self._metrics['lru_evictions'] += 1
                except IOError as e:
                    logger.warning(f"Failed to evict cache file {cache_path}: {e}")
        
        # Update size tracking
        self._size_data -= total_bytes_removed
        
        # Save updated LRU data and metrics
        self._save_lru_data()
        self._save_metrics()
        
        if evicted > 0:
            logger.info(f"LRU eviction: removed {evicted} items ({total_bytes_removed / (1024*1024):.2f} MB)")
        
        return evicted
    
    def _evict_by_size(self, bytes_to_remove: int) -> int:
        """
        Evict cache entries until the specified number of bytes is freed.
        
        Args:
            bytes_to_remove: Number of bytes to remove
            
        Returns:
            int: Number of bytes actually removed
        """
        with self._lru_lock:
            # Create a list of (timestamp, size, key_hash) tuples
            eviction_candidates = []
            for key_hash, timestamp in self._lru_data.items():
                cache_path = self.base_cache_dir / f"{key_hash}.json"
                if cache_path.exists():
                    size = cache_path.stat().st_size
                    eviction_candidates.append((timestamp, size, key_hash))
            
            # Sort by timestamp (oldest first)
            eviction_candidates.sort()
        
        bytes_removed = 0
        evicted = 0
        
        for timestamp, size, key_hash in eviction_candidates:
            cache_path = self.base_cache_dir / f"{key_hash}.json"
            if cache_path.exists():
                try:
                    # Remove the file
                    cache_path.unlink()
                    
                    # Update tracking
                    with self._lru_lock:
                        if key_hash in self._lru_data:
                            del self._lru_data[key_hash]
                    
                    bytes_removed += size
                    evicted += 1
                    self._metrics['size_evictions'] += 1
                    
                    # Stop if we've removed enough
                    if bytes_removed >= bytes_to_remove:
                        break
                except IOError as e:
                    logger.warning(f"Failed to evict cache file {cache_path}: {e}")
        
        # Update size tracking
        self._size_data -= bytes_removed
        
        # Save updated LRU data and metrics
        self._save_lru_data()
        self._save_metrics()
        
        if evicted > 0:
            logger.info(f"Size-based eviction: removed {evicted} items ({bytes_removed / (1024*1024):.2f} MB)")
        
        return bytes_removed
    
    def _remove_expired_caches(self) -> int:
        """
        Remove expired cache entries.
        
        Returns:
            int: Number of cache entries removed
        """
        expiry_time = datetime.now() - self.default_expiry
        expiry_timestamp = expiry_time.timestamp()
        
        removed = 0
        bytes_removed = 0
        
        for cache_file in self.base_cache_dir.glob("*.json"):
            if cache_file.is_file() and not cache_file.name.startswith("_"):
                # Check if file modification time is older than expiry
                if cache_file.stat().st_mtime < expiry_timestamp:
                    try:
                        # Track size being removed
                        size = cache_file.stat().st_size
                        
                        # Remove the file
                        cache_file.unlink()
                        
                        # Update LRU tracking
                        key_hash = cache_file.stem
                        with self._lru_lock:
                            if key_hash in self._lru_data:
                                del self._lru_data[key_hash]
                        
                        removed += 1
                        bytes_removed += size
                        self._metrics['expirations'] += 1
                    except IOError as e:
                        logger.warning(f"Failed to remove expired cache {cache_file}: {e}")
        
        # Update size tracking
        self._size_data -= bytes_removed
        
        # Save updated LRU data and metrics
        if removed > 0:
            self._save_lru_data()
            self._save_metrics()
            logger.info(f"Expired {removed} cache entries ({bytes_removed / (1024*1024):.2f} MB)")
        
        return removed
    
    def _maintenance_loop(self):
        """Background maintenance thread for cache management."""
        logger.info("Starting cache maintenance thread")
        last_expiry_check = 0
        last_lru_check = 0
        last_size_check = 0
        
        try:
            while not self._stop_maintenance:
                now = time.time()
                
                # Check for expired entries (every day)
                if now - last_expiry_check > 86400:  # 24 hours
                    self._remove_expired_caches()
                    last_expiry_check = now
                
                # Check LRU status (every hour by default)
                if now - last_lru_check > self.lru_check_interval:
                    # Only evict if we're using more than 70% of max size
                    if self._size_data > (self.max_size_bytes * 0.7):
                        # Evict 5% of entries
                        total_entries = len(self._lru_data)
                        self._evict_by_lru(max(1, int(total_entries * 0.05)))
                    last_lru_check = now
                
                # Check size limits (every 10 minutes)
                if now - last_size_check > 600:  # 10 minutes
                    self._check_cache_size()
                    last_size_check = now
                
                # Sleep for a while
                time.sleep(60)  # Check every minute
        except Exception as e:
            logger.error(f"Error in cache maintenance thread: {e}")
        finally:
            logger.info("Cache maintenance thread stopped")
    
    def stop_maintenance(self):
        """Stop the maintenance thread."""
        self._stop_maintenance = True
        if self._maintenance_thread.is_alive():
            self._maintenance_thread.join(timeout=2.0)
        
        # Save final state
        self._save_lru_data()
        self._save_metrics()
    
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
            
            # Update LRU tracking (this is valid and accessed data)
            self._update_lru(cache_key)
                
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
                
                # Update LRU tracking
                self._update_lru(cache_key)
                
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
        hashed_key = self._hash_key(cache_key)
        
        try:
            # Track cache set
            self._metrics['sets'] += 1
            
            # Add timing information to the cache
            cache_data = {
                'cache_date': datetime.now().isoformat(),
                'data': data
            }
            
            # Track size before and after to update size monitoring
            old_size = cache_path.stat().st_size if cache_path.exists() else 0
            
            with open(cache_path, 'w') as f:
                json.dump(cache_data, f, indent=2)
            
            # Update size monitoring
            new_size = cache_path.stat().st_size
            self._size_data = self._size_data - old_size + new_size
            
            # Update LRU tracking
            self._update_lru(cache_key)
            
            # Register dependency if a repo URL was provided
            if repo_url:
                self.register_dependency(repo_url, cache_key)
            
            # Save metrics
            self._save_metrics()
            
            # Check cache size after new data and perform size-based eviction if needed
            self._check_cache_size()
                
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
        bytes_removed = 0
        
        for cache_key in list(self._dependency_map[repo_id]):  # Create a copy of the set to safely iterate
            cache_path = self._get_cache_path(cache_key)
            hashed_key = self._hash_key(cache_key)
            
            if cache_path.exists():
                try:
                    # Track size being removed
                    size = cache_path.stat().st_size
                    
                    # Remove the file
                    cache_path.unlink()
                    
                    # Update LRU tracking
                    with self._lru_lock:
                        if hashed_key in self._lru_data:
                            del self._lru_data[hashed_key]
                    
                    invalidated_count += 1
                    bytes_removed += size
                    logger.debug(f"Invalidated cache: {cache_key}")
                except IOError as e:
                    logger.warning(f"Failed to invalidate cache {cache_key}: {e}")
        
        # Update size tracking
        self._size_data -= bytes_removed
        
        # Track invalidations
        if invalidated_count > 0:
            self._metrics['invalidations'] += invalidated_count
            self._save_metrics()
            self._save_lru_data()
            
        # Clear the dependency list after invalidation
        self._dependency_map[repo_id] = set()
        self._save_dependency_map(repo_id)
        
        logger.info(f"Invalidated {invalidated_count} cache items "
                  f"({bytes_removed / (1024*1024):.2f} MB) for repository: {repo_url}")
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
        hashed_key = self._hash_key(cache_key)
        
        if not cache_path.exists():
            return True
            
        try:
            # Track size being removed
            size = cache_path.stat().st_size
            
            # Remove the file
            cache_path.unlink()
            
            # Update size tracking
            self._size_data -= size
            
            # Update LRU tracking
            with self._lru_lock:
                if hashed_key in self._lru_data:
                    del self._lru_data[hashed_key]
            
            # Track invalidation
            self._metrics['invalidations'] += 1
            self._save_metrics()
            self._save_lru_data()
            
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
        total_bytes_removed = 0
        
        for cache_file in self.base_cache_dir.glob("*.json"):
            # Skip special directories
            if cache_file.is_dir() or cache_file.name.startswith("_"):
                continue
                
            try:
                # Track size being removed
                size = cache_file.stat().st_size
                
                # Remove the file
                cache_file.unlink()
                
                count += 1
                total_bytes_removed += size
            except IOError as e:
                logger.warning(f"Failed to clear cache file {cache_file}: {e}")
        
        # Reset LRU tracking
        with self._lru_lock:
            self._lru_data = {}
        self._save_lru_data()
        
        # Reset size tracking
        self._size_data = 0
        
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
        
        logger.info(f"Cleared {count} cache items ({total_bytes_removed / (1024*1024):.2f} MB)")
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
        
        # Add file and size information
        cache_files = list(self.base_cache_dir.glob("*.json"))
        cache_files = [f for f in cache_files if f.is_file() and not f.name.startswith("_")]
        
        metrics['cache_files'] = len(cache_files)
        metrics['dependency_maps'] = len(self._dependency_map)
        metrics['cache_size_mb'] = self._size_data / (1024 * 1024)
        metrics['max_size_mb'] = self.max_size_bytes / (1024 * 1024)
        metrics['size_percentage'] = (
            self._size_data / self.max_size_bytes * 100 if self.max_size_bytes > 0 else 0
        )
        
        return metrics
    
    def get_lru_stats(self) -> Dict[str, Any]:
        """
        Get LRU statistics.
        
        Returns:
            Dict containing LRU statistics
        """
        with self._lru_lock:
            lru_count = len(self._lru_data)
            
            # Get oldest and newest timestamps
            oldest = min(self._lru_data.values()) if self._lru_data else 0
            newest = max(self._lru_data.values()) if self._lru_data else 0
            
            # Calculate age distribution
            now = time.time()
            age_buckets = {
                "1h": 0,   # Less than 1 hour
                "1d": 0,   # Less than 1 day
                "1w": 0,   # Less than 1 week
                "1m": 0,   # Less than 1 month
                "old": 0   # Older than 1 month
            }
            
            for timestamp in self._lru_data.values():
                age_seconds = now - timestamp
                
                if age_seconds < 3600:  # 1 hour
                    age_buckets["1h"] += 1
                elif age_seconds < 86400:  # 1 day
                    age_buckets["1d"] += 1
                elif age_seconds < 604800:  # 1 week
                    age_buckets["1w"] += 1
                elif age_seconds < 2592000:  # 1 month (30 days)
                    age_buckets["1m"] += 1
                else:
                    age_buckets["old"] += 1
            
            return {
                "count": lru_count,
                "oldest": datetime.fromtimestamp(oldest).isoformat() if oldest else None,
                "newest": datetime.fromtimestamp(newest).isoformat() if newest else None,
                "age_distribution": age_buckets
            }
    
    def get_cache_info(self) -> Dict[str, Any]:
        """
        Get comprehensive cache information.
        
        Returns:
            Dict containing all cache info
        """
        metrics = self.get_metrics()
        lru_stats = self.get_lru_stats()
        
        return {
            "metrics": metrics,
            "lru": lru_stats,
            "config": {
                "expiry_days": self.default_expiry.days,
                "max_size_mb": self.max_size_bytes / (1024 * 1024),
                "lru_check_interval": self.lru_check_interval,
                "size_limit_percentage": self.size_limit_percentage
            }
        }

# Create a backwards compatible function to get the enhanced cache manager
def get_enhanced_cache_manager(cache_dir=None, max_size_mb=500):
    """
    Get an instance of the enhanced cache manager.
    
    Args:
        cache_dir: Base directory for cache storage (default: metadata/cache)
        max_size_mb: Maximum cache size in MB (default: 500MB)
    """
    if cache_dir is None:
        cache_dir = os.path.join("metadata", "cache")
    
    return EnhancedCacheManager(cache_dir, max_size_mb=max_size_mb)