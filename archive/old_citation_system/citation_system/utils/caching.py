#!/usr/bin/env python3
"""
Caching utilities for the citation system.
"""

import hashlib
import json
import logging
import time
from pathlib import Path
from typing import Any, Dict, Optional

from ..config import CACHE_DIR, CACHE_EXPIRY

logger = logging.getLogger(__name__)


def generate_cache_key(identifier: str, namespace: str = "") -> str:
    """
    Generate a cache key from an identifier and optional namespace.
    
    Args:
        identifier: Main identifier (e.g. DOI)
        namespace: Optional namespace to avoid key collisions
    
    Returns:
        A hash string to use as a cache key
    """
    if namespace:
        full_key = f"{namespace}:{identifier}"
    else:
        full_key = identifier
    
    # Create a hash of the identifier for use as a filename
    return hashlib.md5(full_key.encode('utf-8')).hexdigest()


def get_cache_path(cache_key: str) -> Path:
    """
    Get the path to a cache file for a given cache key.
    
    Args:
        cache_key: The cache key
    
    Returns:
        Path to the cache file
    """
    return CACHE_DIR / f"{cache_key}.json"


def save_to_cache(data: Dict[str, Any], cache_key: str) -> bool:
    """
    Save data to the cache.
    
    Args:
        data: Data to cache
        cache_key: Cache key
    
    Returns:
        True if saving was successful, False otherwise
    """
    # Add timestamp to cached data
    cache_data = {
        "timestamp": time.time(),
        "data": data
    }
    
    cache_path = get_cache_path(cache_key)
    
    try:
        with open(cache_path, 'w') as f:
            json.dump(cache_data, f, indent=2)
        return True
    except Exception as e:
        logger.warning(f"Failed to save data to cache: {e}")
        return False


def load_from_cache(
    cache_key: str, 
    max_age: Optional[int] = None
) -> Optional[Dict[str, Any]]:
    """
    Load data from the cache if it exists and is not expired.
    
    Args:
        cache_key: Cache key
        max_age: Maximum age of cache in seconds (defaults to CACHE_EXPIRY)
    
    Returns:
        Cached data or None if not found or expired
    """
    cache_path = get_cache_path(cache_key)
    
    if not cache_path.exists():
        return None
    
    try:
        with open(cache_path, 'r') as f:
            cache_data = json.load(f)
        
        # Check if cache is expired
        timestamp = cache_data.get("timestamp", 0)
        max_age = max_age if max_age is not None else CACHE_EXPIRY
        
        if time.time() - timestamp > max_age:
            logger.debug(f"Cache expired for key: {cache_key}")
            return None
        
        return cache_data.get("data")
    except Exception as e:
        logger.warning(f"Failed to load data from cache: {e}")
        return None


def clear_cache(namespace: Optional[str] = None) -> int:
    """
    Clear all cached data or data for a specific namespace.
    
    Args:
        namespace: Optional namespace to clear
    
    Returns:
        Number of cache files deleted
    """
    count = 0
    
    try:
        if namespace:
            # Delete only files for the given namespace
            for cache_file in CACHE_DIR.glob("*.json"):
                try:
                    with open(cache_file, 'r') as f:
                        cache_data = json.load(f)
                    
                    if cache_data.get("namespace") == namespace:
                        cache_file.unlink()
                        count += 1
                except Exception:
                    pass
        else:
            # Delete all cache files
            for cache_file in CACHE_DIR.glob("*.json"):
                cache_file.unlink()
                count += 1
        
        return count
    except Exception as e:
        logger.error(f"Error clearing cache: {e}")
        return count