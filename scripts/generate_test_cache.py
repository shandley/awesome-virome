#!/usr/bin/env python3
"""
Script to generate test data for the cache monitoring system.
This script will create and populate the cache with test entries to simulate real cache activity.
"""

import os
import sys
import time
import random
import json
from pathlib import Path

# Add the parent directory to the path
script_dir = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.dirname(script_dir)
sys.path.append(script_dir)

# Import the cache manager
try:
    from apis.citations_api import cache_manager
    HAS_CACHE_MANAGER = True
except ImportError:
    HAS_CACHE_MANAGER = False
    print("Error: Cache manager not found. Make sure apis/citations_api.py exists.")
    sys.exit(1)

# Cache directory
CACHE_DIR = Path(repo_root) / "metadata" / "cache"
CACHE_DIR.mkdir(exist_ok=True, parents=True)

def generate_test_data():
    """Generate test cache entries for monitoring."""
    print("Generating test cache data...")
    
    # Test repositories to use for cache keys
    test_repos = [
        "https://github.com/deepmind/alphafold",
        "https://github.com/bbuchfink/diamond",
        "https://github.com/hyattpd/Prodigal",
        "https://github.com/ablab/spades",
        "https://github.com/biobakery/MetaPhlAn"
    ]
    
    # Generate cache keys
    cache_keys = []
    for repo in test_repos:
        cache_keys.append(f"citation_data_{repo}")
        cache_keys.append(f"api_response_{repo}")
        cache_keys.append(f"repo_info_{repo}")
    
    # Generate 30 cache entries
    hit_count = 0
    miss_count = 0
    
    for i in range(30):
        key = random.choice(cache_keys)
        
        # Simulate cache miss (first access)
        if not cache_manager.is_valid(key):
            print(f"Cache miss for {key}")
            miss_count += 1
            
            # Store mock data in cache
            data = {
                "test_id": i,
                "timestamp": time.time(),
                "repo": key.split('_', 2)[2] if len(key.split('_')) > 2 else "unknown",
                "data": f"Test data for {key}"
            }
            
            cache_manager.set(key, data)
        else:
            # Simulate cache hit
            print(f"Cache hit for {key}")
            hit_count += 1
            data = cache_manager.get(key)
        
        # Sleep to simulate real usage pattern
        time.sleep(0.1)
    
    # Simulate some invalidations
    invalidation_count = 0
    for i in range(5):
        repo = random.choice(test_repos)
        print(f"Invalidating cache for {repo}")
        count = cache_manager.invalidate_repo_caches(repo)
        invalidation_count += count
    
    # Print summary
    print("\nCache Test Data Generated")
    print("-------------------------")
    print(f"Total entries created: {hit_count + miss_count}")
    print(f"Cache hits: {hit_count}")
    print(f"Cache misses: {miss_count}")
    print(f"Cache invalidations: {invalidation_count}")
    print("-------------------------")
    print("Run python scripts/monitor_cache.py --graphs to see the results")

if __name__ == "__main__":
    generate_test_data()