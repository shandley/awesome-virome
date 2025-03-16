#!/usr/bin/env python3
"""
Utility script for managing the cache in the awesome-virome repository.

This script provides functionality to clear caches, including:
- Clearing all caches
- Clearing caches related to specific repositories
- Listing cache statistics
"""

import os
import sys
import argparse
import logging
from pathlib import Path

# Import the cache manager
try:
    from apis.citations_api import cache_manager
    HAS_CACHE_MANAGER = True
except ImportError:
    HAS_CACHE_MANAGER = False
    print("Error: Cache manager not found. Make sure apis/citations_api.py exists and is properly configured.")
    sys.exit(1)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def clear_all_caches():
    """Clear all caches."""
    count = cache_manager.clear_all()
    logger.info(f"Cleared {count} cache entries")
    return count

def clear_repo_caches(repo_url):
    """Clear caches related to a specific repository."""
    count = cache_manager.invalidate_repo_caches(repo_url)
    logger.info(f"Cleared {count} cache entries for repository: {repo_url}")
    return count

def get_cache_stats():
    """Get statistics about the current cache."""
    base_dir = cache_manager.base_cache_dir
    deps_dir = cache_manager.deps_dir
    
    # Count cache files
    cache_files = list(base_dir.glob("*.json"))
    cache_count = len([f for f in cache_files if not f.is_dir() and not f.name.startswith("_")])
    
    # Count dependency files
    deps_files = list(deps_dir.glob("*.json"))
    deps_count = len(deps_files)
    
    # Calculate total size
    total_size = sum(f.stat().st_size for f in cache_files if not f.is_dir())
    total_size += sum(f.stat().st_size for f in deps_files)
    
    # Get repository stats
    repo_counts = {}
    for repo_id, cache_keys in cache_manager._dependency_map.items():
        repo_counts[repo_id] = len(cache_keys)
    
    # Sort repositories by cache count
    top_repos = sorted(repo_counts.items(), key=lambda x: x[1], reverse=True)[:10]
    
    # Get performance metrics
    performance = cache_manager.get_metrics()
    
    return {
        "cache_count": cache_count,
        "deps_count": deps_count,
        "total_size_bytes": total_size,
        "total_size_mb": total_size / (1024 * 1024),
        "repos_tracked": len(repo_counts),
        "top_repos": top_repos,
        "performance": performance
    }

def print_cache_stats():
    """Print cache statistics."""
    stats = get_cache_stats()
    perf = stats['performance']
    
    print("\n===== CACHE STATISTICS =====")
    print(f"Total cache entries: {stats['cache_count']}")
    print(f"Dependency mappings: {stats['deps_count']}")
    print(f"Total size: {stats['total_size_mb']:.2f} MB")
    print(f"Repositories tracked: {stats['repos_tracked']}")
    
    print("\n----- Performance Metrics -----")
    print(f"Cache hits: {perf['hits']}")
    print(f"Cache misses: {perf['misses']}")
    print(f"Total requests: {perf['hits'] + perf['misses']}")
    
    hit_rate = perf['hits'] / (perf['hits'] + perf['misses']) * 100 if perf['hits'] + perf['misses'] > 0 else 0
    print(f"Hit rate: {hit_rate:.1f}%")
    
    print(f"Cache sets: {perf['sets']}")
    print(f"Cache invalidations: {perf['invalidations']}")
    
    efficiency = (perf['hits'] - perf['invalidations']) / perf['sets'] * 100 if perf['sets'] > 0 else 0
    print(f"Efficiency: {efficiency:.1f}%")
    
    print(f"Metrics since: {perf.get('start_time', 'Unknown')}")
    
    if stats['top_repos']:
        print("\n----- Top Repositories by Cache Entries -----")
        for repo_id, count in stats['top_repos']:
            print(f"  {repo_id}: {count} entries")
    
    print("============================\n")

def main():
    """Main entry point for the script."""
    parser = argparse.ArgumentParser(description="Cache management utility for awesome-virome")
    
    # Add command line arguments
    parser.add_argument('--clear-all', action='store_true', help='Clear all caches')
    parser.add_argument('--clear-repo', type=str, help='Clear caches for a specific repository URL')
    parser.add_argument('--stats', action='store_true', help='Show cache statistics')
    
    args = parser.parse_args()
    
    # No arguments provided, show help
    if not any(vars(args).values()):
        parser.print_help()
        return
    
    # Show statistics
    if args.stats:
        print_cache_stats()
    
    # Clear specific repository caches
    if args.clear_repo:
        clear_repo_caches(args.clear_repo)
    
    # Clear all caches
    if args.clear_all:
        clear_all_caches()

if __name__ == "__main__":
    main()