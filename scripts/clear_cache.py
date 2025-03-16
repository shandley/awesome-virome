#!/usr/bin/env python3
"""
Utility script for managing the cache in the awesome-virome repository.

This script provides functionality to clear caches, including:
- Clearing all caches
- Clearing caches related to specific repositories
- Listing cache statistics
- Running advanced cache monitoring
"""

import os
import sys
import argparse
import logging
import subprocess
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

# Script paths
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
MONITOR_SCRIPT = os.path.join(SCRIPT_DIR, "monitor_cache.py")
CRON_MONITOR_SCRIPT = os.path.join(SCRIPT_DIR, "cron_cache_monitor.sh")

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

def run_monitoring(args_list):
    """Run the cache monitoring script with the given arguments."""
    if not os.path.exists(MONITOR_SCRIPT):
        logger.error(f"Cache monitoring script not found at {MONITOR_SCRIPT}")
        return 1
    
    cmd = [sys.executable, MONITOR_SCRIPT] + args_list
    logger.info(f"Running monitoring: {' '.join(cmd)}")
    
    return subprocess.call(cmd)

def run_scheduled_monitoring(mode):
    """Run the scheduled monitoring script with the given mode."""
    if not os.path.exists(CRON_MONITOR_SCRIPT):
        logger.error(f"Scheduled monitoring script not found at {CRON_MONITOR_SCRIPT}")
        return 1
    
    cmd = [CRON_MONITOR_SCRIPT, mode]
    logger.info(f"Running scheduled monitoring ({mode}): {' '.join(cmd)}")
    
    return subprocess.call(cmd)

def main():
    """Main entry point for the script."""
    parser = argparse.ArgumentParser(description="Cache management utility for awesome-virome")
    
    # Cache management arguments
    group1 = parser.add_argument_group('Cache Management')
    group1.add_argument('--clear-all', action='store_true', help='Clear all caches')
    group1.add_argument('--clear-repo', type=str, help='Clear caches for a specific repository URL')
    group1.add_argument('--stats', action='store_true', help='Show cache statistics')
    
    # Monitoring arguments
    group2 = parser.add_argument_group('Cache Monitoring')
    group2.add_argument('--monitor', action='store_true', help='Run basic cache monitoring')
    group2.add_argument('--continuous', action='store_true', help='Run continuous monitoring')
    group2.add_argument('--duration', type=int, help='Duration of continuous monitoring in seconds')
    group2.add_argument('--interval', type=int, default=60, help='Monitoring interval in seconds')
    group2.add_argument('--graphs', action='store_true', help='Generate performance graphs')
    group2.add_argument('--export-csv', action='store_true', help='Export metrics to CSV')
    group2.add_argument('--csv-path', type=str, help='Path for CSV export')
    
    # Scheduled monitoring arguments
    group3 = parser.add_argument_group('Scheduled Monitoring')
    group3.add_argument('--scheduled', choices=['snapshot', 'hourly', 'daily', 'weekly'], 
                      help='Run scheduled monitoring with specified mode')
    
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
    
    # Run monitoring
    if args.monitor or args.continuous or args.graphs or args.export_csv:
        monitor_args = []
        
        if args.continuous:
            monitor_args.append('--continuous')
        
        if args.duration:
            monitor_args.extend(['--duration', str(args.duration)])
        
        if args.interval:
            monitor_args.extend(['--interval', str(args.interval)])
        
        if args.graphs:
            monitor_args.append('--graphs')
        
        if args.export_csv:
            monitor_args.append('--export-csv')
            
        if args.csv_path:
            monitor_args.extend(['--csv-path', args.csv_path])
        
        run_monitoring(monitor_args)
    
    # Run scheduled monitoring
    if args.scheduled:
        run_scheduled_monitoring(args.scheduled)

if __name__ == "__main__":
    main()