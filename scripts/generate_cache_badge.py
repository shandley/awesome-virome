#!/usr/bin/env python3
"""
Generate cache status badge for shields.io endpoint.

This script creates a JSON file in the .github/badges directory that shields.io
can use to display a badge showing the current cache health status.
"""

import os
import sys
import json
import argparse
from pathlib import Path

# Try to import the enhanced cache monitor
try:
    sys.path.insert(0, os.path.dirname(__file__))
    from apis.enhanced_cache import get_enhanced_cache_manager
    from enhanced_monitor_cache import EnhancedCacheMonitor
    HAS_CACHE_MONITOR = True
except ImportError:
    HAS_CACHE_MONITOR = False
    print("Error: Enhanced cache monitor not found.")
    sys.exit(1)

# Constants
BADGE_DIR = Path(os.path.join(".github", "badges"))
BADGE_FILE = BADGE_DIR / "cache-status.json"
CACHE_DIR = Path(os.path.join("metadata", "cache"))

def get_cache_health():
    """
    Get the current cache health by running a monitoring snapshot.
    
    Returns:
        tuple: (status, metrics) where status is GOOD, FAIR, or POOR
              and metrics is a dict of key cache metrics
    """
    # Get the cache manager
    cache_manager = get_enhanced_cache_manager(CACHE_DIR)
    
    # Create the monitor
    monitor = EnhancedCacheMonitor(cache_manager)
    
    # Take a snapshot and analyze
    snapshot, analysis = monitor.monitor_once(report=False, graphs=False)
    
    # Get the health status
    status = analysis['health'].upper()
    
    # Extract key metrics
    metrics = {
        'entries': snapshot['metrics'].get('cache_files', 0),
        'size_mb': snapshot['metrics'].get('cache_size_mb', 0),
        'size_percentage': snapshot['metrics'].get('size_percentage', 0),
        'hit_rate': snapshot['metrics'].get('hit_rate', 0) * 100,
        'efficiency': snapshot['metrics'].get('cache_efficiency', 0) * 100,
        'warnings': len(analysis['warnings']),
        'recommendations': len(analysis['recommendations'])
    }
    
    return status, metrics

def generate_badge(status, metrics):
    """
    Generate a badge JSON file for shields.io endpoint.
    
    Args:
        status: GOOD, FAIR, or POOR
        metrics: Dict of cache metrics
    
    Returns:
        Path to the generated badge file
    """
    # Create the badge directory if it doesn't exist
    BADGE_DIR.mkdir(exist_ok=True, parents=True)
    
    # Map status to color
    color_map = {
        'GOOD': 'brightgreen',
        'FAIR': 'yellow',
        'POOR': 'red',
    }
    
    color = color_map.get(status, 'lightgrey')
    
    # Format badge data
    badge_data = {
        'schemaVersion': 1,
        'label': 'Cache',
        'message': status,
        'color': color
    }
    
    # Add cache metrics as custom fields (not used by shields.io but useful for reference)
    for key, value in metrics.items():
        if isinstance(value, float):
            badge_data[key] = f"{value:.1f}"
        else:
            badge_data[key] = value
    
    # Write the badge file
    with open(BADGE_FILE, 'w') as f:
        json.dump(badge_data, f, indent=2)
    
    return BADGE_FILE

def main():
    """Main entry point for the script."""
    parser = argparse.ArgumentParser(description="Generate cache status badge")
    parser.add_argument('--output', type=str, default=str(BADGE_FILE),
                        help=f'Output file path (default: {BADGE_FILE})')
    
    args = parser.parse_args()
    
    try:
        # Get cache health
        status, metrics = get_cache_health()
        
        # Generate badge
        badge_file = generate_badge(status, metrics)
        
        print(f"Generated cache status badge at {badge_file}")
        print(f"Status: {status}")
        print(f"Hit rate: {metrics['hit_rate']:.1f}%")
        print(f"Efficiency: {metrics['efficiency']:.1f}%")
        print(f"Warnings: {metrics['warnings']}")
        
        # For GitHub Actions, set output variables
        if os.environ.get('GITHUB_ACTIONS') == 'true':
            with open(os.environ.get('GITHUB_ENV', ''), 'a') as env_file:
                env_file.write(f"CACHE_STATUS={status}\n")
                env_file.write(f"CACHE_COLOR={badge_data['color']}\n")
    
    except Exception as e:
        print(f"Error generating badge: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()