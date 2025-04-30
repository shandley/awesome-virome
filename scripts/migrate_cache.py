#!/usr/bin/env python3
"""
Migration script to convert from the old cache system to the enhanced one.
This script preserves all existing cache data while enabling the new functionality.
"""

import os
import sys
import time
import json
import shutil
import logging
import argparse
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Import the cache managers
try:
    sys.path.insert(0, os.path.join(os.path.dirname(__file__)))
    from apis.citations_api import cache_manager as old_cache_manager
    from apis.enhanced_cache import get_enhanced_cache_manager
    HAS_CACHE_MANAGERS = True
except ImportError as e:
    HAS_CACHE_MANAGERS = False
    logger.error(f"Error importing cache managers: {e}")
    logger.error("Make sure apis/citations_api.py and apis/enhanced_cache.py exist")
    sys.exit(1)

def backup_cache_directory(cache_dir, backup_suffix="_backup"):
    """
    Create a backup of the cache directory.
    
    Args:
        cache_dir: Path to the cache directory
        backup_suffix: Suffix to add to the backup directory name
        
    Returns:
        Path to the backup directory
    """
    cache_dir = Path(cache_dir)
    backup_dir = Path(str(cache_dir) + backup_suffix)
    
    logger.info(f"Creating backup of cache directory: {backup_dir}")
    
    # Remove existing backup if it exists
    if backup_dir.exists():
        logger.warning(f"Removing existing backup directory: {backup_dir}")
        shutil.rmtree(backup_dir)
    
    # Create backup
    shutil.copytree(cache_dir, backup_dir)
    logger.info(f"Backup created at: {backup_dir}")
    
    return backup_dir

def migrate_cache(cache_dir, max_size_mb=500):
    """
    Migrate from the old cache system to the enhanced one.
    
    Args:
        cache_dir: Path to the cache directory
        max_size_mb: Maximum cache size in MB for the new cache system
        
    Returns:
        True if successful, False otherwise
    """
    # Create a backup first
    backup_dir = backup_cache_directory(cache_dir)
    
    try:
        # Create the enhanced cache manager
        enhanced_cache = get_enhanced_cache_manager(cache_dir, max_size_mb)
        
        # Create LRU directory structure
        lru_dir = Path(cache_dir) / "_lru"
        lru_dir.mkdir(exist_ok=True)
        
        # Get existing cache metrics
        old_metrics = old_cache_manager.get_metrics()
        
        # Add migration info to metrics
        migration_info = {
            "migration_date": time.strftime("%Y-%m-%d %H:%M:%S"),
            "backup_location": str(backup_dir),
            "original_cache_size": old_metrics.get('cache_files', 0),
            "max_size_mb": max_size_mb
        }
        
        # Save migration info
        with open(Path(cache_dir) / "_metrics" / "migration_info.json", 'w') as f:
            json.dump(migration_info, f, indent=2)
        
        logger.info(f"Migration completed successfully")
        logger.info(f"New cache system configured with {max_size_mb} MB maximum size")
        logger.info(f"Backup is available at: {backup_dir}")
        
        return True
        
    except Exception as e:
        logger.error(f"Migration failed: {e}")
        logger.error(f"Consider restoring from backup: {backup_dir}")
        return False

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Migrate from the old cache system to the enhanced one"
    )
    
    parser.add_argument(
        "--cache-dir",
        default=os.path.join("metadata", "cache"),
        help="Path to the cache directory (default: metadata/cache)"
    )
    
    parser.add_argument(
        "--max-size",
        type=int,
        default=500,
        help="Maximum cache size in MB (default: 500)"
    )
    
    parser.add_argument(
        "--no-backup",
        action="store_true",
        help="Skip creating a backup (not recommended)"
    )
    
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Simulate migration without making changes"
    )
    
    return parser.parse_args()

def main():
    """Main entry point."""
    args = parse_args()
    
    # Check if the cache directory exists
    cache_dir = Path(args.cache_dir)
    if not cache_dir.exists():
        logger.error(f"Cache directory does not exist: {cache_dir}")
        return 1
    
    # Dry run - just show what would happen
    if args.dry_run:
        logger.info(f"DRY RUN - No changes will be made")
        logger.info(f"Would migrate cache directory: {cache_dir}")
        logger.info(f"Would set maximum cache size to: {args.max_size} MB")
        return 0
    
    # Perform migration
    if args.no_backup:
        logger.warning("Skipping backup as requested (not recommended)")
        try:
            # Create the enhanced cache manager
            enhanced_cache = get_enhanced_cache_manager(str(cache_dir), args.max_size)
            
            # Create LRU directory structure
            lru_dir = cache_dir / "_lru"
            lru_dir.mkdir(exist_ok=True)
            
            logger.info(f"Migration completed successfully")
            return 0
        except Exception as e:
            logger.error(f"Migration failed: {e}")
            return 1
    else:
        if migrate_cache(str(cache_dir), args.max_size):
            return 0
        else:
            return 1

if __name__ == "__main__":
    sys.exit(main())