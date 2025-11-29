# Enhanced Cache Management System

This document describes the enhanced cache management system for the awesome-virome repository.

## Overview

The enhanced cache system provides all the features of the original system plus:

1. **Automatic Expiration**: Automatically removes cache entries older than the specified expiry time
2. **LRU Eviction**: Maintains a "least recently used" tracking system to intelligently remove stale entries
3. **Size-Based Limits**: Enforces a maximum cache size with automatic eviction when limits are approached
4. **Performance Metrics**: Tracks detailed metrics including eviction counts and effectiveness
5. **Background Maintenance**: Runs cache cleanup in a background thread to maintain performance
6. **Advanced Monitoring**: Provides comprehensive monitoring tools with visualizations

## Configuration

The enhanced cache system is configured with these parameters:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `cache_dir` | "metadata/cache" | Base directory for cache storage |
| `expiry_days` | 30 | Number of days before entries automatically expire |
| `max_size_mb` | 500 | Maximum cache size in MB |
| `lru_check_interval` | 3600 | How often to check LRU status (seconds) |
| `size_limit_percentage` | 0.9 | When to start cleanup (% of max size) |

## How It Works

### Automatic Expiration

Cache entries older than `expiry_days` are automatically removed during maintenance cycles. This ensures the cache doesn't accumulate stale data.

### LRU (Least Recently Used) Tracking

The system tracks when each cache entry was last accessed. When cache cleanup is needed, the least recently used entries are removed first. This preserves frequently accessed data while removing entries that haven't been used recently.

### Size-Based Eviction

When the cache size exceeds `max_size_mb * size_limit_percentage`, the system automatically removes enough least recently used entries to bring the size down to 80% of the maximum.

### Background Maintenance

The cache system runs a background thread that periodically:
- Checks for and removes expired entries
- Performs LRU-based cleanup if needed
- Ensures the cache stays within size limits

## Directory Structure

```
metadata/cache/
  |- *.json                   # Cache data files
  |- _dependencies/           # Repository dependency mappings
  |- _metrics/                # Performance metrics
  |- _lru/                    # LRU tracking data
  |- _monitoring/             # Monitoring history and graphs
      |- enhanced/            # Enhanced monitoring data
          |- graphs/          # Performance visualizations
```

## Monitoring Tools

### Enhanced Monitor Script

The `enhanced_monitor_cache.py` script provides comprehensive monitoring:

```bash
# Basic monitoring
python scripts/enhanced_monitor_cache.py

# Continuous monitoring
python scripts/enhanced_monitor_cache.py --continuous --interval 300

# Generate performance graphs
python scripts/enhanced_monitor_cache.py --graphs

# Export metrics to CSV
python scripts/enhanced_monitor_cache.py --export-csv
```

### Maintenance Script

The `cache_maintenance.sh` script handles regular maintenance:

```bash
# Daily maintenance
./scripts/cache_maintenance.sh daily

# Weekly comprehensive maintenance
./scripts/cache_maintenance.sh weekly

# Quick snapshot
./scripts/cache_maintenance.sh snapshot
```

## Migrating from the Old System

The `migrate_cache.py` script handles migration from the old cache system:

```bash
# Migrate with default settings
python scripts/migrate_cache.py

# Migrate with custom max size
python scripts/migrate_cache.py --max-size 1000
```

## GitHub Actions Integration

The repository includes a GitHub Actions workflow (`cache-maintenance.yml`) that automatically runs cache maintenance:

- Daily maintenance runs every day
- Weekly comprehensive maintenance runs on Sundays
- Generates and uploads monitoring graphs
- Updates a cache status badge

## Performance Impact

The enhanced cache system improves overall performance by:

- Reducing the number of API requests through better cache retention
- Preventing unbounded cache growth
- Automatically refreshing stale cache entries
- Optimizing disk usage by removing unused entries

## Troubleshooting

If you encounter issues with the cache system:

1. Check the maintenance logs in the `logs/` directory
2. Run a snapshot monitoring to check cache health:
   ```bash
   python scripts/enhanced_monitor_cache.py --summary
   ```
3. If needed, manually clear the cache:
   ```bash
   python scripts/clear_cache.py --clear-all
   ```
4. Restore from a backup if available (backups are created during migration)