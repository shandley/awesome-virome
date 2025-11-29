# Cache System Documentation

Comprehensive documentation for the awesome-virome cache management and monitoring system.

## Table of Contents

- [Overview](#overview)
- [Core Components](#core-components)
- [Configuration](#configuration)
- [Usage](#usage)
- [Monitoring](#monitoring)
- [Maintenance](#maintenance)
- [Troubleshooting](#troubleshooting)
- [Best Practices](#best-practices)

## Overview

The awesome-virome cache system provides intelligent, efficient caching for API responses and metadata with:

### Key Features

1. **Dependency Tracking**: Links cached data to specific repositories for intelligent invalidation
2. **Selective Invalidation**: Automatically invalidates related cache entries when repositories update
3. **Automatic Expiration**: Removes entries older than the configured expiry time
4. **LRU Eviction**: Removes least recently used entries when size limits are approached
5. **Size-Based Limits**: Enforces maximum cache size with automatic cleanup
6. **Performance Metrics**: Tracks hits, misses, efficiency, and more
7. **Rate Limit Awareness**: Dynamically adjusts behavior based on API quotas
8. **Background Maintenance**: Runs cleanup automatically to maintain performance
9. **Comprehensive Monitoring**: Real-time metrics, trend analysis, and visualizations

## Core Components

### CacheManager Class

Located in `scripts/apis/citations_api.py`, provides:
- Cache storage and retrieval
- Dependency mapping between repositories and cached data
- Cache invalidation based on repository updates
- Freshness validation with configurable expiration
- Performance metrics tracking

### RateLimiter Class

Located in `scripts/apis/citations_api.py`, provides:
- Basic rate limiting to space out API requests
- Dynamic adjustment based on API rate limit headers
- Intelligent backoff when approaching limits
- Automatic retry logic for rate-limited requests
- Support for various API header formats

### Monitoring Tools

- **enhanced_monitor_cache.py**: Comprehensive monitoring with real-time metrics, trend analysis, and visualizations
- **cache_maintenance.sh**: Scheduled maintenance script for daily/weekly cleanup
- **clear_cache.py**: Manual cache management and statistics

## Configuration

### Cache Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `cache_dir` | "metadata/cache" | Base directory for cache storage |
| `expiry_days` | 30 | Days before entries automatically expire |
| `max_size_mb` | 500 | Maximum cache size in MB |
| `lru_check_interval` | 3600 | LRU check frequency (seconds) |
| `size_limit_percentage` | 0.9 | Cleanup trigger (% of max size) |

### Directory Structure

```
metadata/cache/
  |- *.json                   # Cache data files
  |- _dependencies/           # Repository dependency mappings
  |- _metrics/                # Performance metrics
  |- _lru/                    # LRU tracking data
  |- _monitoring/             # Monitoring history
      |- enhanced/            # Enhanced monitoring data
          |- graphs/          # Performance visualizations
```

## Usage

### API Methods

```python
# Store data in cache
cache_manager.set(cache_key, data, repo_url=optional_repo_url)

# Retrieve data from cache
data = cache_manager.get(cache_key)

# Check if cache is valid
is_valid = cache_manager.is_valid(cache_key)

# Invalidate all caches related to a repository
cache_manager.invalidate_repo_caches(repo_url)

# Invalidate a specific cache entry
cache_manager.invalidate(cache_key)

# Clear all caches
cache_manager.clear_all()
```

### Command-Line Usage

#### View Statistics

```bash
# Show cache statistics and performance metrics
python scripts/clear_cache.py --stats
```

Example output:
```
===== CACHE STATISTICS =====
Total cache entries: 42
Dependency mappings: 15
Total size: 2.34 MB
Repositories tracked: 10

----- Performance Metrics -----
Cache hits: 156
Cache misses: 58
Total requests: 214
Hit rate: 72.9%
Efficiency: 146.9%
============================
```

#### Clear Cache

```bash
# Clear all caches
python scripts/clear_cache.py --clear-all

# Clear caches for a specific repository
python scripts/clear_cache.py --clear-repo "https://github.com/username/repo"
```

## Monitoring

### Enhanced Monitoring

```bash
# Basic monitoring snapshot
python scripts/enhanced_monitor_cache.py

# Continuous monitoring
python scripts/enhanced_monitor_cache.py --continuous --interval 300

# Generate performance graphs
python scripts/enhanced_monitor_cache.py --graphs

# Export metrics to CSV
python scripts/enhanced_monitor_cache.py --export-csv

# Summary status (used in workflows)
python scripts/enhanced_monitor_cache.py --summary
```

### Scheduled Monitoring

Set up cron for automated monitoring:

```bash
# Hourly during work hours
0 9-17 * * 1-5 /path/to/scripts/cron_cache_monitor.sh hourly

# Daily summary at midnight
0 0 * * * /path/to/scripts/cron_cache_monitor.sh daily

# Weekly comprehensive on Sunday
0 1 * * 0 /path/to/scripts/cron_cache_monitor.sh weekly
```

### Key Metrics

- **Hit Rate**: % of requests served from cache
- **Miss Rate**: % of requests requiring fresh data
- **Efficiency**: `(hits - invalidations) / sets * 100%`
- **Invalidations**: Number of entries invalidated
- **Cache Size**: Total number of entries
- **Growth Rate**: Cache growth velocity

### Health Assessment

| Hit Rate | Status | Recommendation |
|----------|--------|----------------|
| >70% | Good | Cache working effectively |
| 50-70% | Fair | Monitor, room for improvement |
| <50% | Poor | Review configuration |

| Efficiency | Status | Meaning |
|------------|--------|---------|
| >70% | Good | Entries reused effectively |
| 30-70% | Fair | Moderate reuse |
| <30% | Poor | Many wasted entries |

## Maintenance

### Automated Maintenance

GitHub Actions workflow (`cache-maintenance.yml`):
- Daily maintenance every day
- Weekly comprehensive on Sundays
- Generates monitoring graphs
- Updates cache status badge

### Manual Maintenance

```bash
# Daily maintenance
./scripts/cache_maintenance.sh daily

# Weekly comprehensive
./scripts/cache_maintenance.sh weekly

# Quick snapshot
./scripts/cache_maintenance.sh snapshot
```

### Migration

Migrate from older cache systems:

```bash
# Migrate with defaults
python scripts/migrate_cache.py

# Migrate with custom max size
python scripts/migrate_cache.py --max-size 1000
```

Note: Backups are created automatically during migration.

## Troubleshooting

### Low Hit Rate

If hit rate is consistently low (<50%):

1. Check cache expiration times (may be too short)
2. Verify repository dependency tracking configuration
3. Ensure API clients use cache manager properly
4. Review invalidation patterns

### High Invalidation Rate

If invalidations are excessive (>10/hour):

1. Review repository dependency mapping logic
2. Consider more selective dependency tracking
3. Check for cascading invalidations from repo updates
4. Reduce unnecessary dependencies

### Disk Space Issues

If cache is growing too rapidly:

1. Implement more aggressive expiration for low-priority data
2. Use selective cache clearing for oldest entries
3. Review which API responses are being cached
4. Adjust `max_size_mb` and `size_limit_percentage`

### Quick Diagnostics

```bash
# Check cache health
python scripts/enhanced_monitor_cache.py --summary

# View detailed stats
python scripts/clear_cache.py --stats

# Check logs
tail -f logs/cache_monitor_*.log
```

## Best Practices

### Development

1. **Associate with repositories**: Always link cached data to repositories when applicable
2. **Validate before use**: Check `is_valid()` before using cached data
3. **Handle cache misses**: Always have a fallback for cache misses
4. **Use rate limiter**: Integrate rate limiter for all API calls

### Operations

1. **Regular monitoring**: Run daily monitoring to track trends
2. **Baseline comparison**: Compare current metrics to historical baselines
3. **Act on recommendations**: Implement health analysis recommendations
4. **Tune parameters**: Adjust expiration and size limits based on patterns
5. **Review dependencies**: Optimize repository dependency tracking

### Performance

1. **Smart invalidation**: Only invalidate what's truly dependent
2. **Batch operations**: Group cache operations when possible
3. **Monitor efficiency**: Track efficiency metric closely
4. **Prewarm cache**: Use `cache_warming.py` for predictable access patterns

## Advanced Topics

### Dependency Tracking

The system creates mappings between repositories and cache keys:

1. When data is cached with `repo_url`, a dependency is recorded
2. Dependencies persist across script runs
3. When a repository updates, all dependent caches are invalidated
4. Dependencies are stored in `metadata/cache/_dependencies/`

### Rate Limiting Strategy

The enhanced rate limiter:

1. Extracts rate limit info from API response headers
2. Tracks remaining quota and reset times
3. Dynamically adjusts spacing when quota is low
4. Implements smart retry with exponential backoff
5. Supports multiple header formats (GitHub, CrossRef, etc.)

### LRU Eviction

Least Recently Used tracking:

1. Access times are recorded for each cache entry
2. When cleanup is needed, oldest entries are removed first
3. Frequently accessed data is preserved longer
4. LRU data stored in `metadata/cache/_lru/`

## Testing

### Unit Tests

```bash
# Run cache manager tests
python -m pytest scripts/tests/test_cache_manager.py -v

# With coverage
python -m pytest scripts/tests/test_cache_manager.py --cov=scripts/apis/citations_api -v
```

### Integration Tests

GitHub Actions workflow `test_caching.yml` runs automated tests in CI.

---

*For questions or issues, see `.github/workflows/README.md` or check the logs in `logs/` directory.*

**Consolidated Documentation**: This document replaces `CACHING.md`, `CACHE_MANAGEMENT.md`, and `CACHE_MONITORING.md` (archived November 2025).
