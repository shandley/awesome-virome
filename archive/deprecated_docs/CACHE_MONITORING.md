# Cache Monitoring System for awesome-virome

This document describes the cache monitoring system implemented for the awesome-virome repository's smart caching infrastructure.

## Overview

The cache monitoring system provides:

1. **Real-time Performance Metrics**: Collects and analyzes cache performance in real-time
2. **Historical Trend Analysis**: Tracks changes in cache metrics over time
3. **Automatic Health Checks**: Evaluates cache health and identifies potential issues
4. **Visualization**: Generates graphs of key performance indicators
5. **Scheduled Monitoring**: Supports automated monitoring via cron
6. **Reporting**: Exports metrics in various formats (JSON, CSV)
7. **Alerting**: Detects concerning metrics and can trigger alerts

## Key Components

### `monitor_cache.py` Script

The core monitoring tool that provides:

- On-demand cache metrics snapshots
- Continuous monitoring with configurable intervals
- Historical metrics tracking and comparison
- Cache health analysis with recommendations
- Performance graph generation
- Metrics export

### `cron_cache_monitor.sh` Script

A shell script for scheduled monitoring:

- Designed to be run by cron for regular monitoring
- Supports different monitoring modes (snapshot, hourly, daily, weekly)
- Includes automatic cleanup of old reports
- Can trigger alerts when cache health issues are detected

### History and Metrics Storage

- Metrics are stored in `metadata/cache/_monitoring/`
- Historical metrics enable tracking changes over time
- Time series data allows for trend analysis
- Performance graphs visualize key metrics

## Monitoring Metrics

The system tracks these key metrics:

- **Hit Rate**: Percentage of cache requests that are served from cache
- **Miss Rate**: Percentage of requests that require fetching fresh data
- **Efficiency**: Measure of how effectively the cache is being used
- **Invalidations**: Number of cache entries that have been invalidated
- **Cache Size**: Total number of entries in the cache
- **Growth Rate**: How quickly the cache is growing over time
- **Repository Dependencies**: How many cache entries depend on each repository

## Usage

### Basic Monitoring

Run a one-time cache check with detailed report:

```bash
python scripts/monitor_cache.py
```

### Continuous Monitoring

Monitor the cache continuously for a specified duration:

```bash
# Monitor for 1 hour (3600 seconds) with 60-second intervals
python scripts/monitor_cache.py --continuous --duration 3600 --interval 60
```

### Generate Performance Graphs

Create visualizations of cache performance:

```bash
python scripts/monitor_cache.py --graphs
```

### Export Metrics to CSV

Export historical metrics to CSV for further analysis:

```bash
python scripts/monitor_cache.py --export-csv --csv-path /path/to/metrics.csv
```

### Scheduled Monitoring with Cron

Set up automated monitoring using the cron script:

```bash
# Run a quick snapshot
./scripts/cron_cache_monitor.sh

# Run hourly monitoring with graphs
./scripts/cron_cache_monitor.sh hourly

# Run daily monitoring with extended metrics
./scripts/cron_cache_monitor.sh daily

# Run weekly monitoring with comprehensive analysis
./scripts/cron_cache_monitor.sh weekly
```

Example crontab entries:

```
# Run hourly monitoring during work hours
0 9-17 * * 1-5 /path/to/awesome-virome/scripts/cron_cache_monitor.sh hourly

# Run daily summary at midnight
0 0 * * * /path/to/awesome-virome/scripts/cron_cache_monitor.sh daily

# Run weekly comprehensive report on Sunday
0 1 * * 0 /path/to/awesome-virome/scripts/cron_cache_monitor.sh weekly
```

## Health Analysis

The monitoring system analyzes cache health based on multiple factors:

- Hit rate below 30% is considered poor
- Hit rate between 30-50% is considered fair
- Hit rate above 50% is considered good
- Negative cache efficiency is a critical warning
- High invalidation rates may indicate excessive dependency invalidation
- Rapidly declining hit rates suggest cache configuration issues

When health issues are detected, the system provides specific recommendations for improving cache performance.

## Visualization and Reporting

The monitoring system generates several types of visualizations:

1. **Hit Rate vs Miss Rate**: Shows the balance between cache hits and misses
2. **Cache Efficiency**: Tracks how efficiently the cache is being utilized
3. **Invalidations**: Shows cache invalidation patterns over time
4. **Cache Size**: Tracks growth in the number of cached items

Reports include:

- Current metrics snapshot
- Change compared to baseline
- Health assessment
- Specific warnings for potential issues
- Actionable recommendations to improve performance
- Insights on positive aspects of cache performance

## Interpreting Metrics

### Hit Rate

- **Good**: >70% - Cache is working effectively
- **Fair**: 50-70% - Reasonable performance but room for improvement
- **Poor**: <50% - Cache may need configuration adjustments

### Efficiency

Calculated as: `(hits - invalidations) / sets * 100%`

- **Good**: >70% - Cache entries are being reused effectively
- **Fair**: 30-70% - Moderate reuse with some wasted entries
- **Poor**: <30% - Many cached entries are being invalidated before reuse

### Invalidation Rate

High invalidation rates (>10/hour) may indicate:

- Over-aggressive repository dependency tracking
- Too many dependencies between repositories and cache entries
- Frequent repository updates causing cascade invalidations

## Best Practices

1. **Regular Monitoring**: Run daily monitoring to track performance trends
2. **Baseline Comparison**: Compare current metrics to historical baselines
3. **Act on Recommendations**: Implement the recommendations provided by the analyzer
4. **Adjust Cache Expiration**: Tune expiration times based on hit rate patterns
5. **Review Dependencies**: Optimize repository dependency tracking if invalidation rates are high

## Troubleshooting

### Low Hit Rate

If hit rate is consistently low:

- Check if cache expiration times are too short
- Ensure repository dependency tracking is correctly configured
- Verify that API clients are using the cache manager properly

### High Invalidation Rate

If invalidation rate is excessive:

- Review repository dependency mapping logic
- Consider more selective dependency tracking
- Check for repository update patterns causing cascading invalidations

### Disk Space Concerns

If cache size is growing rapidly:

- Implement more aggressive expiration for less important data
- Consider selective cache clearing for oldest entries
- Review which API responses are being cached