# awesome-virome Scripts

This directory contains scripts used for maintaining and enhancing the awesome-virome repository.

## Scripts Overview

- `apis/` - API integration modules for GitHub, Semantic Scholar, Zenodo, CrossRef, etc.
- `enhance_metadata.py` - Enhance repository metadata with additional information
- `bioinformatics_metadata.py` - Generate bioinformatics-specific metadata
- `academic_impact.py` - Collect academic impact metrics for tools
- `clear_cache.py` - Cache management utilities
- `comprehensive_citation_data.py` - Authoritative script for generating impact_data.json
- `auto_fix_dois.py` - Fixes DOI format issues across repositories
- `monitor_cache.py` - Advanced cache monitoring system
- `cron_cache_monitor.sh` - Scheduled cache monitoring for production environments
- `test_cache_monitoring.sh` - Test script for the cache monitoring system
- `citation_report.py` - Generate citation reports for tools

## Citation Data System

The awesome-virome repository uses a streamlined citation data collection system:

1. **Data Collection** - Citation data is collected from PubMed, CrossRef, and other sources
2. **DOI Validation** - DOI format is validated and inconsistencies are fixed
3. **Report Generation** - Citation reports are generated for analysis
4. **Impact Data** - All citation data is consolidated in impact_data.json

### Citation Data Scripts

- **comprehensive_citation_data.py** - The authoritative script for generating `impact_data.json`
  - Collects data from multiple sources: academic_impact, bioinformatics, and citation reports
  - Uses only real citation data with no synthetic data generation
  - Consolidates all available citation information

- **auto_fix_dois.py** - The authoritative script for DOI format fixing
  - Used by the `fix-dois.yml` workflow
  - Automatically fixes common DOI format issues

- **update_citation_counts.py** - Updates citation counts from CrossRef API
- **pubmed_citations.py** - Collects citation data from PubMed
- **citation_report.py** - Generates reports based on citation data
- **validate_citations.py** - Validates DOI consistency across the repository

## Cache System

The awesome-virome repository uses a smart caching system that provides:

1. **Smart Invalidation** - Automatic invalidation of cache entries when repositories change
2. **Dependency Tracking** - Track which cached data depends on which repositories
3. **Rate Limit Awareness** - Adaptive behavior based on API rate limits
4. **Performance Metrics** - Tracking of cache efficiency and performance

For more details, see [CACHING.md](CACHING.md).

## Cache Monitoring System

The cache monitoring system provides real-time and historical metrics on cache performance, with features for alerting, visualization, and optimization recommendations.

Features include:
- Real-time performance monitoring
- Historical trend analysis
- Cache health checks
- Performance visualization
- CSV export of metrics
- Scheduled monitoring via cron

For more information, see [CACHE_MONITORING.md](CACHE_MONITORING.md).

## Usage

### Metadata Enhancement

```bash
# Enhance metadata for all repositories
python enhance_metadata.py

# Enhance metadata with additional bioinformatics information
python bioinformatics_metadata.py

# Collect academic impact information
python academic_impact.py
```

### Cache Management

```bash
# View cache statistics
python clear_cache.py --stats

# Clear all caches
python clear_cache.py --clear-all

# Clear caches for a specific repository
python clear_cache.py --clear-repo "https://github.com/username/repo"
```

### Cache Monitoring

```bash
# Run basic monitoring
python monitor_cache.py

# Run continuous monitoring for 1 hour
python monitor_cache.py --continuous --duration 3600

# Generate performance graphs
python monitor_cache.py --graphs

# Export metrics to CSV
python monitor_cache.py --export-csv
```

### Scheduled Monitoring

```bash
# Run a quick snapshot
./cron_cache_monitor.sh

# Run hourly monitoring
./cron_cache_monitor.sh hourly

# Run daily monitoring
./cron_cache_monitor.sh daily

# Run weekly comprehensive monitoring
./cron_cache_monitor.sh weekly
```

### Testing

```bash
# Run tests for the cache monitoring system
./test_cache_monitoring.sh
```

## Development

When contributing to this repository, please ensure all scripts maintain compatibility with Python 3.6+ and follow the repository's code style.