# Cache System Documentation

Documentation for the awesome-virome API-response cache used during metadata collection.

## Overview

The cache stores API responses (GitHub, CrossRef, and others) so repeated metadata runs do not re-fetch data that has not changed. It provides:

1. Dependency tracking: links cached data to specific repositories for targeted invalidation
2. Selective invalidation: invalidates related entries when a repository updates
3. Automatic expiration: removes entries older than the configured expiry time
4. LRU eviction: removes least recently used entries when the size limit is approached
5. Size-based limits: enforces a maximum cache size with automatic cleanup
6. Rate-limit awareness: adjusts request spacing based on API quota headers

## Core components

### CacheManager

Located in `scripts/apis/citations_api.py`. Handles storage and retrieval, dependency mapping between repositories and cached data, invalidation on repository updates, and freshness validation with configurable expiration. It is imported by the metadata scripts (`enhance_metadata.py`, `bioinformatics_metadata.py`, `update_check.py`).

### RateLimiter

Also in `scripts/apis/citations_api.py`. Spaces out API requests, reads rate-limit headers to adjust dynamically, backs off as quota runs low, and retries rate-limited requests. Supports the GitHub and CrossRef header formats.

## Configuration

| Parameter | Default | Description |
|-----------|---------|-------------|
| `cache_dir` | `metadata/cache` | Base directory for cache storage |
| `expiry_days` | 30 | Days before entries expire |
| `max_size_mb` | 500 | Maximum cache size in MB |
| `lru_check_interval` | 3600 | LRU check frequency (seconds) |
| `size_limit_percentage` | 0.9 | Cleanup trigger (fraction of max size) |

Directory layout:

```
metadata/cache/
  |- *.json            # cache data files
  |- _dependencies/    # repository dependency mappings
  |- _metrics/         # performance metrics
  |- _lru/             # LRU tracking data
```

## Usage

### API methods

```python
cache_manager.set(cache_key, data, repo_url=optional_repo_url)   # store
data = cache_manager.get(cache_key)                              # retrieve
is_valid = cache_manager.is_valid(cache_key)                     # freshness check
cache_manager.invalidate_repo_caches(repo_url)                   # invalidate a repo's entries
cache_manager.invalidate(cache_key)                             # invalidate one entry
cache_manager.clear_all()                                       # clear everything
```

### Command-line tools

```bash
# Show cache statistics
python scripts/clear_cache.py --stats

# Clear all caches
python scripts/clear_cache.py --clear-all

# Clear caches for a specific repository
python scripts/clear_cache.py --clear-repo "https://github.com/username/repo"

# Pre-warm the cache for predictable access patterns
python scripts/cache_warming.py
```

## How invalidation works

- When data is cached with a `repo_url`, a dependency is recorded under `metadata/cache/_dependencies/`.
- When a repository updates, its dependent cache entries are invalidated.
- Access times are tracked under `metadata/cache/_lru/`; on cleanup, the oldest entries are removed first.

## Testing

```bash
python -m pytest scripts/tests/test_cache_manager.py -v
```

---

*History: this document replaced `CACHING.md`, `CACHE_MANAGEMENT.md`, and `CACHE_MONITORING.md` (archived November 2025). The separate cache-monitoring and badge tooling was removed in July 2026 along with the cache-maintenance workflow.*
