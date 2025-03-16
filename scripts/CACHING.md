# Smart Caching System for awesome-virome

This document describes the smart caching system implemented for the awesome-virome repository's metadata collection tools.

## Overview

The smart caching system provides:

1. **Dependency Tracking**: Tracks which cached data is associated with which repositories
2. **Selective Invalidation**: When a repository is updated, all its cached data is automatically invalidated
3. **Improved Performance**: Reduces API calls to external services like GitHub, Semantic Scholar, etc.
4. **Reduced Rate Limiting**: Helps avoid hitting API rate limits by efficiently using cached data
5. **Performance Metrics**: Tracks cache hits, misses, and efficiency metrics
6. **Rate Limit Awareness**: Dynamically adjusts behavior based on remaining API quota

## Key Components

### `CacheManager` Class

The core of the smart caching system is the `CacheManager` class in `scripts/apis/citations_api.py`. This class:

- Manages cache storage and retrieval
- Maintains dependency maps between repositories and cached data
- Provides methods for invalidating caches based on repository updates
- Validates cache freshness using configurable expiration policies
- Tracks performance metrics like hits, misses, and invalidations
- Calculates cache efficiency statistics

### `RateLimiter` Class

The enhanced rate limiter in `scripts/apis/citations_api.py` provides:

- Basic rate limiting to space out API requests
- Dynamic adjustment based on API rate limit headers
- Intelligent backoff when approaching rate limits
- Automatic retry logic for rate limit errors
- Header parsing for various API rate limit formats

### Integration Points

The caching system is integrated at several key points:

1. **Repository Metadata Updates**: In `enhance_metadata.py`, when repository metadata is updated, associated caches are automatically invalidated
2. **API Clients**: All API clients (GitHub, Semantic Scholar, etc.) use the cache manager for data storage and retrieval
3. **Academic Impact Collection**: The academic impact collector uses the cache manager to efficiently store and retrieve citation data

### Utilities

A utility script (`clear_cache.py`) is provided for manual cache management:

- View cache statistics
- Clear all caches
- Clear caches for specific repositories

## Usage

### API Methods

The cache manager provides these core methods:

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

### Command-line Usage

Cache management via command line:

```bash
# Show cache statistics and performance metrics
python scripts/clear_cache.py --stats

# Clear all caches
python scripts/clear_cache.py --clear-all

# Clear caches for a specific repository
python scripts/clear_cache.py --clear-repo "https://github.com/username/repo"
```

Example statistics output:

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
Cache sets: 98
Cache invalidations: 12
Efficiency: 146.9%

----- Top Repositories by Cache Entries -----
  5a8d7e2b1c9f: 12 entries
  3f7b2a6c9d4e: 8 entries
  ...
============================
```

## How Repository Dependency Tracking Works

1. When data is cached, it can be associated with a repository URL
2. Internally, the cache manager creates a dependency mapping between the repository and the cache key
3. When a repository is updated, the cache manager looks up all dependent cache entries and invalidates them
4. Dependencies are stored persistently, so they survive across different script runs

## Performance Considerations

- Metadata about repository dependencies is stored in the `_dependencies` subdirectory
- Performance metrics are stored in the `_metrics` subdirectory
- Cached data is hashed and stored in the base cache directory
- Cache expiration is based on file modification times for efficiency
- The system minimizes disk I/O by only loading mappings when needed

## Rate Limit Awareness

The enhanced rate limiter:

1. Extracts rate limit information from API response headers
2. Tracks remaining quota and reset times
3. Dynamically adjusts request spacing when quota is low
4. Implements smart retry logic for rate-limited requests
5. Supports different header formats used by various APIs

## Testing

The caching system includes comprehensive tests to verify functionality:

- **Unit Tests**: Located in `scripts/tests/test_cache_manager.py`
- **GitHub Actions**: Workflow in `.github/workflows/test_caching.yml` runs tests in CI
- **Performance Metrics**: Measured during actual usage to verify effectiveness

To run the tests:

```bash
# Run unit tests
python -m pytest scripts/tests/test_cache_manager.py -v

# Run with coverage measurement
python -m pytest scripts/tests/test_cache_manager.py -v --cov=scripts/apis/citations_api
```

See `scripts/tests/README.md` for more details on testing.