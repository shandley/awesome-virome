# Testing the Awesome-Virome Cache System

This directory contains tests for the various components of the Awesome-Virome caching system.

## Running the Tests

### Unit Tests

To run the unit tests for the CacheManager:

```bash
# Run from the repository root
python -m pytest scripts/tests/test_cache_manager.py -v

# With coverage
python -m pytest scripts/tests/test_cache_manager.py -v --cov=scripts/apis/citations_api
```

### Integration Testing

To test the cache system in a live environment, you can use:

```bash
# Run the enhance_metadata.py script with a small limit
python scripts/enhance_metadata.py --limit 5

# Check the cache statistics
python scripts/clear_cache.py --stats

# Clear a specific repository's cache
python scripts/clear_cache.py --clear-repo "https://github.com/user/repo"

# Clear all caches
python scripts/clear_cache.py --clear-all
```

## Cache Features

The caching system includes the following features:

1. **Smart Invalidation**: Cache items are only invalidated when their associated repositories change, preventing unnecessary API requests

2. **Rate Limit Aware Caching**: The system automatically adjusts cache expiration based on remaining API rate limits

3. **Performance Metrics**: Tracks hits, misses, invalidations, and efficiency of the cache

## Adding New Tests

When adding new API clients or cache features, please add corresponding tests that verify:

1. Cache hits and misses are properly handled
2. Cache dependencies are correctly tracked
3. Cache invalidation works as expected
4. Rate limiting is respected

## GitHub Actions

The GitHub Actions workflow in `.github/workflows/test_caching.yml` runs both the unit tests and a simplified integration test to ensure the caching system works as expected.