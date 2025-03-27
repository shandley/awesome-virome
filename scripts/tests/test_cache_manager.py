#!/usr/bin/env python3
"""
Unit tests for the CacheManager class in citations_api.py
"""
import os
import json
import tempfile
import unittest
import shutil
from pathlib import Path
from unittest.mock import patch, MagicMock
from datetime import datetime, timedelta

# Add parent directory to path
import sys
script_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(script_dir)

# Import the module explicitly for coverage analysis
import apis.citations_api
from apis.citations_api import CacheManager

class TestCacheManager(unittest.TestCase):
    """Test cases for the CacheManager class."""
    
    def setUp(self):
        """Set up a temporary directory for cache storage."""
        self.temp_dir = tempfile.mkdtemp()
        self.cache_manager = CacheManager(self.temp_dir, expiry_days=1)
        
        # Test data
        self.test_data = {"key": "value", "nested": {"data": 123}}
        self.repo_url = "https://github.com/example/repo"
        self.cache_key = "test_cache_key"
        
    def tearDown(self):
        """Clean up temporary directory."""
        shutil.rmtree(self.temp_dir)
        
    def test_set_and_get(self):
        """Test setting and retrieving cache data."""
        # Set data in cache
        result = self.cache_manager.set(self.cache_key, self.test_data)
        self.assertTrue(result)
        
        # Get data from cache
        cached_data = self.cache_manager.get(self.cache_key)
        self.assertEqual(cached_data, self.test_data)
        
    def test_cache_expiration(self):
        """Test that expired caches are properly invalidated."""
        # Set data in cache
        result = self.cache_manager.set(self.cache_key, self.test_data)
        self.assertTrue(result)
        
        # Verify it's valid
        self.assertTrue(self.cache_manager.is_valid(self.cache_key))
        
        # Mock datetime.now to return a future date beyond expiry
        future_time = datetime.now() + timedelta(days=2)
        with patch('apis.citations_api.datetime') as mock_datetime:
            mock_datetime.now.return_value = future_time
            mock_datetime.fromtimestamp.return_value = datetime.now()  # Keep file modified time as is
            
            # Verify it's now invalid
            self.assertFalse(self.cache_manager.is_valid(self.cache_key))
            
            # Attempt to get the data (should return None)
            self.assertIsNone(self.cache_manager.get(self.cache_key))
            
    def test_dependency_tracking(self):
        """Test that dependencies between repos and caches are properly tracked."""
        # Set data with dependency on repo
        self.cache_manager.set(self.cache_key, self.test_data, self.repo_url)
        
        # Verify the dependency map was created
        repo_id = self.cache_manager._hash_key(self.repo_url)
        self.assertIn(repo_id, self.cache_manager._dependency_map)
        self.assertIn(self.cache_key, self.cache_manager._dependency_map[repo_id])
        
        # Check the dependency file was created
        deps_file = self.cache_manager.deps_dir / f"{repo_id}.json"
        self.assertTrue(deps_file.exists())
        
        # Verify content of dependency file
        with open(deps_file, 'r') as f:
            deps_data = json.load(f)
            self.assertIn(self.cache_key, deps_data)
            
    def test_invalidate_repo_caches(self):
        """Test invalidating all caches for a specific repository."""
        # Create multiple cache entries for the same repo
        keys = ["cache1", "cache2", "cache3"]
        for key in keys:
            self.cache_manager.set(key, {"data": key}, self.repo_url)
            
        # Create a cache entry for a different repo
        other_repo = "https://github.com/other/repo"
        self.cache_manager.set("other_cache", {"data": "other"}, other_repo)
        
        # Invalidate caches for the first repo
        count = self.cache_manager.invalidate_repo_caches(self.repo_url)
        
        # Verify the correct number of caches were invalidated
        self.assertEqual(count, 3)
        
        # Verify the caches for the first repo are gone
        for key in keys:
            self.assertFalse(self.cache_manager.is_valid(key))
        
        # Verify the cache for the other repo still exists
        self.assertTrue(self.cache_manager.is_valid("other_cache"))
        
    def test_clear_all(self):
        """Test clearing all caches and dependency maps."""
        # Create multiple cache entries
        self.cache_manager.set("key1", {"data": 1}, "repo1")
        self.cache_manager.set("key2", {"data": 2}, "repo2")
        self.cache_manager.set("key3", {"data": 3}, "repo1")
        
        # Clear all caches
        count = self.cache_manager.clear_all()
        
        # Verify the correct number of caches were cleared
        self.assertEqual(count, 3)
        
        # Verify dependency maps are cleared
        self.assertEqual(self.cache_manager._dependency_map, {})
        
        # Verify all caches are invalid
        self.assertFalse(self.cache_manager.is_valid("key1"))
        self.assertFalse(self.cache_manager.is_valid("key2"))
        self.assertFalse(self.cache_manager.is_valid("key3"))
        
    def test_invalid_cache_data(self):
        """Test handling of invalid cache data."""
        # Create a cache file with invalid JSON
        invalid_path = self.cache_manager._get_cache_path("invalid_json")
        with open(invalid_path, 'w') as f:
            f.write("{ invalid json")
            
        # Attempt to get the data
        result = self.cache_manager.get("invalid_json")
        self.assertIsNone(result)
        
        # Create a cache file with valid JSON but invalid structure
        wrong_structure_path = self.cache_manager._get_cache_path("wrong_structure")
        with open(wrong_structure_path, 'w') as f:
            json.dump({"not_cache_date": "something"}, f)
            
        # Verify it's considered invalid
        self.assertFalse(self.cache_manager.is_valid("wrong_structure"))

if __name__ == "__main__":
    unittest.main()