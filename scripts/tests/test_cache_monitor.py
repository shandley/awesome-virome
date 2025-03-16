#!/usr/bin/env python3
"""
Unit tests for the CacheMonitor class in monitor_cache.py
"""
import os
import sys
import json
import unittest
import tempfile
import datetime
import shutil
import importlib.util
from pathlib import Path
from unittest.mock import patch, MagicMock

# Add parent directory to path
script_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(script_dir)

# Import the necessary modules
sys.path.insert(0, script_dir)
from apis.citations_api import CacheManager

# Import the module to test
try:
    # Import with dynamic path resolution
    monitor_path = os.path.join(script_dir, "monitor_cache.py")
    spec = importlib.util.spec_from_file_location("monitor_cache", monitor_path)
    monitor_cache = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(monitor_cache)
    CacheMonitor = monitor_cache.CacheMonitor
except ImportError:
    # Fallback to direct import
    try:
        from monitor_cache import CacheMonitor
    except ImportError:
        print("Error: Could not import CacheMonitor. Make sure monitor_cache.py exists.")
        sys.exit(1)

class TestCacheMonitor(unittest.TestCase):
    """Test cases for the CacheMonitor class."""
    
    def setUp(self):
        """Set up a temporary directory for testing."""
        self.temp_dir = tempfile.mkdtemp()
        self.cache_dir = os.path.join(self.temp_dir, "cache")
        self.monitoring_dir = os.path.join(self.cache_dir, "_monitoring")
        os.makedirs(self.monitoring_dir, exist_ok=True)
        
        # Create a mock cache manager
        self.cache_manager = MagicMock(spec=CacheManager)
        
        # Set up mock metrics
        self.mock_metrics = {
            'hits': 100,
            'misses': 50,
            'invalidations': 10,
            'sets': 120,
            'hit_rate': 0.67,
            'cache_efficiency': 0.75,
            'cache_files': 80,
            'dependency_maps': 15
        }
        
        self.cache_manager.get_metrics.return_value = self.mock_metrics
        
        # Create the monitor
        self.monitor = CacheMonitor(
            self.cache_manager,
            interval=1,  # Short interval for testing
            history_dir=Path(self.monitoring_dir)
        )
        
    def tearDown(self):
        """Clean up temporary directory."""
        shutil.rmtree(self.temp_dir)
        
    def test_snapshot_metrics(self):
        """Test that snapshot_metrics correctly creates a metrics snapshot."""
        # Take a snapshot
        snapshot = self.monitor.snapshot_metrics()
        
        # Verify snapshot structure
        self.assertIn('timestamp', snapshot)
        self.assertIn('date', snapshot)
        self.assertIn('hour', snapshot)
        self.assertIn('metrics', snapshot)
        
        # Verify metrics data
        metrics = snapshot['metrics']
        self.assertEqual(metrics['hits'], 100)
        self.assertEqual(metrics['misses'], 50)
        self.assertEqual(metrics['hit_rate'], 0.67)
        self.assertEqual(metrics['cache_efficiency'], 0.75)
        
        # Verify time series was updated
        self.assertEqual(len(self.monitor.time_series['timestamps']), 1)
        self.assertEqual(len(self.monitor.time_series['hit_rate']), 1)
        self.assertEqual(self.monitor.time_series['hit_rate'][0], 67.0)  # 0.67 * 100
        
        # Verify history was saved
        history_file = os.path.join(self.monitoring_dir, "metrics_history.json")
        self.assertTrue(os.path.exists(history_file))
        
        # Verify history content
        with open(history_file, 'r') as f:
            history = json.load(f)
            self.assertIn('daily_metrics', history)
            self.assertIn('hourly_metrics', history)
            self.assertIn('peak_values', history)
            self.assertIn('baseline', history)
            
            # Check peak values
            self.assertEqual(history['peak_values']['hit_rate'], 0.67)
            self.assertEqual(history['peak_values']['efficiency'], 0.75)
            
    def test_analyze_snapshot(self):
        """Test that analyze_snapshot correctly analyzes a metrics snapshot."""
        # Create a snapshot with known metrics
        snapshot = {
            'timestamp': datetime.datetime.now().isoformat(),
            'metrics': self.mock_metrics,
            'derived': {
                'hit_rate_change': 0.05,
                'invalidation_rate': 2,
                'growth_rate': 10
            }
        }
        
        # Analyze the snapshot
        analysis = self.monitor.analyze_snapshot(snapshot)
        
        # Verify analysis structure
        self.assertIn('health', analysis)
        self.assertIn('warnings', analysis)
        self.assertIn('recommendations', analysis)
        self.assertIn('insights', analysis)
        
        # Verify analysis content
        self.assertEqual(analysis['health'], 'good')
        self.assertTrue(any("hit rate" in insight.lower() for insight in analysis['insights']))
        self.assertTrue(any("efficiency" in insight.lower() for insight in analysis['insights']))
        
        # Test with poor metrics
        snapshot['metrics']['hit_rate'] = 0.25
        snapshot['metrics']['cache_efficiency'] = -0.1
        snapshot['derived']['invalidation_rate'] = 15
        
        analysis = self.monitor.analyze_snapshot(snapshot)
        
        self.assertEqual(analysis['health'], 'poor')
        self.assertTrue(len(analysis['warnings']) >= 2)
        self.assertTrue(len(analysis['recommendations']) >= 2)
        
    def test_generate_report(self):
        """Test that generate_report creates a readable report."""
        # Create a snapshot and analysis
        snapshot = {
            'timestamp': datetime.datetime.now().isoformat(),
            'metrics': self.mock_metrics
        }
        
        analysis = {
            'health': 'good',
            'warnings': [],
            'recommendations': [],
            'insights': [
                'Hit rate is healthy at 67.0%',
                'Cache efficiency is healthy at 75.0%'
            ]
        }
        
        # Generate the report
        report = self.monitor.generate_report(snapshot, analysis)
        
        # Verify report content
        self.assertIn("CACHE MONITORING REPORT", report)
        self.assertIn("Status: GOOD", report)
        self.assertIn("Hit rate: 67.0%", report)
        self.assertIn("Efficiency: 75.0%", report)
        self.assertIn("💡 Hit rate is healthy", report)
        
    def test_export_metrics_to_csv(self):
        """Test exporting metrics to CSV."""
        # Set up daily metrics in history
        today = datetime.datetime.now().strftime('%Y-%m-%d')
        yesterday = (datetime.datetime.now() - datetime.timedelta(days=1)).strftime('%Y-%m-%d')
        
        self.monitor.history['daily_metrics'] = {
            today: {
                'snapshots': 24,
                'avg_hit_rate': 0.7,
                'avg_efficiency': 0.8,
                'total_invalidations': 12,
                'max_cache_size': 85
            },
            yesterday: {
                'snapshots': 24,
                'avg_hit_rate': 0.65,
                'avg_efficiency': 0.75,
                'total_invalidations': 10,
                'max_cache_size': 80
            }
        }
        
        # Export to CSV
        csv_path = os.path.join(self.monitoring_dir, "test_export.csv")
        result_path = self.monitor.export_metrics_to_csv(csv_path)
        
        # Verify CSV was created
        self.assertEqual(result_path, csv_path)
        self.assertTrue(os.path.exists(csv_path))
        
        # Verify CSV content
        with open(csv_path, 'r') as f:
            lines = f.readlines()
            self.assertEqual(len(lines), 3)  # Header + 2 days
            self.assertIn("Date,Snapshots,Avg Hit Rate,Avg Efficiency,Total Invalidations,Max Cache Size", lines[0])
    
if __name__ == "__main__":
    import importlib.util
    unittest.main()