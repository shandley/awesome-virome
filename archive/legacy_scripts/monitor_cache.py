#!/usr/bin/env python3
"""
Advanced cache monitoring system for the awesome-virome repository.

This script provides comprehensive monitoring for the cache system, including:
- Real-time performance metrics collection
- Trend analysis across multiple runs
- API rate limit monitoring
- Cache efficiency recommendations
- Export of metrics to various formats (JSON, CSV)
- System health checks
"""

import os
import sys
import json
import time
import argparse
import logging
import datetime
import csv
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Any, Optional, Tuple

# Optional imports for visualization - make these not required for basic operation
MATPLOTLIB_AVAILABLE = False
NUMPY_AVAILABLE = False
try:
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.ticker import MaxNLocator
    MATPLOTLIB_AVAILABLE = True
    NUMPY_AVAILABLE = True
except ImportError:
    pass

# Import the cache manager
try:
    from apis.citations_api import cache_manager, RateLimiter
    HAS_CACHE_MANAGER = True
except ImportError:
    HAS_CACHE_MANAGER = False
    print("Error: Cache manager not found. Make sure apis/citations_api.py exists and is properly configured.")
    sys.exit(1)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("cache_monitoring.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Constants
DEFAULT_MONITORING_INTERVAL = 60  # seconds
HISTORY_DIR = Path(os.path.join("metadata", "cache", "_monitoring"))
HISTORY_DIR.mkdir(exist_ok=True, parents=True)

class CacheMonitor:
    """Advanced monitoring system for the cache manager."""
    
    def __init__(self, 
                cache_manager, 
                interval: int = DEFAULT_MONITORING_INTERVAL,
                history_dir: Path = HISTORY_DIR):
        """
        Initialize the cache monitor.
        
        Args:
            cache_manager: The cache manager to monitor
            interval: Monitoring interval in seconds
            history_dir: Directory to store historical metrics
        """
        self.cache_manager = cache_manager
        self.interval = interval
        self.history_dir = history_dir
        self.history_dir.mkdir(exist_ok=True)
        
        # Load historical metrics
        self.history = self._load_history()
        
        # Time series data
        self.time_series = {
            'timestamps': [],
            'hit_rate': [],
            'miss_rate': [],
            'efficiency': [],
            'invalidations': [],
            'cache_size': []
        }
        
    def _load_history(self) -> Dict[str, Any]:
        """Load historical metrics from storage."""
        history_file = self.history_dir / "metrics_history.json"
        if history_file.exists():
            try:
                with open(history_file, 'r') as f:
                    return json.load(f)
            except (json.JSONDecodeError, IOError) as e:
                logger.warning(f"Failed to load metrics history: {e}")
        
        # Default history structure
        return {
            'daily_metrics': {},
            'hourly_metrics': {},
            'peak_values': {
                'hit_rate': 0,
                'miss_rate': 0,
                'efficiency': 0,
                'cache_size': 0,
                'invalidations': 0
            },
            'baseline': None,
            'last_updated': None
        }
    
    def _save_history(self):
        """Save historical metrics to storage."""
        history_file = self.history_dir / "metrics_history.json"
        try:
            # Update the last updated timestamp
            self.history['last_updated'] = datetime.datetime.now().isoformat()
            
            with open(history_file, 'w') as f:
                json.dump(self.history, f, indent=2)
        except IOError as e:
            logger.warning(f"Failed to save metrics history: {e}")
    
    def snapshot_metrics(self) -> Dict[str, Any]:
        """Take a snapshot of current cache metrics."""
        metrics = self.cache_manager.get_metrics()
        
        # Calculate rate of change for key metrics
        now = datetime.datetime.now()
        timestamp = now.isoformat()
        today = now.strftime('%Y-%m-%d')
        hour = now.strftime('%Y-%m-%d %H:00')
        
        # Get previous metrics for comparison
        prev_metrics = self.history.get('baseline')
        
        # Create snapshot with derived metrics
        snapshot = {
            'timestamp': timestamp,
            'date': today,
            'hour': hour,
            'metrics': metrics,
            'derived': {}
        }
        
        # Add time series data
        self.time_series['timestamps'].append(timestamp)
        self.time_series['hit_rate'].append(metrics.get('hit_rate', 0) * 100)
        self.time_series['miss_rate'].append(100 - (metrics.get('hit_rate', 0) * 100))
        self.time_series['efficiency'].append(metrics.get('cache_efficiency', 0) * 100)
        self.time_series['invalidations'].append(metrics.get('invalidations', 0))
        self.time_series['cache_size'].append(metrics.get('cache_files', 0))
        
        # Calculate rates of change if we have previous metrics
        if prev_metrics:
            time_diff = (datetime.datetime.fromisoformat(timestamp) - 
                         datetime.datetime.fromisoformat(prev_metrics['timestamp'])).total_seconds()
            
            if time_diff > 0:
                # Calculate hourly rates of change
                hours = time_diff / 3600
                snapshot['derived']['hit_rate_change'] = (
                    (metrics.get('hit_rate', 0) - prev_metrics['metrics'].get('hit_rate', 0)) / hours
                )
                snapshot['derived']['invalidation_rate'] = (
                    (metrics.get('invalidations', 0) - prev_metrics['metrics'].get('invalidations', 0)) / hours
                )
                snapshot['derived']['growth_rate'] = (
                    (metrics.get('cache_files', 0) - prev_metrics['metrics'].get('cache_files', 0)) / hours
                )
        
        # Update daily and hourly aggregates
        if today not in self.history['daily_metrics']:
            self.history['daily_metrics'][today] = {
                'snapshots': 0,
                'avg_hit_rate': 0,
                'avg_efficiency': 0,
                'total_invalidations': 0,
                'max_cache_size': 0
            }
        
        daily = self.history['daily_metrics'][today]
        daily['snapshots'] += 1
        daily['avg_hit_rate'] = (
            (daily['avg_hit_rate'] * (daily['snapshots'] - 1) + metrics.get('hit_rate', 0)) / 
            daily['snapshots']
        )
        daily['avg_efficiency'] = (
            (daily['avg_efficiency'] * (daily['snapshots'] - 1) + metrics.get('cache_efficiency', 0)) /
            daily['snapshots']
        )
        daily['total_invalidations'] = max(
            daily['total_invalidations'],
            metrics.get('invalidations', 0)
        )
        daily['max_cache_size'] = max(
            daily['max_cache_size'],
            metrics.get('cache_files', 0)
        )
        
        # Update hourly metrics
        if hour not in self.history['hourly_metrics']:
            self.history['hourly_metrics'][hour] = {
                'snapshots': 0,
                'avg_hit_rate': 0,
                'avg_efficiency': 0,
                'total_invalidations': 0,
                'max_cache_size': 0
            }
        
        hourly = self.history['hourly_metrics'][hour]
        hourly['snapshots'] += 1
        hourly['avg_hit_rate'] = (
            (hourly['avg_hit_rate'] * (hourly['snapshots'] - 1) + metrics.get('hit_rate', 0)) / 
            hourly['snapshots']
        )
        hourly['avg_efficiency'] = (
            (hourly['avg_efficiency'] * (hourly['snapshots'] - 1) + metrics.get('cache_efficiency', 0)) /
            hourly['snapshots']
        )
        hourly['total_invalidations'] = max(
            hourly['total_invalidations'],
            metrics.get('invalidations', 0)
        )
        hourly['max_cache_size'] = max(
            hourly['max_cache_size'],
            metrics.get('cache_files', 0)
        )
        
        # Update peak values
        self.history['peak_values']['hit_rate'] = max(
            self.history['peak_values']['hit_rate'],
            metrics.get('hit_rate', 0)
        )
        self.history['peak_values']['miss_rate'] = max(
            self.history['peak_values']['miss_rate'],
            1 - metrics.get('hit_rate', 0) if metrics.get('hit_rate') is not None else 0
        )
        self.history['peak_values']['efficiency'] = max(
            self.history['peak_values']['efficiency'],
            metrics.get('cache_efficiency', 0)
        )
        self.history['peak_values']['cache_size'] = max(
            self.history['peak_values']['cache_size'],
            metrics.get('cache_files', 0)
        )
        self.history['peak_values']['invalidations'] = max(
            self.history['peak_values']['invalidations'],
            metrics.get('invalidations', 0)
        )
        
        # Set as baseline if not set
        if not self.history['baseline']:
            self.history['baseline'] = snapshot
        
        # Save updated history
        self._save_history()
        
        return snapshot
    
    def analyze_snapshot(self, snapshot: Dict[str, Any]) -> Dict[str, Any]:
        """
        Analyze a metrics snapshot and generate insights.
        
        Args:
            snapshot: A metrics snapshot
            
        Returns:
            Dict of analysis results and recommendations
        """
        metrics = snapshot['metrics']
        derived = snapshot.get('derived', {})
        
        analysis = {
            'health': 'good',
            'warnings': [],
            'recommendations': [],
            'insights': []
        }
        
        # Check hit rate (below 50% is concerning)
        hit_rate = metrics.get('hit_rate', 0) * 100
        if hit_rate < 30:
            analysis['health'] = 'poor'
            analysis['warnings'].append(f"Hit rate is very low ({hit_rate:.1f}%)")
            analysis['recommendations'].append(
                "Consider examining cache expiration settings and dependency maps"
            )
        elif hit_rate < 50:
            analysis['health'] = 'fair'
            analysis['warnings'].append(f"Hit rate is below average ({hit_rate:.1f}%)")
            analysis['recommendations'].append(
                "Consider selective adjustment of cache expiration times"
            )
        else:
            analysis['insights'].append(f"Hit rate is healthy at {hit_rate:.1f}%")
        
        # Check cache efficiency (should be positive)
        efficiency = metrics.get('cache_efficiency', 0) * 100
        if efficiency < 0:
            analysis['health'] = 'poor'
            analysis['warnings'].append(f"Cache efficiency is negative ({efficiency:.1f}%)")
            analysis['recommendations'].append(
                "Too many invalidations compared to cache hits. Review repository dependency tracking."
            )
        elif efficiency < 50:
            if analysis['health'] != 'poor':
                analysis['health'] = 'fair'
            analysis['warnings'].append(f"Cache efficiency is low ({efficiency:.1f}%)")
            analysis['recommendations'].append(
                "Consider reviewing repository dependency tracking for over-invalidation"
            )
        else:
            analysis['insights'].append(f"Cache efficiency is healthy at {efficiency:.1f}%")
        
        # Check invalidation rate trends
        if 'invalidation_rate' in derived and derived['invalidation_rate'] > 10:
            analysis['warnings'].append(
                f"High invalidation rate: {derived['invalidation_rate']:.1f} invalidations/hour"
            )
            analysis['recommendations'].append(
                "Consider optimizing repository dependency tracking to reduce unnecessary invalidations"
            )
        
        # Check cache growth rate
        if 'growth_rate' in derived and derived['growth_rate'] > 50:
            analysis['warnings'].append(
                f"High cache growth rate: {derived['growth_rate']:.1f} entries/hour"
            )
            analysis['recommendations'].append(
                "Monitor disk usage and consider more aggressive expiration for less important data"
            )
        
        # Check hit rate trends
        if 'hit_rate_change' in derived:
            if derived['hit_rate_change'] < -0.05:  # 5% decrease per hour
                analysis['warnings'].append(
                    f"Declining hit rate: {derived['hit_rate_change'] * 100:.1f}% change per hour"
                )
                analysis['recommendations'].append(
                    "Investigate recent changes in cache usage patterns"
                )
            elif derived['hit_rate_change'] > 0.05:  # 5% increase per hour
                analysis['insights'].append(
                    f"Improving hit rate: +{derived['hit_rate_change'] * 100:.1f}% per hour"
                )
        
        return analysis
    
    def generate_report(self, snapshot: Dict[str, Any], analysis: Dict[str, Any]) -> str:
        """
        Generate a human-readable report of cache metrics and analysis.
        
        Args:
            snapshot: A metrics snapshot
            analysis: Analysis of the snapshot
            
        Returns:
            Formatted report string
        """
        metrics = snapshot['metrics']
        timestamp = snapshot['timestamp']
        
        # Calculate values for report
        total_requests = metrics.get('hits', 0) + metrics.get('misses', 0)
        hit_rate = metrics.get('hit_rate', 0) * 100
        efficiency = metrics.get('cache_efficiency', 0) * 100
        cache_files = metrics.get('cache_files', 0)
        
        # Format the report
        report = [
            "===== CACHE MONITORING REPORT =====",
            f"Timestamp: {timestamp}",
            f"Status: {analysis['health'].upper()}",
            "",
            "----- Cache Performance -----",
            f"Hit rate: {hit_rate:.1f}% ({metrics.get('hits', 0)} hits, {metrics.get('misses', 0)} misses)",
            f"Efficiency: {efficiency:.1f}%",
            f"Total cache entries: {cache_files}",
            f"Cache sets: {metrics.get('sets', 0)}",
            f"Invalidations: {metrics.get('invalidations', 0)}",
            "",
        ]
        
        # Add historical comparison if available
        baseline = self.history.get('baseline')
        if baseline and baseline['timestamp'] != timestamp:
            baseline_metrics = baseline['metrics']
            baseline_hit_rate = baseline_metrics.get('hit_rate', 0) * 100
            baseline_efficiency = baseline_metrics.get('cache_efficiency', 0) * 100
            baseline_files = baseline_metrics.get('cache_files', 0)
            
            report.extend([
                "----- Change Since Baseline -----",
                f"Hit rate: {hit_rate - baseline_hit_rate:+.1f}% points",
                f"Efficiency: {efficiency - baseline_efficiency:+.1f}% points",
                f"Cache entries: {cache_files - baseline_files:+d}",
                f"Invalidations: +{metrics.get('invalidations', 0) - baseline_metrics.get('invalidations', 0)}",
                "",
            ])
        
        # Add warnings
        if analysis['warnings']:
            report.append("----- Warnings -----")
            for warning in analysis['warnings']:
                report.append(f"⚠️ {warning}")
            report.append("")
        
        # Add recommendations
        if analysis['recommendations']:
            report.append("----- Recommendations -----")
            for recommendation in analysis['recommendations']:
                report.append(f"📋 {recommendation}")
            report.append("")
        
        # Add insights
        if analysis['insights']:
            report.append("----- Insights -----")
            for insight in analysis['insights']:
                report.append(f"💡 {insight}")
            report.append("")
        
        report.append("===============================")
        
        return "\n".join(report)
    
    def export_metrics_to_csv(self, output_path: str = None):
        """
        Export metrics history to CSV format.
        
        Args:
            output_path: Path to save the CSV file (default: metrics_export.csv in history dir)
        """
        if output_path is None:
            output_path = self.history_dir / "metrics_export.csv"
        
        try:
            # Sort and keep most recent entries if too many (max 10000 rows)
            daily_metrics = sorted(
                [(date, metrics) for date, metrics in self.history['daily_metrics'].items()],
                key=lambda x: x[0],
                reverse=True
            )[:10000]
            
            with open(output_path, 'w', newline='') as csvfile:
                writer = csv.writer(csvfile)
                # Write header
                writer.writerow([
                    'Date', 'Snapshots', 'Avg Hit Rate', 'Avg Efficiency', 
                    'Total Invalidations', 'Max Cache Size'
                ])
                
                # Write data
                for date, metrics in daily_metrics:
                    writer.writerow([
                        date,
                        metrics['snapshots'],
                        f"{metrics['avg_hit_rate'] * 100:.2f}%",
                        f"{metrics['avg_efficiency'] * 100:.2f}%",
                        metrics['total_invalidations'],
                        metrics['max_cache_size']
                    ])
            
            logger.info(f"Exported metrics history to {output_path}")
            return output_path
        except Exception as e:
            logger.error(f"Failed to export metrics to CSV: {e}")
            return None
    
    def generate_performance_graphs(self, output_dir: str = None):
        """
        Generate performance graphs from time series data.
        
        Args:
            output_dir: Directory to save the graphs (default: history_dir / "graphs")
        """
        if not MATPLOTLIB_AVAILABLE or not NUMPY_AVAILABLE:
            logger.warning("Matplotlib and/or NumPy not available. Skipping graph generation.")
            return None
            
        if not output_dir:
            output_dir = self.history_dir / "graphs"
        
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True, parents=True)
        
        try:
            # Convert timestamps to readable format
            timestamps = [datetime.datetime.fromisoformat(ts) for ts in self.time_series['timestamps']]
            
            if not timestamps:
                logger.warning("No time series data available for graphing")
                return None
            
            # Set up the figures
            plt.style.use('ggplot')
            
            # 1. Hit Rate and Miss Rate
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.plot(timestamps, self.time_series['hit_rate'], 'g-', label='Hit Rate (%)')
            ax.plot(timestamps, self.time_series['miss_rate'], 'r-', label='Miss Rate (%)')
            ax.set_title('Cache Hit Rate vs Miss Rate')
            ax.set_xlabel('Time')
            ax.set_ylabel('Percentage (%)')
            ax.grid(True)
            ax.legend()
            fig.autofmt_xdate()
            fig.tight_layout()
            hit_rate_path = output_path / "hit_rate.png"
            fig.savefig(hit_rate_path)
            plt.close(fig)
            
            # 2. Cache Efficiency
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.plot(timestamps, self.time_series['efficiency'], 'b-', label='Efficiency (%)')
            ax.set_title('Cache Efficiency Over Time')
            ax.set_xlabel('Time')
            ax.set_ylabel('Efficiency (%)')
            ax.grid(True)
            ax.legend()
            fig.autofmt_xdate()
            fig.tight_layout()
            efficiency_path = output_path / "efficiency.png"
            fig.savefig(efficiency_path)
            plt.close(fig)
            
            # 3. Invalidations
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.plot(timestamps, self.time_series['invalidations'], 'r-', label='Invalidations')
            ax.set_title('Cache Invalidations Over Time')
            ax.set_xlabel('Time')
            ax.set_ylabel('Count')
            ax.grid(True)
            ax.legend()
            fig.autofmt_xdate()
            fig.tight_layout()
            invalidations_path = output_path / "invalidations.png"
            fig.savefig(invalidations_path)
            plt.close(fig)
            
            # 4. Cache Size
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.plot(timestamps, self.time_series['cache_size'], 'b-', label='Cache Size')
            ax.set_title('Cache Size Over Time')
            ax.set_xlabel('Time')
            ax.set_ylabel('Number of Cache Files')
            ax.grid(True)
            ax.legend()
            ax.yaxis.set_major_locator(MaxNLocator(integer=True))
            fig.autofmt_xdate()
            fig.tight_layout()
            cache_size_path = output_path / "cache_size.png"
            fig.savefig(cache_size_path)
            plt.close(fig)
            
            logger.info(f"Generated performance graphs in {output_path}")
            
            # Return paths to the generated graphs
            return {
                'hit_rate': str(hit_rate_path),
                'efficiency': str(efficiency_path),
                'invalidations': str(invalidations_path),
                'cache_size': str(cache_size_path)
            }
        except Exception as e:
            logger.error(f"Failed to generate performance graphs: {e}")
            return None
    
    def monitor_once(self, report=True, graphs=False):
        """
        Run a single monitoring cycle.
        
        Args:
            report: Whether to print the report
            graphs: Whether to generate graphs
            
        Returns:
            Tuple of (snapshot, analysis)
        """
        snapshot = self.snapshot_metrics()
        analysis = self.analyze_snapshot(snapshot)
        
        if report:
            report_text = self.generate_report(snapshot, analysis)
            print(report_text)
        
        if graphs and MATPLOTLIB_AVAILABLE and NUMPY_AVAILABLE:
            self.generate_performance_graphs()
        
        return snapshot, analysis
    
    def monitor_continuously(self, duration=None, report=True, graphs=False):
        """
        Run continuous monitoring for a specified duration.
        
        Args:
            duration: Duration in seconds (None for indefinite monitoring)
            report: Whether to print reports
            graphs: Whether to generate graphs periodically
            
        Returns:
            List of snapshots collected during monitoring
        """
        snapshots = []
        start_time = time.time()
        graphs_interval = 10  # Generate graphs every 10 monitoring cycles
        cycle_count = 0
        
        try:
            logger.info(f"Starting continuous monitoring with {self.interval}s interval")
            
            while True:
                # Check if we've reached the specified duration
                if duration and (time.time() - start_time) >= duration:
                    break
                
                # Take a snapshot and analyze
                snapshot, analysis = self.monitor_once(report=report, graphs=False)
                snapshots.append(snapshot)
                
                # Generate graphs periodically
                cycle_count += 1
                if graphs and MATPLOTLIB_AVAILABLE and NUMPY_AVAILABLE and cycle_count % graphs_interval == 0:
                    self.generate_performance_graphs()
                
                # Wait for the next monitoring cycle
                time.sleep(self.interval)
                
        except KeyboardInterrupt:
            logger.info("Monitoring stopped by user")
        
        # Generate final graphs
        if graphs and MATPLOTLIB_AVAILABLE and NUMPY_AVAILABLE:
            self.generate_performance_graphs()
        
        logger.info(f"Monitoring completed. Collected {len(snapshots)} snapshots.")
        return snapshots

def main():
    """Main entry point for the script."""
    parser = argparse.ArgumentParser(description="Advanced Cache Monitoring Tool")
    
    # Add command line arguments
    parser.add_argument('--continuous', action='store_true', help='Run continuous monitoring')
    parser.add_argument('--interval', type=int, default=DEFAULT_MONITORING_INTERVAL, 
                        help=f'Monitoring interval in seconds (default: {DEFAULT_MONITORING_INTERVAL})')
    parser.add_argument('--duration', type=int, help='Duration of continuous monitoring in seconds')
    parser.add_argument('--graphs', action='store_true', help='Generate performance graphs')
    parser.add_argument('--export-csv', action='store_true', help='Export metrics to CSV')
    parser.add_argument('--csv-path', type=str, help='Path for CSV export')
    parser.add_argument('--silent', action='store_true', help='Run without printing reports')
    parser.add_argument('--summary', action='store_true', help='Print a short summary of the cache state')
    
    args = parser.parse_args()
    
    # No arguments provided, show help
    if len(sys.argv) == 1:
        parser.print_help()
        return
    
    # Create the cache monitor
    monitor = CacheMonitor(cache_manager, interval=args.interval)
    
    # Handle summary request (simple output specifically for CI workflows)
    if args.summary:
        snapshot, analysis = monitor.monitor_once(report=False, graphs=False)
        metrics = snapshot['metrics']
        print("===== CACHE SUMMARY =====")
        print(f"Status: {analysis['health']}")
        print(f"Cache entries: {metrics.get('cache_files', 0)}")
        print(f"Hit rate: {metrics.get('hit_rate', 0) * 100:.1f}%")
        print(f"Efficiency: {metrics.get('cache_efficiency', 0) * 100:.1f}%")
        print(f"Warnings: {len(analysis['warnings'])}")
        print("========================")
        return
    
    # Run monitoring based on args
    if args.continuous:
        monitor.monitor_continuously(
            duration=args.duration,
            report=not args.silent,
            graphs=args.graphs
        )
    else:
        monitor.monitor_once(report=not args.silent, graphs=args.graphs)
    
    # Export metrics to CSV if requested
    if args.export_csv:
        monitor.export_metrics_to_csv(args.csv_path)

if __name__ == "__main__":
    main()