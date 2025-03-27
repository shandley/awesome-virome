#!/usr/bin/env python3
"""
Script to store test metrics history in a dedicated branch.
This script takes results from various GitHub Actions workflow steps,
stores them in a structured format, and commits them to a metrics-history branch.
"""

import os
import sys
import json
import datetime
import subprocess
from pathlib import Path

# Configuration
METRICS_DIR = "metrics_history"
PERFORMANCE_DIR = os.path.join(METRICS_DIR, "performance")
VALIDATION_DIR = os.path.join(METRICS_DIR, "validation")
LINKS_DIR = os.path.join(METRICS_DIR, "links")
SUMMARY_FILE = os.path.join(METRICS_DIR, "metrics_summary.json")

def setup_git_config():
    """Configure git for the GitHub Actions environment."""
    subprocess.run(["git", "config", "user.name", "github-actions"])
    subprocess.run(["git", "config", "user.email", "github-actions@github.com"])

def ensure_directories_exist():
    """Create necessary directories if they don't exist."""
    for directory in [METRICS_DIR, PERFORMANCE_DIR, VALIDATION_DIR, LINKS_DIR]:
        Path(directory).mkdir(parents=True, exist_ok=True)

def load_or_create_summary():
    """Load the metrics summary file or create it if it doesn't exist."""
    if os.path.exists(SUMMARY_FILE):
        with open(SUMMARY_FILE, 'r') as f:
            try:
                return json.load(f)
            except json.JSONDecodeError:
                print(f"Warning: {SUMMARY_FILE} was not valid JSON. Creating new summary.")
    
    # Initialize new summary structure
    return {
        "last_updated": "",
        "runs": [],
        "performance_trends": {
            "data_load_time": [],
            "readme_parse_time": []
        },
        "validation_stats": {
            "total_runs": 0,
            "successful_runs": 0
        },
        "link_health": {
            "history": []
        }
    }

def update_performance_metrics(summary, run_id, benchmark_results):
    """Update performance metrics with the latest benchmark results."""
    # Parse the benchmark results string
    # Example format: "Data load: 0.1234s, README parse: 0.5678s"
    try:
        data_load = float(benchmark_results.split("Data load: ")[1].split("s,")[0])
        readme_parse = float(benchmark_results.split("README parse: ")[1].split("s")[0])
    except (IndexError, ValueError) as e:
        print(f"Error parsing benchmark results: {e}")
        return
    
    # Create the timestamped metrics entry
    timestamp = datetime.datetime.now().isoformat()
    
    # Add to the summary trends
    summary["performance_trends"]["data_load_time"].append({
        "timestamp": timestamp,
        "run_id": run_id,
        "value": data_load
    })
    
    summary["performance_trends"]["readme_parse_time"].append({
        "timestamp": timestamp,
        "run_id": run_id,
        "value": readme_parse
    })
    
    # Save a detailed performance record
    performance_file = os.path.join(PERFORMANCE_DIR, f"{run_id}.json")
    with open(performance_file, 'w') as f:
        json.dump({
            "timestamp": timestamp,
            "run_id": run_id,
            "data_load_time": data_load,
            "readme_parse_time": readme_parse,
            "raw_results": benchmark_results
        }, f, indent=2)
    
    return {
        "data_load_time": data_load,
        "readme_parse_time": readme_parse
    }

def update_link_metrics(summary, run_id, link_report_path):
    """Update link health metrics from the Lychee link checker report."""
    if not os.path.exists(link_report_path):
        print(f"Warning: Link report file {link_report_path} not found")
        return None
    
    # Parse the Lychee markdown output to count issues
    with open(link_report_path, 'r') as f:
        content = f.read()
    
    # Simple parsing of the markdown report
    broken_count = content.count("✖")
    excluded_count = content.count("➖")
    success_count = content.count("✓")
    
    timestamp = datetime.datetime.now().isoformat()
    
    # Add entry to link health history
    link_entry = {
        "timestamp": timestamp,
        "run_id": run_id,
        "broken_links": broken_count,
        "successful_links": success_count,
        "excluded_links": excluded_count,
        "total_links": broken_count + success_count + excluded_count
    }
    
    summary["link_health"]["history"].append(link_entry)
    
    # Save detailed link report
    links_file = os.path.join(LINKS_DIR, f"{run_id}.json")
    with open(links_file, 'w') as f:
        json.dump(link_entry, f, indent=2)
    
    return link_entry

def update_validation_metrics(summary, run_id, workflow_status):
    """Update validation metrics based on workflow job statuses."""
    timestamp = datetime.datetime.now().isoformat()
    
    validation_status = {
        "timestamp": timestamp,
        "run_id": run_id,
        "status": workflow_status,
        "jobs": {}
    }
    
    # Update the validation stats
    summary["validation_stats"]["total_runs"] += 1
    if workflow_status == "success":
        summary["validation_stats"]["successful_runs"] += 1
    
    # Save detailed validation report
    validation_file = os.path.join(VALIDATION_DIR, f"{run_id}.json")
    with open(validation_file, 'w') as f:
        json.dump(validation_status, f, indent=2)
    
    return validation_status

def update_summary(summary, run_id, perf_metrics, link_metrics, validation_metrics):
    """Update the main summary file with metrics from this run."""
    timestamp = datetime.datetime.now().isoformat()
    summary["last_updated"] = timestamp
    
    # Add this run to the runs list
    run_entry = {
        "run_id": run_id,
        "timestamp": timestamp,
        "performance": perf_metrics,
        "links": link_metrics,
        "validation": validation_metrics
    }
    
    # Add to the beginning of the runs list
    summary["runs"].insert(0, run_entry)
    
    # Keep only the last 100 runs to avoid the file growing too large
    if len(summary["runs"]) > 100:
        summary["runs"] = summary["runs"][:100]
    
    # Keep only the last 100 entries in each trend
    for trend in summary["performance_trends"]:
        if len(summary["performance_trends"][trend]) > 100:
            summary["performance_trends"][trend] = summary["performance_trends"][trend][-100:]
    
    # Keep only the last 100 link health entries
    if len(summary["link_health"]["history"]) > 100:
        summary["link_health"]["history"] = summary["link_health"]["history"][-100:]
    
    # Save the updated summary
    with open(SUMMARY_FILE, 'w') as f:
        json.dump(summary, f, indent=2)

def commit_and_push():
    """Commit the updated metrics files and push to the metrics-history branch."""
    subprocess.run(["git", "add", METRICS_DIR])
    commit_message = f"Update metrics history - {datetime.datetime.now().isoformat()}"
    subprocess.run(["git", "commit", "-m", commit_message])
    subprocess.run(["git", "push", "origin", "metrics-history"])

def main():
    # Get inputs from environment
    run_id = os.environ.get("GITHUB_RUN_ID", f"local-{int(datetime.datetime.now().timestamp())}")
    benchmark_results = os.environ.get("BENCHMARK_RESULTS", "Data load: 0.0s, README parse: 0.0s")
    workflow_status = os.environ.get("WORKFLOW_STATUS", "unknown")
    link_report_path = os.environ.get("LINK_REPORT_PATH", "./lychee/out.md")
    
    # Setup git config
    setup_git_config()
    
    # Ensure directories exist
    ensure_directories_exist()
    
    # Load or create summary
    summary = load_or_create_summary()
    
    # Update metrics
    perf_metrics = update_performance_metrics(summary, run_id, benchmark_results)
    link_metrics = update_link_metrics(summary, run_id, link_report_path)
    validation_metrics = update_validation_metrics(summary, run_id, workflow_status)
    
    # Update the summary
    update_summary(summary, run_id, perf_metrics, link_metrics, validation_metrics)
    
    # Commit and push changes
    commit_and_push()
    
    print(f"Metrics history updated for run {run_id}")
    return 0

if __name__ == "__main__":
    sys.exit(main())