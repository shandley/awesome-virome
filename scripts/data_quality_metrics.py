#!/usr/bin/env python3
"""
Data Quality Metrics Collection

This script analyzes the data quality of the awesome-virome repository,
tracking metrics over time to ensure data consistency and completeness.
"""

import os
import sys
import json
import logging
import argparse
import datetime
from pathlib import Path
from collections import defaultdict, Counter
from typing import Dict, List, Any, Optional, Tuple

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("data_quality.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Constants
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(SCRIPT_DIR)
DATA_JSON_PATH = os.path.join(ROOT_DIR, "data.json")
METADATA_DIR = os.path.join(ROOT_DIR, "metadata")
REPORTS_DIR = os.path.join(ROOT_DIR, "reports")
METRICS_HISTORY_FILE = os.path.join(REPORTS_DIR, "data_quality_history.json")

# Define expected fields and their importance
CRITICAL_FIELDS = ["name", "description", "repository", "id"]
IMPORTANT_FIELDS = ["homepage", "language", "license", "stars", "last_updated"]
METADATA_FIELDS = [
    "stars", "forks", "open_issues", "license", "created_at", "updated_at",
    "is_archived", "default_branch", "last_commit", "contributors_count",
    "citation", "doi", "publication"
]
BIOINFORMATICS_FIELDS = [
    "package_managers", "container_images", "installation_methods",
    "input_formats", "output_formats", "dependencies", "language_versions"
]
ACADEMIC_FIELDS = [
    "citation_count", "first_published", "has_preprint", "journal",
    "citation_trend"
]

def load_data_json() -> Dict[str, Any]:
    """Load the current data.json file"""
    try:
        with open(DATA_JSON_PATH, 'r') as f:
            return json.load(f)
    except (json.JSONDecodeError, FileNotFoundError) as e:
        logger.error(f"Error loading data.json: {e}")
        return {"tools": []}

def load_metadata_for_tool(tool_id: str) -> Dict[str, Any]:
    """Load the metadata file for a specific tool"""
    metadata_path = os.path.join(METADATA_DIR, f"{tool_id}.json")
    try:
        if os.path.exists(metadata_path):
            with open(metadata_path, 'r') as f:
                return json.load(f)
        return {}
    except (json.JSONDecodeError, FileNotFoundError) as e:
        logger.warning(f"Error loading metadata for {tool_id}: {e}")
        return {}

def load_metrics_history() -> Dict[str, Any]:
    """Load historical metrics data"""
    os.makedirs(REPORTS_DIR, exist_ok=True)
    
    if os.path.exists(METRICS_HISTORY_FILE):
        try:
            with open(METRICS_HISTORY_FILE, 'r') as f:
                return json.load(f)
        except (json.JSONDecodeError, FileNotFoundError) as e:
            logger.warning(f"Error loading metrics history: {e}")
    
    # Default structure if file doesn't exist
    return {
        "history": [],
        "first_recorded": datetime.datetime.now().isoformat()
    }

def save_metrics_history(history: Dict[str, Any]):
    """Save historical metrics data"""
    try:
        with open(METRICS_HISTORY_FILE, 'w') as f:
            json.dump(history, f, indent=2)
    except (IOError, OSError) as e:
        logger.error(f"Error saving metrics history: {e}")

def analyze_data_quality() -> Dict[str, Any]:
    """
    Analyze the quality of data in the repository.
    
    Returns:
        Dict containing various quality metrics
    """
    data = load_data_json()
    tools = data.get("tools", [])
    
    if not tools:
        logger.error("No tools found in data.json")
        return {}
    
    metrics = {
        "timestamp": datetime.datetime.now().isoformat(),
        "total_tools": len(tools),
        "fields": defaultdict(int),
        "critical_fields_completion": 0,
        "important_fields_completion": 0,
        "metadata_completion": 0,
        "bioinformatics_completion": 0,
        "academic_completion": 0,
        "tools_with_metadata": 0,
        "tools_with_github_repo": 0,
        "tools_with_citation": 0,
        "tools_by_language": defaultdict(int),
        "tools_by_license": defaultdict(int),
        "avg_description_length": 0,
        "total_stars": 0,
        "tools_updated_recently": 0,  # last 30 days
    }
    
    # Field counters
    for field in CRITICAL_FIELDS + IMPORTANT_FIELDS + METADATA_FIELDS + BIOINFORMATICS_FIELDS + ACADEMIC_FIELDS:
        metrics["fields"][field] = 0
    
    # Quality indicators
    description_lengths = []
    tool_metadata_completion = []
    critical_field_count = 0
    important_field_count = 0
    metadata_field_count = 0
    bioinformatics_field_count = 0
    academic_field_count = 0
    recent_cutoff = datetime.datetime.now() - datetime.timedelta(days=30)
    
    # Process each tool
    for tool in tools:
        # Load metadata
        tool_id = tool.get("id")
        metadata = {}
        if tool_id:
            metadata = load_metadata_for_tool(tool_id)
            if metadata:
                metrics["tools_with_metadata"] += 1
        
        # Merge tool and metadata for analysis
        combined_data = {**tool, **metadata}
        
        # Count fields
        for field in CRITICAL_FIELDS + IMPORTANT_FIELDS + METADATA_FIELDS + BIOINFORMATICS_FIELDS + ACADEMIC_FIELDS:
            if field in combined_data and combined_data[field]:
                metrics["fields"][field] += 1
                
                # Count by category
                if field in CRITICAL_FIELDS:
                    critical_field_count += 1
                if field in IMPORTANT_FIELDS:
                    important_field_count += 1
                if field in METADATA_FIELDS:
                    metadata_field_count += 1
                if field in BIOINFORMATICS_FIELDS:
                    bioinformatics_field_count += 1
                if field in ACADEMIC_FIELDS:
                    academic_field_count += 1
        
        # Repository analysis
        repo_url = combined_data.get("repository", "")
        if repo_url and "github.com" in repo_url:
            metrics["tools_with_github_repo"] += 1
        
        # Description analysis
        if "description" in combined_data and combined_data["description"]:
            description_lengths.append(len(combined_data["description"]))
        
        # Language and license counts
        language = combined_data.get("language")
        if language:
            metrics["tools_by_language"][language] += 1
            
        license_name = combined_data.get("license")
        if license_name:
            metrics["tools_by_license"][license_name] += 1
        
        # Citation analysis
        if "citation" in combined_data and combined_data["citation"]:
            metrics["tools_with_citation"] += 1
        
        # Star count
        stars = combined_data.get("stars", 0)
        if stars and isinstance(stars, (int, float)):
            metrics["total_stars"] += stars
        
        # Recent updates
        last_updated = combined_data.get("last_updated") or combined_data.get("updated_at")
        if last_updated:
            try:
                update_date = datetime.datetime.fromisoformat(last_updated.replace('Z', '+00:00'))
                if update_date >= recent_cutoff:
                    metrics["tools_updated_recently"] += 1
            except (ValueError, TypeError):
                pass
    
    # Calculate average description length
    if description_lengths:
        metrics["avg_description_length"] = sum(description_lengths) / len(description_lengths)
    
    # Calculate completion percentages
    if tools:
        metrics["critical_fields_completion"] = (critical_field_count / (len(tools) * len(CRITICAL_FIELDS))) * 100
        metrics["important_fields_completion"] = (important_field_count / (len(tools) * len(IMPORTANT_FIELDS))) * 100
        metrics["metadata_completion"] = (metadata_field_count / (len(tools) * len(METADATA_FIELDS))) * 100
        metrics["bioinformatics_completion"] = (bioinformatics_field_count / (len(tools) * len(BIOINFORMATICS_FIELDS))) * 100
        metrics["academic_completion"] = (academic_field_count / (len(tools) * len(ACADEMIC_FIELDS))) * 100
    
    # Convert defaultdicts to regular dicts for JSON serialization
    metrics["fields"] = dict(metrics["fields"])
    metrics["tools_by_language"] = dict(metrics["tools_by_language"])
    metrics["tools_by_license"] = dict(metrics["tools_by_license"])
    
    # Sort language and license by frequency
    metrics["top_languages"] = sorted(
        metrics["tools_by_language"].items(), 
        key=lambda x: x[1], 
        reverse=True
    )[:10]
    
    metrics["top_licenses"] = sorted(
        metrics["tools_by_license"].items(), 
        key=lambda x: x[1], 
        reverse=True
    )[:10]
    
    return metrics

def compare_with_history(current_metrics: Dict[str, Any], history: Dict[str, Any]) -> Dict[str, Any]:
    """
    Compare current metrics with historical data.
    
    Returns:
        Dict containing trend information
    """
    if not history.get("history"):
        return {"first_measurement": True}
    
    # Get the most recent entry
    previous = history["history"][-1]
    
    # Calculate changes in key metrics
    changes = {
        "total_tools": current_metrics["total_tools"] - previous["total_tools"],
        "critical_fields_change": current_metrics["critical_fields_completion"] - previous["critical_fields_completion"],
        "metadata_completion_change": current_metrics["metadata_completion"] - previous["metadata_completion"],
        "tools_with_citation_change": current_metrics["tools_with_citation"] - previous["tools_with_citation"],
        "tools_with_metadata_change": current_metrics["tools_with_metadata"] - previous["tools_with_metadata"],
        "tools_updated_recently_change": current_metrics["tools_updated_recently"] - previous["tools_updated_recently"],
    }
    
    # Calculate the time difference
    current_time = datetime.datetime.fromisoformat(current_metrics["timestamp"].replace('Z', '+00:00'))
    previous_time = datetime.datetime.fromisoformat(previous["timestamp"].replace('Z', '+00:00'))
    time_diff_days = (current_time - previous_time).total_seconds() / (60 * 60 * 24)
    
    changes["days_since_last_measurement"] = time_diff_days
    
    return changes

def format_metrics_report(metrics: Dict[str, Any], comparison: Dict[str, Any]) -> str:
    """Format metrics into a readable report"""
    report = ["# Data Quality Metrics Report", ""]
    report.append(f"Generated: {metrics['timestamp']}")
    report.append("")
    
    report.append("## Summary")
    report.append(f"- Total Tools: {metrics['total_tools']}")
    report.append(f"- Tools with Metadata: {metrics['tools_with_metadata']} ({metrics['tools_with_metadata']/metrics['total_tools']*100:.1f}%)")
    report.append(f"- Tools with Citations: {metrics['tools_with_citation']} ({metrics['tools_with_citation']/metrics['total_tools']*100:.1f}%)")
    report.append(f"- Tools Updated in Last 30 Days: {metrics['tools_updated_recently']} ({metrics['tools_updated_recently']/metrics['total_tools']*100:.1f}%)")
    report.append(f"- Total GitHub Stars: {metrics['total_stars']}")
    report.append("")
    
    report.append("## Data Completeness")
    report.append(f"- Critical Fields: {metrics['critical_fields_completion']:.1f}%")
    report.append(f"- Important Fields: {metrics['important_fields_completion']:.1f}%")
    report.append(f"- Metadata Fields: {metrics['metadata_completion']:.1f}%")
    report.append(f"- Bioinformatics Fields: {metrics['bioinformatics_completion']:.1f}%")
    report.append(f"- Academic Impact Fields: {metrics['academic_completion']:.1f}%")
    report.append("")
    
    report.append("## Top Languages")
    for language, count in metrics["top_languages"][:5]:
        report.append(f"- {language}: {count}")
    report.append("")
    
    report.append("## Top Licenses")
    for license_name, count in metrics["top_licenses"][:5]:
        report.append(f"- {license_name}: {count}")
    report.append("")
    
    if "first_measurement" not in comparison:
        report.append("## Changes Since Last Measurement")
        report.append(f"- Days Since Last Measurement: {comparison['days_since_last_measurement']:.1f}")
        report.append(f"- Tools Added: {comparison['total_tools']}")
        report.append(f"- Critical Field Completion: {comparison['critical_fields_change']:.1f}% change")
        report.append(f"- Metadata Completion: {comparison['metadata_completion_change']:.1f}% change")
        report.append(f"- Tools with Citations: {comparison['tools_with_citation_change']} change")
        report.append("")
    
    report.append("## Field Completion")
    
    report.append("### Most Complete Fields (>80%)")
    over_80 = [(field, count) for field, count in metrics["fields"].items() 
               if count/metrics["total_tools"] > 0.8]
    for field, count in sorted(over_80, key=lambda x: x[1], reverse=True)[:10]:
        report.append(f"- {field}: {count/metrics['total_tools']*100:.1f}%")
    report.append("")
    
    report.append("### Fields Needing Attention (<50%)")
    under_50 = [(field, count) for field, count in metrics["fields"].items() 
                if count/metrics["total_tools"] < 0.5]
    for field, count in sorted(under_50, key=lambda x: x[1])[:10]:
        report.append(f"- {field}: {count/metrics['total_tools']*100:.1f}%")
    
    return "\n".join(report)

def main():
    """Main entry point for the script."""
    parser = argparse.ArgumentParser(description="Data quality metrics collection and analysis")
    parser.add_argument("--output-file", type=str, default="metrics_report.json",
                        help="Path to save the metrics report (default: metrics_report.json)")
    parser.add_argument("--report", action="store_true",
                        help="Generate a human-readable report")
    parser.add_argument("--report-file", type=str, default="metrics_report.md",
                        help="Path to save the human-readable report (default: metrics_report.md)")
    args = parser.parse_args()
    
    # Make sure the reports directory exists
    os.makedirs(REPORTS_DIR, exist_ok=True)
    
    # Analyze data quality
    logger.info("Starting data quality analysis...")
    metrics = analyze_data_quality()
    
    if not metrics:
        logger.error("Failed to generate metrics. Exiting.")
        return 1
    
    # Load history and add current metrics
    history = load_metrics_history()
    comparison = compare_with_history(metrics, history)
    history["history"].append(metrics)
    
    # Save updated history
    save_metrics_history(history)
    
    # Save current metrics
    try:
        with open(os.path.join(REPORTS_DIR, args.output_file), 'w') as f:
            json.dump(metrics, f, indent=2)
        logger.info(f"Metrics saved to {args.output_file}")
    except (IOError, OSError) as e:
        logger.error(f"Error saving metrics: {e}")
    
    # Generate human-readable report if requested
    if args.report:
        report = format_metrics_report(metrics, comparison)
        try:
            with open(os.path.join(REPORTS_DIR, args.report_file), 'w') as f:
                f.write(report)
            logger.info(f"Report saved to {args.report_file}")
        except (IOError, OSError) as e:
            logger.error(f"Error saving report: {e}")
    
    logger.info("Data quality analysis completed successfully.")
    return 0

if __name__ == "__main__":
    sys.exit(main())