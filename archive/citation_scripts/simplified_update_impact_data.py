#!/usr/bin/env python3
"""
Simplified script to update impact_data.json with real citation data.

This script combines citation data from multiple sources:
1. Most cited tools from most_cited_report.json
2. Citation trends over time from citation_trends_report.json

The generated impact_data.json file is used by the citation dashboard
to display real citation metrics instead of mock data.
"""

import json
import os
import sys
from datetime import datetime
from pathlib import Path

# Set up paths
BASE_DIR = Path(__file__).resolve().parent.parent
REPORTS_DIR = BASE_DIR / "reports" / "citations"
OUTPUT_FILE = BASE_DIR / "impact_data.json"

# Input files
MOST_CITED_FILE = REPORTS_DIR / "most_cited_report.json"
TRENDS_FILE = REPORTS_DIR / "citation_trends_report.json"
SUMMARY_FILE = REPORTS_DIR / "citation_reports_summary.json"

def main():
    """Generate impact_data.json from real citation data."""
    print(f"Updating impact_data.json with real citation data...")
    
    # Check if input files exist
    for file_path in [MOST_CITED_FILE, TRENDS_FILE, SUMMARY_FILE]:
        if not file_path.exists():
            print(f"Error: Required file not found: {file_path}")
            sys.exit(1)
    
    try:
        # Load citation data from reports
        with open(MOST_CITED_FILE, 'r') as f:
            most_cited_data = json.load(f)
        
        with open(TRENDS_FILE, 'r') as f:
            trends_data = json.load(f)
            
        with open(SUMMARY_FILE, 'r') as f:
            summary_data = json.load(f)
        
        # Create impact data structure
        impact_data = {
            "last_updated": datetime.now().isoformat(),
            "total_tools": summary_data.get("total_tools_analyzed", 0),
            "tools_with_citations": summary_data.get("total_tools_with_citations", 0),
            "total_citations": summary_data.get("total_citations", 0),
            "average_citations": most_cited_data.get("average_citations_per_tool", 0),
            "tools": []
        }
        
        # Process tool data
        for tool in most_cited_data.get("most_cited_tools", []):
            # Get yearly citation data if available
            citations_by_year = {}
            if tool["name"] in trends_data.get("tool_yearly_data", {}):
                citations_by_year = trends_data["tool_yearly_data"][tool["name"]]
            
            # Create tool entry
            tool_entry = {
                "name": tool["name"],
                "url": tool.get("url", ""),
                "doi": tool.get("doi", ""),
                "total_citations": tool.get("total_citations", 0),
                "influential_citations": tool.get("influential_citations", 0),
                "citations_by_year": citations_by_year,
                "category": get_tool_category(tool["name"])
            }
            
            impact_data["tools"].append(tool_entry)
        
        # Write to impact_data.json
        with open(OUTPUT_FILE, 'w') as f:
            json.dump(impact_data, f, indent=2)
        
        print(f"Successfully updated {OUTPUT_FILE}")
        print(f"Included {len(impact_data['tools'])} tools with citation data")
        print(f"Total citations: {impact_data['total_citations']}")
        
    except Exception as e:
        print(f"Error updating impact_data.json: {e}")
        sys.exit(1)

def get_tool_category(tool_name):
    """
    Determine category for a tool based on its name.
    """
    # Simple heuristic based on tool name
    if any(term in tool_name.lower() for term in ["phage", "prophage"]):
        return "Phage Analysis"
    elif "host" in tool_name.lower():
        return "Host Prediction"
    elif any(term in tool_name.lower() for term in ["sort", "class", "tax"]):
        return "Taxonomy"
    elif any(term in tool_name.lower() for term in ["find", "ident", "detect"]):
        return "Virus Identification"
    else:
        return "Other Tools"

if __name__ == "__main__":
    main()