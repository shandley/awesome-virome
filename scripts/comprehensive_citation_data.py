#!/usr/bin/env python3
"""
Comprehensive script to update impact_data.json with citation data for all tools.

This script gathers citation data from multiple sources:
1. Academic impact files in metadata/academic_impact/
2. Bioinformatics metadata files in metadata/bioinformatics/
3. Existing citation reports in reports/citations/

The goal is to create a comprehensive impact_data.json file that includes
citation data for all tools, not just the most-cited ones.
"""

import json
import os
import sys
import glob
from datetime import datetime
from pathlib import Path
from collections import defaultdict

# Set up paths
BASE_DIR = Path(__file__).resolve().parent.parent
ACADEMIC_IMPACT_DIR = BASE_DIR / "metadata" / "academic_impact"
BIOINFORMATICS_DIR = BASE_DIR / "metadata" / "bioinformatics"
REPORTS_DIR = BASE_DIR / "reports" / "citations"
OUTPUT_FILE = BASE_DIR / "impact_data.json"

def main():
    """Generate comprehensive impact_data.json with citation data for all tools."""
    print(f"Updating impact_data.json with comprehensive citation data...")
    
    # Collect citation data from all possible sources
    all_tools_data = collect_all_citation_data()
    
    # Calculate summary metrics
    total_citations = sum(tool.get("total_citations", 0) for tool in all_tools_data.values())
    tools_with_citations = sum(1 for tool in all_tools_data.values() if tool.get("total_citations", 0) > 0)
    
    if tools_with_citations > 0:
        average_citations = round(total_citations / tools_with_citations, 2)
    else:
        average_citations = 0
    
    # Create impact data structure
    impact_data = {
        "last_updated": datetime.now().isoformat(),
        "total_tools": len(all_tools_data),
        "tools_with_citations": tools_with_citations,
        "total_citations": total_citations,
        "average_citations": average_citations,
        "tools": []
    }
    
    # Add tool entries to impact_data
    for tool_name, tool_data in all_tools_data.items():
        # Only include tools with citation data in the final output
        if tool_data.get("total_citations", 0) > 0 or tool_data.get("citations_by_year"):
            impact_data["tools"].append(tool_data)
    
    # Sort tools by citation count (descending)
    impact_data["tools"].sort(key=lambda x: x.get("total_citations", 0), reverse=True)
    
    # Write to impact_data.json
    with open(OUTPUT_FILE, 'w') as f:
        json.dump(impact_data, f, indent=2)
    
    print(f"Successfully updated {OUTPUT_FILE}")
    print(f"Included {len(impact_data['tools'])} tools with citation data")
    print(f"Total citations: {impact_data['total_citations']}")

def collect_all_citation_data():
    """Collect citation data from all available sources."""
    all_tools_data = {}
    
    # 1. Collect data from most_cited_report.json and citation_trends_report.json
    collect_from_citation_reports(all_tools_data)
    
    # 2. Collect data from metadata/academic_impact/ directory
    collect_from_academic_impact(all_tools_data)
    
    # 3. Collect data from metadata/bioinformatics/ directory
    collect_from_bioinformatics(all_tools_data)
    
    return all_tools_data

def collect_from_citation_reports(all_tools_data):
    """Collect citation data from citation report files."""
    # Check if report files exist
    most_cited_file = REPORTS_DIR / "most_cited_report.json"
    trends_file = REPORTS_DIR / "citation_trends_report.json"
    
    if most_cited_file.exists() and trends_file.exists():
        # Load most cited tools
        with open(most_cited_file, 'r') as f:
            most_cited_data = json.load(f)
        
        # Load citation trends
        with open(trends_file, 'r') as f:
            trends_data = json.load(f)
        
        # Process most cited tools
        for tool in most_cited_data.get("most_cited_tools", []):
            name = tool["name"]
            
            # Create or update tool entry
            if name not in all_tools_data:
                all_tools_data[name] = {
                    "name": name,
                    "url": tool.get("url", ""),
                    "doi": tool.get("doi", ""),
                    "total_citations": tool.get("total_citations", 0),
                    "influential_citations": tool.get("influential_citations", 0),
                    "category": get_tool_category(name)
                }
            
            # Get yearly citation data if available
            if name in trends_data.get("tool_yearly_data", {}):
                all_tools_data[name]["citations_by_year"] = trends_data["tool_yearly_data"][name]

def collect_from_academic_impact(all_tools_data):
    """Collect citation data from academic_impact JSON files."""
    # Get all JSON files in academic_impact directory
    academic_files = glob.glob(str(ACADEMIC_IMPACT_DIR / "*.json"))
    
    for file_path in academic_files:
        try:
            with open(file_path, 'r') as f:
                tool_data = json.load(f)
                
            name = tool_data.get("name")
            if not name:
                continue
                
            doi = tool_data.get("doi")
            citation_metrics = tool_data.get("citation_metrics", {})
            
            # Get citation counts
            total_citations = citation_metrics.get("total_citations", 0)
            if isinstance(total_citations, str) and total_citations.isdigit():
                total_citations = int(total_citations)
            elif not isinstance(total_citations, (int, float)):
                total_citations = 0
                
            influential_citations = citation_metrics.get("influential_citations", 0)
            if isinstance(influential_citations, str) and influential_citations.isdigit():
                influential_citations = int(influential_citations)
            elif not isinstance(influential_citations, (int, float)):
                influential_citations = 0
                
            # Get citations by year
            citations_by_year = citation_metrics.get("citations_by_year", {})
            if not isinstance(citations_by_year, dict):
                citations_by_year = {}
            
            # Create or update tool entry
            if name not in all_tools_data:
                all_tools_data[name] = {
                    "name": name,
                    "url": tool_data.get("url", ""),
                    "doi": doi,
                    "category": get_tool_category(name)
                }
            
            # Update citation data if available
            if total_citations > 0:
                all_tools_data[name]["total_citations"] = total_citations
                
            if influential_citations > 0:
                all_tools_data[name]["influential_citations"] = influential_citations
                
            if citations_by_year:
                all_tools_data[name]["citations_by_year"] = citations_by_year
                
        except (json.JSONDecodeError, IOError) as e:
            print(f"Error reading {file_path}: {e}")
            continue

def collect_from_bioinformatics(all_tools_data):
    """Collect citation data from bioinformatics metadata files."""
    # Get all JSON files in bioinformatics directory
    bioinformatics_files = glob.glob(str(BIOINFORMATICS_DIR / "*.json"))
    
    for file_path in bioinformatics_files:
        try:
            with open(file_path, 'r') as f:
                tool_data = json.load(f)
                
            name = tool_data.get("name")
            if not name:
                continue
                
            academic_impact = tool_data.get("academic_impact", {})
            if not academic_impact:
                continue
                
            doi = academic_impact.get("doi")
            total_citations = academic_impact.get("total_citations", 0)
            influential_citations = academic_impact.get("influential_citations", 0)
            citations_by_year = academic_impact.get("citations_by_year", {})
            
            # Create or update tool entry
            if name not in all_tools_data:
                all_tools_data[name] = {
                    "name": name,
                    "url": tool_data.get("url", ""),
                    "doi": doi,
                    "category": get_tool_category(name)
                }
            
            # Update citation data if available and not already set
            if total_citations > 0 and "total_citations" not in all_tools_data[name]:
                all_tools_data[name]["total_citations"] = total_citations
                
            if influential_citations > 0 and "influential_citations" not in all_tools_data[name]:
                all_tools_data[name]["influential_citations"] = influential_citations
                
            if citations_by_year and "citations_by_year" not in all_tools_data[name]:
                all_tools_data[name]["citations_by_year"] = citations_by_year
                
        except (json.JSONDecodeError, IOError) as e:
            print(f"Error reading {file_path}: {e}")
            continue

def get_tool_category(tool_name):
    """
    Determine category for a tool based on its name.
    """
    # Normalize tool name for comparison
    tool_lower = tool_name.lower()
    
    # Simple heuristic based on tool name
    if any(term in tool_lower for term in ["phage", "prophage", "phate"]):
        return "Phage Analysis"
    elif "host" in tool_lower:
        return "Host Prediction"
    elif any(term in tool_lower for term in ["sort", "class", "tax", "bertax"]):
        return "Taxonomy"
    elif any(term in tool_lower for term in ["find", "ident", "detect", "seeker"]):
        return "Virus Identification"
    elif any(term in tool_lower for term in ["assembl", "spades"]):
        return "Genome Assembly"
    elif any(term in tool_lower for term in ["annot", "gene", "protein"]):
        return "Genome Annotation"
    elif any(term in tool_lower for term in ["meta", "community"]):
        return "Metagenomics"
    elif any(term in tool_lower for term in ["struct", "fold", "alphafold"]):
        return "Structural Analysis"
    else:
        return "Other Tools"

if __name__ == "__main__":
    main()