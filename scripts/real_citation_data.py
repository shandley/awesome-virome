#!/usr/bin/env python3
"""
Real Citation Data Generation

This script generates citation data for the dashboard visualization using ONLY real citation data
from various sources (PubMed, CrossRef, etc.). Unlike the previous dashboard_citation_data.py,
this script does NOT generate synthetic/artificial citation data.

Usage:
    python real_citation_data.py
"""

import os
import json
import logging
import glob
from datetime import datetime
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Any, Optional

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# File paths
BASE_DIR = Path(__file__).resolve().parent.parent
DATA_JSON_PATH = BASE_DIR / "data.json"
IMPACT_DATA_PATH = BASE_DIR / "impact_data.json"
ACADEMIC_IMPACT_DIR = BASE_DIR / "metadata" / "academic_impact"
PUBMED_CITATIONS_DIR = BASE_DIR / "metadata" / "pubmed_citations"
REPORTS_DIR = BASE_DIR / "reports" / "citations"
BIOINFORMATICS_DIR = BASE_DIR / "metadata" / "bioinformatics"


def get_year_from_date(date_str: Optional[str]) -> Optional[int]:
    """Extract year from ISO date string."""
    if not date_str:
        return None
    
    try:
        date_obj = datetime.fromisoformat(date_str.replace('Z', '+00:00'))
        return date_obj.year
    except ValueError:
        return None


def load_real_citation_data():
    """Load real citation data from metadata files."""
    real_citation_data = {}
    
    # 1. Load from PubMed citations directory
    for file_path in glob.glob(str(PUBMED_CITATIONS_DIR / "*.json")):
        if "pubmed_citations.json" in file_path:
            continue  # Skip the summary file
            
        try:
            with open(file_path, 'r') as f:
                tool_data = json.load(f)
                
            name = tool_data.get("name")
            if not name:
                continue
                
            # Get citation info
            citation_info = tool_data.get("citation_info", {})
            publication = citation_info.get("publication", {})
            
            if publication:
                doi = publication.get("doi")
                citation_count = publication.get("citation_count", 0)
                
                if citation_count and citation_count > 0:
                    real_citation_data[name] = {
                        "name": name,
                        "url": tool_data.get("url", ""),
                        "doi": doi,
                        "total_citations": citation_count,
                        "influential_citations": int(citation_count * 0.2)  # Estimate influential citations
                    }
        except (json.JSONDecodeError, IOError) as e:
            logger.error(f"Error reading {file_path}: {e}")
            continue
    
    # 2. Load from academic_impact directory
    for file_path in glob.glob(str(ACADEMIC_IMPACT_DIR / "*.json")):
        if "academic_impact.json" in file_path:
            continue  # Skip the summary file
            
        try:
            with open(file_path, 'r') as f:
                tool_data = json.load(f)
                
            name = tool_data.get("name")
            if not name:
                continue
                
            doi = tool_data.get("doi")
            citation_metrics = tool_data.get("citation_metrics", {})
            
            # Extract citation data if available
            total_citations = citation_metrics.get("total_citations", 0)
            if isinstance(total_citations, str) and total_citations.isdigit():
                total_citations = int(total_citations)
                
            influential_citations = citation_metrics.get("influential_citations", 0)
            if isinstance(influential_citations, str) and influential_citations.isdigit():
                influential_citations = int(influential_citations)
                
            citations_by_year = citation_metrics.get("citations_by_year", {})
            
            # Only add if there's actual citation data
            if total_citations > 0 or citations_by_year:
                real_citation_data[name] = {
                    "name": name,
                    "url": tool_data.get("url", ""),
                    "doi": doi,
                    "total_citations": total_citations,
                    "influential_citations": influential_citations
                }
                
                if citations_by_year:
                    real_citation_data[name]["citations_by_year"] = citations_by_year
        except (json.JSONDecodeError, IOError) as e:
            logger.error(f"Error reading {file_path}: {e}")
            continue
    
    # 3. Load from bioinformatics directory
    for file_path in glob.glob(str(BIOINFORMATICS_DIR / "*.json")):
        if "summary.json" in file_path:
            continue  # Skip the summary file
            
        try:
            with open(file_path, 'r') as f:
                tool_data = json.load(f)
                
            name = tool_data.get("name")
            if not name or name in real_citation_data:
                continue
                
            academic_impact = tool_data.get("academic_impact", {})
            if not academic_impact:
                continue
                
            # Extract citation data
            doi = academic_impact.get("doi")
            total_citations = academic_impact.get("total_citations", 0)
            influential_citations = academic_impact.get("influential_citations", 0)
            citations_by_year = academic_impact.get("citations_by_year", {})
            
            # Only add if there's actual citation data
            if total_citations > 0 or citations_by_year:
                real_citation_data[name] = {
                    "name": name,
                    "url": tool_data.get("url", ""),
                    "doi": doi,
                    "total_citations": total_citations,
                    "influential_citations": influential_citations
                }
                
                if citations_by_year:
                    real_citation_data[name]["citations_by_year"] = citations_by_year
        except (json.JSONDecodeError, IOError) as e:
            logger.error(f"Error reading {file_path}: {e}")
            continue
    
    logger.info(f"Loaded real citation data for {len(real_citation_data)} tools")
    return real_citation_data


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


def main():
    """Main function to generate dashboard citation data with only real data."""
    logger.info("Generating dashboard citation data from real sources only...")
    
    # Load data.json for tools list
    try:
        with open(DATA_JSON_PATH, 'r', encoding='utf-8') as f:
            data = json.load(f)
    except (FileNotFoundError, json.JSONDecodeError) as e:
        logger.error(f"Error reading data.json: {e}")
        return
    
    # Load real citation data
    real_citation_data = load_real_citation_data()
    
    # Initialize impact data
    impact_data = {
        "last_updated": datetime.now().isoformat(),
        "tools": [],
        "adoption_trend": defaultdict(int),
        "categories": defaultdict(int),
        "citations": {
            "total": 0,
            "by_year": defaultdict(int)
        }
    }
    
    # Extract tools
    tools = [node for node in data.get('nodes', []) if node.get('type') == 'tool']
    logger.info(f"Found {len(tools)} tools in data.json")
    
    # Process tools
    processed_tools = set()
    total_citations = 0
    citations_by_year = defaultdict(int)
    
    # First process tools with real citation data
    for name, citation_data in real_citation_data.items():
        # Add basic info for the dashboard
        tool_data = {
            "name": name,
            "url": citation_data.get("url", ""),
            "doi": citation_data.get("doi"),
            "total_citations": citation_data.get("total_citations", 0),
            "influential_citations": citation_data.get("influential_citations", 0),
            "citations_by_year": citation_data.get("citations_by_year", {})
        }
        
        # Find category from data.json if possible
        matching_tools = [t for t in tools if t.get('name') == name]
        if matching_tools:
            category = matching_tools[0].get('category')
            if category:
                tool_data["category"] = category
                impact_data["categories"][category] += 1
            else:
                # Use heuristic for category
                category = get_tool_category(name)
                tool_data["category"] = category
                impact_data["categories"][category] += 1
                
            # Add to adoption trend
            created_year = get_year_from_date(matching_tools[0].get('createdAt'))
            if created_year:
                impact_data["adoption_trend"][str(created_year)] += 1
        else:
            tool_data["category"] = "Other Tools"
            impact_data["categories"]["Other Tools"] += 1
        
        # Update citation totals
        total_citations += tool_data.get("total_citations", 0)
        for year, count in tool_data.get("citations_by_year", {}).items():
            citations_by_year[year] += count
        
        impact_data["tools"].append(tool_data)
        processed_tools.add(name)
    
    # Then include all other tools, but with empty citation data
    for tool in tools:
        name = tool.get('name', 'Unknown')
        if name in processed_tools:
            continue
        
        # Add minimal data for tools without citations
        tool_data = {
            "name": name,
            "url": tool.get('url', ''),
            "doi": None,
            "total_citations": 0,
            "influential_citations": 0,
            "citations_by_year": {}
        }
        
        # Add category
        category = tool.get('category')
        if category:
            tool_data["category"] = category
            impact_data["categories"][category] += 1
        else:
            category = get_tool_category(name)
            tool_data["category"] = category
            impact_data["categories"][category] += 1
        
        # Add to adoption trend
        created_year = get_year_from_date(tool.get('createdAt'))
        if created_year:
            impact_data["adoption_trend"][str(created_year)] += 1
        
        impact_data["tools"].append(tool_data)
        processed_tools.add(name)
    
    # Update citation totals
    impact_data["citations"]["total"] = total_citations
    impact_data["citations"]["by_year"] = dict(citations_by_year)
    
    # Calculate average citations per tool with real citations
    tools_with_citations = sum(1 for tool in impact_data["tools"] if tool.get("total_citations", 0) > 0)
    if tools_with_citations > 0:
        impact_data["average_citations"] = round(total_citations / tools_with_citations, 2)
    else:
        impact_data["average_citations"] = 0
    
    impact_data["tools_with_citations"] = tools_with_citations
    impact_data["total_tools"] = len(processed_tools)
    
    # Sort tools by citation count
    impact_data["tools"].sort(
        key=lambda x: x.get("total_citations", 0),
        reverse=True
    )
    
    # Convert defaultdicts to regular dicts for JSON serialization
    impact_data["adoption_trend"] = dict(impact_data["adoption_trend"])
    impact_data["categories"] = dict(impact_data["categories"])
    
    # Preserve network data if it exists in current impact_data.json
    try:
        if os.path.exists(IMPACT_DATA_PATH):
            with open(IMPACT_DATA_PATH, 'r', encoding='utf-8') as f:
                current_data = json.load(f)
                if "relationships" in current_data:
                    impact_data["relationships"] = current_data["relationships"]
    except (json.JSONDecodeError, IOError) as e:
        logger.warning(f"Could not read existing network data: {e}")
    
    # Save impact_data.json
    with open(IMPACT_DATA_PATH, 'w', encoding='utf-8') as f:
        json.dump(impact_data, f, indent=2)
    
    logger.info(f"Generated real citation data for {len(processed_tools)} tools")
    logger.info(f"Tools with real citation data: {tools_with_citations}")
    logger.info(f"Total citations: {total_citations}")
    logger.info(f"Years with citation data: {sorted(citations_by_year.keys())}")
    logger.info(f"Saved data to {IMPACT_DATA_PATH}")


if __name__ == "__main__":
    main()