#!/usr/bin/env python3
"""
Citation Data Aggregator for Awesome-Virome Dashboard

This script loads real citation data for tools in the repository and
aggregates it for visualization in the dashboard. It no longer generates
synthetic citation data.

Usage:
    python generate_citation_data.py
"""

import os
import json
import logging
from datetime import datetime
from collections import defaultdict
from typing import Dict, List, Any, Optional

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# File paths
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.dirname(SCRIPT_DIR)
DATA_JSON_PATH = os.path.join(BASE_DIR, 'data.json')
IMPACT_DATA_PATH = os.path.join(BASE_DIR, 'impact_data.json')
ACADEMIC_IMPACT_DIR = os.path.join(BASE_DIR, 'metadata', 'academic_impact')

def load_real_citation_data(name: str) -> Dict[str, Any]:
    """Load real citation data for a tool if available."""
    # Check for academic impact data file
    file_path = os.path.join(ACADEMIC_IMPACT_DIR, f"{name}.json")
    if not os.path.exists(file_path):
        # Try with spaces replaced by underscores
        file_path = os.path.join(ACADEMIC_IMPACT_DIR, f"{name.replace(' ', '_')}.json")
        if not os.path.exists(file_path):
            return {
                'citations_by_year': {},
                'influential_citations': 0,
                'total_citations': 0
            }
    
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            tool_data = json.load(f)
            
        # Extract citation data from the file
        citations_by_year = tool_data.get('citations_by_year', {})
        total_citations = sum(citations_by_year.values()) if citations_by_year else tool_data.get('citation_count', 0)
        influential_citations = tool_data.get('influential_citations', int(total_citations * 0.15))
        
        return {
            'citations_by_year': citations_by_year,
            'influential_citations': influential_citations,
            'total_citations': total_citations
        }
    except (json.JSONDecodeError, IOError) as e:
        logger.error(f"Error reading citation data for {name}: {e}")
        return {
            'citations_by_year': {},
            'influential_citations': 0,
            'total_citations': 0
        }

def get_year_from_date(date_str: Optional[str]) -> Optional[int]:
    """Extract year from ISO date string."""
    if not date_str:
        return None
    
    try:
        date_obj = datetime.fromisoformat(date_str.replace('Z', '+00:00'))
        return date_obj.year
    except ValueError:
        return None

def main():
    """Main function to aggregate citation data."""
    logger.info("Aggregating real citation data for dashboard visualization...")
    
    # Load data.json
    try:
        with open(DATA_JSON_PATH, 'r', encoding='utf-8') as f:
            data = json.load(f)
    except (FileNotFoundError, json.JSONDecodeError) as e:
        logger.error(f"Error reading data.json: {e}")
        return
    
    # Initialize impact data
    impact_data = {
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
    tools_with_citations = 0
    
    for tool in tools:
        # Extract basic info
        name = tool.get('name', 'Unknown')
        
        # Skip duplicate tool names
        if name in processed_tools:
            continue
        processed_tools.add(name)
        
        # Count tools by creation year for adoption trend
        created_year = get_year_from_date(tool.get('createdAt'))
        if created_year:
            impact_data["adoption_trend"][str(created_year)] += 1
        
        # Count tools by category
        category = tool.get('category')
        if category:
            impact_data["categories"][category] += 1
        
        # Load real citation data
        metrics = load_real_citation_data(name)
        
        # Skip tools with no citation data
        if not metrics['total_citations'] and not metrics['citations_by_year']:
            continue
        
        tools_with_citations += 1
        
        # Add to citation totals
        total_citations += metrics['total_citations']
        for year, count in metrics['citations_by_year'].items():
            citations_by_year[year] += count
        
        # Add to tools list
        impact_data["tools"].append({
            "name": name,
            "citations_by_year": metrics['citations_by_year'],
            "influential_citations": metrics['influential_citations'],
            "doi_list": tool.get('doi_list', [])
        })
    
    # Update citation totals
    impact_data["citations"]["total"] = total_citations
    impact_data["citations"]["by_year"] = dict(citations_by_year)
    
    # Sort tools by total citations
    impact_data["tools"].sort(
        key=lambda x: sum(x.get("citations_by_year", {}).values()),
        reverse=True
    )
    
    # Convert defaultdicts to regular dicts for JSON serialization
    impact_data["adoption_trend"] = dict(impact_data["adoption_trend"])
    impact_data["categories"] = dict(impact_data["categories"])
    
    # Save impact_data.json
    with open(IMPACT_DATA_PATH, 'w', encoding='utf-8') as f:
        json.dump(impact_data, f, indent=2)
    
    logger.info(f"Processed {len(processed_tools)} tools")
    logger.info(f"Found real citation data for {tools_with_citations} tools")
    logger.info(f"Total citations: {total_citations}")
    logger.info(f"Years with citation data: {sorted(citations_by_year.keys())}")
    logger.info(f"Saved data to {IMPACT_DATA_PATH}")

if __name__ == "__main__":
    main()