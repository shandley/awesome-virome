#!/usr/bin/env python3
"""
Generate Citation Data for Awesome-Virome Dashboard

This script generates realistic sample citation data for all tools
in the repository, suitable for visualization in the dashboard.
It doesn't rely on external APIs or packages.

Usage:
    python generate_citation_data.py
"""

import os
import json
import random
import time
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

def generate_citation_metrics(name: str, created_year: Optional[int] = None) -> Dict[str, Any]:
    """Generate realistic sample citation metrics for a tool."""
    # Use hash of name for consistent but random-looking numbers
    seed = sum(ord(c) for c in name)
    random.seed(seed)
    
    # Determine start year for citations
    current_year = datetime.now().year
    default_start = max(2014, current_year - 10)
    start_year = created_year if created_year and created_year >= 2014 else default_start
    start_year = min(start_year, current_year - 2)  # Ensure at least 2 years of data
    
    # Generate years
    years = list(range(start_year, current_year + 1))
    
    # Base citations depends on name length (just for variety)
    base_citations = max(3, min(20, len(name) + random.randint(1, 10)))
    
    # Generate citation pattern
    citations_by_year = {}
    total_citations = 0
    
    # Different citation growth patterns based on hash
    pattern_type = seed % 5
    
    for i, year in enumerate(years):
        year_str = str(year)
        
        if pattern_type == 0:
            # Linear growth
            citations = base_citations + i * random.randint(1, 3)
        elif pattern_type == 1:
            # Exponential growth
            citations = base_citations * (1.2 + 0.1 * random.random()) ** i
        elif pattern_type == 2:
            # Early spike then plateau
            if i < len(years) // 3:
                citations = base_citations * (i + 1)
            else:
                citations = base_citations * (len(years) // 3) * (1 + 0.1 * random.random())
        elif pattern_type == 3:
            # Slow start, rapid middle growth, then slower end
            if i < len(years) // 3:
                citations = base_citations * (1 + 0.2 * i)
            elif i < 2 * len(years) // 3:
                citations = base_citations * (1 + 0.5 * i)
            else:
                citations = base_citations * (1 + 0.3 * i)
        else:
            # Steady citations
            citations = base_citations * (1 + 0.1 * random.random())
        
        # Add some randomness
        citations = int(citations * (0.8 + 0.4 * random.random()))
        
        # For newest tools, reduce citations
        if created_year and created_year >= current_year - 3:
            citations = int(citations * 0.4)
        
        # For most recent year, reduce citations (partial year)
        if year == current_year:
            citations = int(citations * 0.4)
        
        citations_by_year[year_str] = max(1, citations)
        total_citations += citations_by_year[year_str]
    
    # Calculate influential citations
    influential_citations = int(total_citations * (0.1 + 0.1 * random.random()))
    
    return {
        'citations_by_year': citations_by_year,
        'influential_citations': influential_citations,
        'total_citations': total_citations
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
    """Main function to generate citation data."""
    logger.info("Generating citation data for dashboard visualization...")
    
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
        
        # Generate citation metrics
        metrics = generate_citation_metrics(name, created_year)
        
        # Add to citation totals
        total_citations += metrics['total_citations']
        for year, count in metrics['citations_by_year'].items():
            citations_by_year[year] += count
        
        # Add to tools list
        impact_data["tools"].append({
            "name": name,
            "citations_by_year": metrics['citations_by_year'],
            "influential_citations": metrics['influential_citations']
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
    
    logger.info(f"Generated citation data for {len(processed_tools)} tools")
    logger.info(f"Total citations: {total_citations}")
    logger.info(f"Years with citation data: {sorted(citations_by_year.keys())}")
    logger.info(f"Saved data to {IMPACT_DATA_PATH}")

if __name__ == "__main__":
    main()