#!/usr/bin/env python3
"""
Generate statistics about citation sources.

This script analyzes the citation sources in impact_data.json
and generates detailed statistics about where the citation data is coming from.

Usage:
    python citation_statistics.py
"""

import os
import json
import logging
from pathlib import Path
from collections import defaultdict, Counter
from typing import Dict, Any, List, Optional

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Constants
BASE_DIR = Path(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
IMPACT_DATA_PATH = BASE_DIR / "impact_data.json"


def load_impact_data() -> Dict[str, Any]:
    """Load the impact_data.json file."""
    try:
        with open(IMPACT_DATA_PATH, 'r') as f:
            data = json.load(f)
        return data
    except (FileNotFoundError, json.JSONDecodeError) as e:
        logger.error(f"Error loading impact data: {e}")
        return {}


def analyze_citation_sources() -> Dict[str, Any]:
    """Analyze citation sources in impact_data.json."""
    # Load impact data
    impact_data = load_impact_data()
    if not impact_data:
        logger.error("Failed to load impact data")
        return {}
    
    # Get tools
    tools = impact_data.get('tools', [])
    if not tools:
        logger.error("No tools found in impact_data.json")
        return {}
    
    # Initialize statistics
    stats = {
        "total_tools": len(tools),
        "tools_with_citations": 0,
        "total_citations": impact_data.get("total_citations", 0),
        "citation_sources": Counter(),
        "yearly_citation_sources": Counter(),
        "citation_counts_by_source": defaultdict(list),
        "category_sources": defaultdict(Counter),
        "synthetic_data": 0,
        "real_data": 0,
        "top_tools": []
    }
    
    # Process all tools
    for tool in tools:
        name = tool.get('name', 'Unknown')
        citations = tool.get('total_citations', 0)
        
        if citations <= 0:
            continue
            
        stats["tools_with_citations"] += 1
        
        # Track citation sources
        source = tool.get('citation_source', 'unknown')
        yearly_source = tool.get('yearly_citation_source', 'unknown')
        
        stats["citation_sources"][source] += 1
        stats["yearly_citation_sources"][yearly_source] += 1
        
        # Track citations by source
        stats["citation_counts_by_source"][source].append(citations)
        
        # Track by category
        category = tool.get('category', 'Other Tools')
        stats["category_sources"][category][source] += 1
        
        # Track synthetic vs real
        if source == 'synthetic':
            stats["synthetic_data"] += 1
        else:
            stats["real_data"] += 1
        
        # Add to top tools
        stats["top_tools"].append({
            "name": name,
            "citations": citations,
            "source": source,
            "category": category
        })
    
    # Sort top tools by citations
    stats["top_tools"].sort(key=lambda x: x["citations"], reverse=True)
    stats["top_tools"] = stats["top_tools"][:10]  # Keep top 10
    
    # Calculate average citations per source
    stats["avg_citations_by_source"] = {}
    for source, counts in stats["citation_counts_by_source"].items():
        if counts:
            stats["avg_citations_by_source"][source] = sum(counts) / len(counts)
    
    return stats


def print_statistics(stats: Dict[str, Any]):
    """Print statistics about citation sources."""
    print("\n===== CITATION SOURCE STATISTICS =====\n")
    
    tools_with_citations = stats["tools_with_citations"]
    total_tools = stats["total_tools"]
    percent_with_citations = tools_with_citations / total_tools * 100 if total_tools > 0 else 0
    
    print(f"Total tools: {total_tools}")
    print(f"Tools with citations: {tools_with_citations} ({percent_with_citations:.1f}%)")
    print(f"Total citation count: {stats['total_citations']}")
    
    avg_citations = stats['total_citations'] / max(1, tools_with_citations)
    print(f"Average citations per tool: {avg_citations:.1f}")
    
    print("\n--- Real vs Synthetic Data ---")
    real_percent = stats["real_data"] / max(1, tools_with_citations) * 100
    synthetic_percent = stats["synthetic_data"] / max(1, tools_with_citations) * 100
    print(f"Real citation data: {stats['real_data']} tools ({real_percent:.1f}%)")
    print(f"Synthetic citation data: {stats['synthetic_data']} tools ({synthetic_percent:.1f}%)")
    
    print("\n--- Citation Sources ---")
    for source, count in stats["citation_sources"].most_common():
        avg = stats["avg_citations_by_source"].get(source, 0)
        percent = count / max(1, tools_with_citations) * 100
        print(f"{source}: {count} tools ({percent:.1f}%), avg {avg:.1f} citations/tool")
    
    print("\n--- Yearly Citation Sources ---")
    for source, count in stats["yearly_citation_sources"].most_common():
        percent = count / max(1, tools_with_citations) * 100
        print(f"{source}: {count} tools ({percent:.1f}%)")
    
    print("\n--- Categories by Source ---")
    for category, sources in sorted(stats["category_sources"].items(), 
                                     key=lambda x: sum(x[1].values()), 
                                     reverse=True):
        total = sum(sources.values())
        print(f"{category}: {total} tools")
        for source, count in sources.most_common(3):  # Top 3 sources
            source_percent = count / total * 100
            print(f"  - {source}: {count} tools ({source_percent:.1f}%)")
    
    print("\n--- Top 10 Most Cited Tools ---")
    for i, tool in enumerate(stats["top_tools"]):
        print(f"{i+1}. {tool['name']}: {tool['citations']} citations (Source: {tool['source']}, Category: {tool['category']})")


def main():
    """Main function to run the analysis."""
    logger.info("Analyzing citation statistics in impact_data.json...")
    
    # Analyze citation sources
    stats = analyze_citation_sources()
    if not stats:
        logger.error("Failed to analyze citation statistics")
        return 1
    
    # Print statistics
    print_statistics(stats)
    
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())