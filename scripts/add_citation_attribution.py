#!/usr/bin/env python3
"""
Add citation source attribution to impact_data.json.

This script adds citation source attribution fields to impact_data.json if they're missing.
It's designed to be run before generating visualizations to ensure source attribution
data is available regardless of how the impact data was generated.

Usage:
    python add_citation_attribution.py
"""

import os
import json
import logging
from pathlib import Path
from typing import Dict, Any, List

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# Constants
DEFAULT_SOURCE = "pubmed"  # Default source for legacy data
IMPACT_DATA_PATH = "impact_data.json"
KNOWN_SOURCES = ["scopus", "wos", "icite", "pubmed", "crossref"]


def add_attribution_to_impact_data(impact_data_path: str = IMPACT_DATA_PATH) -> bool:
    """
    Add citation source attribution to impact_data.json.
    
    Args:
        impact_data_path: Path to the impact_data.json file
        
    Returns:
        True if changes were made, False otherwise
    """
    logger.info(f"Adding attribution to {impact_data_path}")
    
    try:
        # Load impact data
        with open(impact_data_path, 'r') as f:
            impact_data = json.load(f)
        
        # Track changes
        changed = False
        citation_sources = {}
        yearly_citation_sources = {}
        
        # Process all tools
        tools = impact_data.get('tools', [])
        logger.info(f"Processing {len(tools)} tools in impact data")
        tools_with_citations = 0
        tools_needing_attribution = 0
        
        for tool in tools:
            if tool.get('total_citations', 0) <= 0:
                continue
                
            tools_with_citations += 1
            
            # Add citation source attribution
            if 'citation_source' not in tool:
                # Determine source
                source = None
                
                # Try to get source from existing fields
                if 'primary_source' in tool:
                    source = tool['primary_source']
                else:
                    # Default for tools from the pubmed_citations directory
                    source = DEFAULT_SOURCE
                
                tool['citation_source'] = source
                citation_sources[source] = citation_sources.get(source, 0) + 1
                changed = True
                tools_needing_attribution += 1
                logger.debug(f"Added citation_source={source} to {tool.get('name')}")
            else:
                # Count existing sources
                source = tool['citation_source']
                citation_sources[source] = citation_sources.get(source, 0) + 1
            
            # Add yearly citation source attribution
            if 'yearly_citation_source' not in tool and tool.get('citations_by_year'):
                # Try to get source from existing fields
                source = (
                    tool.get('yearly_data_source') or 
                    tool.get('citation_source') or 
                    DEFAULT_SOURCE
                )
                
                tool['yearly_citation_source'] = source
                yearly_citation_sources[source] = yearly_citation_sources.get(source, 0) + 1
                changed = True
                logger.debug(f"Added yearly_citation_source={source} to {tool.get('name')}")
            elif tool.get('citations_by_year'):
                # Count existing sources
                source = tool['yearly_citation_source']
                yearly_citation_sources[source] = yearly_citation_sources.get(source, 0) + 1
        
        logger.info(f"Found {tools_with_citations} tools with citations")
        logger.info(f"Added attribution to {tools_needing_attribution} tools")
        
        # Add citation sources summary if needed
        if (changed or 'citation_sources' not in impact_data) and (citation_sources or yearly_citation_sources):
            impact_data['citation_sources'] = {
                'total_citations': citation_sources,
                'yearly_citations': yearly_citation_sources
            }
            changed = True
            logger.info(f"Added citation_sources summary: {citation_sources}")
            logger.info(f"Added yearly_citation_sources summary: {yearly_citation_sources}")
        
        # Save changes if needed
        if changed:
            with open(impact_data_path, 'w') as f:
                json.dump(impact_data, f, indent=2)
            logger.info(f"Updated {impact_data_path} with attribution data")
            return True
        else:
            logger.info(f"No changes needed for {impact_data_path}")
            return False
    
    except Exception as e:
        logger.error(f"Error adding attribution to {impact_data_path}: {e}")
        return False


def main():
    """Add attribution to impact_data.json."""
    return 0 if add_attribution_to_impact_data() else 1


if __name__ == "__main__":
    import sys
    sys.exit(main())