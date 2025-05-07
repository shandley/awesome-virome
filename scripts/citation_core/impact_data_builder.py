#!/usr/bin/env python3
"""
Impact Data Builder for Awesome-Virome

This module converts collected citation data into the final impact_data.json format
used by the Awesome-Virome project. It:

1. Loads citation data collected from iCite
2. Merges it with tool metadata from data.json
3. Calculates impact metrics and statistics
4. Generates the impact_data.json file

The builder strictly adheres to the principle of only using real citation data
and never generating synthetic/mock data for visualization purposes.
"""

import os
import json
import logging
from datetime import datetime
from typing import Dict, List, Any, Optional
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler("impact_data_builder.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger("impact_data_builder")

# Constants
DEFAULT_DATA_PATH = "data.json"
DEFAULT_CITATION_DIR = "metadata/citations"
DEFAULT_OUTPUT_PATH = "impact_data.json"
CITATION_RESULTS_FILE = "citation_data.json"
CITATION_SUMMARY_FILE = "citation_summary.json"


class ImpactDataBuilder:
    """Builder for impact_data.json from citation data."""
    
    def __init__(self, data_path: str = DEFAULT_DATA_PATH,
                 citation_dir: str = DEFAULT_CITATION_DIR,
                 output_path: str = DEFAULT_OUTPUT_PATH):
        """
        Initialize the impact data builder.
        
        Args:
            data_path: Path to data.json
            citation_dir: Directory containing citation data
            output_path: Path for output impact_data.json
        """
        self.data_path = data_path
        self.citation_dir = citation_dir
        self.output_path = output_path
        
        logger.info(f"Initialized impact data builder "
                  f"(data_path={data_path}, citation_dir={citation_dir}, output_path={output_path})")
    
    def load_tools(self) -> List[Dict[str, Any]]:
        """
        Load tool data from data.json.
        
        Returns:
            List of tool metadata dictionaries
        """
        try:
            with open(self.data_path, 'r', encoding='utf-8') as f:
                data = json.load(f)
            
            # Handle different data structures
            if "tools" in data:
                tools = data["tools"]
            elif "nodes" in data:
                # Filter for nodes of type 'tool'
                tools = [node for node in data.get("nodes", []) if node.get("type") == "tool"]
            else:
                tools = []
            
            logger.info(f"Loaded {len(tools)} tools from {self.data_path}")
            return tools
            
        except (json.JSONDecodeError, FileNotFoundError) as e:
            logger.error(f"Error loading {self.data_path}: {e}")
            return []
    
    def load_citation_data(self) -> Dict[str, Dict[str, Any]]:
        """
        Load citation data from citation_data.json.
        
        Returns:
            Dictionary mapping tool names to citation data
        """
        citation_path = os.path.join(self.citation_dir, CITATION_RESULTS_FILE)
        
        try:
            with open(citation_path, 'r', encoding='utf-8') as f:
                data = json.load(f)
            
            logger.info(f"Loaded citation data for {len(data)} tools from {citation_path}")
            return data
            
        except (json.JSONDecodeError, FileNotFoundError) as e:
            logger.error(f"Error loading {citation_path}: {e}")
            return {}
    
    def load_citation_summary(self) -> Dict[str, Any]:
        """
        Load citation summary data.
        
        Returns:
            Dictionary containing citation summary
        """
        summary_path = os.path.join(self.citation_dir, CITATION_SUMMARY_FILE)
        
        try:
            with open(summary_path, 'r', encoding='utf-8') as f:
                data = json.load(f)
            
            logger.info(f"Loaded citation summary from {summary_path}")
            return data
            
        except (json.JSONDecodeError, FileNotFoundError) as e:
            logger.error(f"Error loading {summary_path}: {e}")
            return {}
    
    def build_impact_data(self) -> Dict[str, Any]:
        """
        Build the impact_data.json structure.
        
        Returns:
            Dictionary containing impact data
        """
        # Load tools and citation data
        tools = self.load_tools()
        citation_data = self.load_citation_data()
        citation_summary = self.load_citation_summary()
        
        if not tools:
            logger.error("No tools found, cannot build impact data")
            return {}
        
        # Create tool name to metadata lookup
        tool_lookup = {tool.get("name", ""): tool for tool in tools if tool.get("name")}
        
        # Statistics
        total_tools = len(tools)
        tools_with_citations = citation_summary.get("tools_with_citations", 0)
        total_citations = citation_summary.get("total_citations", 0)
        
        # Calculate average citations
        average_citations = 0
        if tools_with_citations > 0:
            average_citations = round(total_citations / tools_with_citations, 2)
        
        # Build impact data structure
        impact_tools = []
        
        # First add tools with citation data
        for tool_name, data in citation_data.items():
            # Get tool metadata
            tool_meta = tool_lookup.get(tool_name, {})
            
            # Create impact tool entry
            impact_tool = {
                "name": tool_name,
                "url": data.get("url") or tool_meta.get("url", ""),
                "doi": data.get("doi", ""),
                "total_citations": data.get("total_citations", 0),
                "citation_source": "icite",
                "citations_by_year": data.get("citations_by_year", {})
            }
            
            # Add category if available
            if "category" in tool_meta:
                impact_tool["category"] = tool_meta["category"]
            elif "category" in data:
                impact_tool["category"] = data["category"]
            else:
                impact_tool["category"] = "Other Tools"
            
            # Add publication metadata if available
            if "publication" in data:
                impact_tool["publication"] = data["publication"]
            
            impact_tools.append(impact_tool)
        
        # Then add tools without citation data with zero citations
        tools_with_data = set(citation_data.keys())
        
        for tool in tools:
            tool_name = tool.get("name", "")
            
            if not tool_name or tool_name in tools_with_data:
                continue
            
            # Create impact tool entry with zero citations
            impact_tool = {
                "name": tool_name,
                "url": tool.get("url", ""),
                "doi": tool.get("doi", ""),
                "total_citations": 0,
                "citation_source": "none",
                "citations_by_year": {}
            }
            
            # Add category if available
            if "category" in tool:
                impact_tool["category"] = tool["category"]
            else:
                impact_tool["category"] = "Other Tools"
            
            impact_tools.append(impact_tool)
        
        # Sort tools by total citations (descending)
        impact_tools.sort(key=lambda x: x.get("total_citations", 0), reverse=True)
        
        # Create impact data structure
        impact_data = {
            "last_updated": datetime.now().isoformat(),
            "tools": impact_tools,
            "total_tools": total_tools,
            "tools_with_citations": tools_with_citations,
            "total_citations": total_citations,
            "average_citations": average_citations,
            "citation_sources": {
                "primary": "icite",
                "total_citations": {
                    "icite": tools_with_citations,
                    "none": total_tools - tools_with_citations
                }
            }
        }
        
        # Add citations by year if available
        if "citations_by_year" in citation_summary:
            impact_data["citations_by_year"] = citation_summary["citations_by_year"]
        
        return impact_data
    
    def save_impact_data(self, impact_data: Dict[str, Any]) -> bool:
        """
        Save impact data to impact_data.json.
        
        Args:
            impact_data: Dictionary containing impact data
            
        Returns:
            True if saved successfully, False otherwise
        """
        try:
            with open(self.output_path, 'w', encoding='utf-8') as f:
                json.dump(impact_data, f, indent=2)
            
            logger.info(f"Saved impact data to {self.output_path}")
            return True
            
        except IOError as e:
            logger.error(f"Error saving impact data: {e}")
            return False
    
    def build_and_save(self) -> bool:
        """
        Build and save impact data.
        
        Returns:
            True if successful, False otherwise
        """
        impact_data = self.build_impact_data()
        
        if not impact_data:
            logger.error("Failed to build impact data")
            return False
        
        # Log statistics
        total_tools = impact_data.get("total_tools", 0)
        tools_with_citations = impact_data.get("tools_with_citations", 0)
        total_citations = impact_data.get("total_citations", 0)
        
        logger.info(f"Built impact data with:")
        logger.info(f"- {total_tools} total tools")
        logger.info(f"- {tools_with_citations} tools with citations ({(tools_with_citations/total_tools*100):.1f}%)")
        logger.info(f"- {total_citations} total citations")
        
        # Save impact data
        return self.save_impact_data(impact_data)


def main():
    """Main function to run the impact data builder."""
    import argparse
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Build impact_data.json from citation data")
    parser.add_argument("--data", default=DEFAULT_DATA_PATH, help="Path to data.json")
    parser.add_argument("--citations", default=DEFAULT_CITATION_DIR, help="Directory containing citation data")
    parser.add_argument("--output", default=DEFAULT_OUTPUT_PATH, help="Path for output impact_data.json")
    args = parser.parse_args()
    
    # Run builder
    builder = ImpactDataBuilder(
        data_path=args.data,
        citation_dir=args.citations,
        output_path=args.output
    )
    
    if builder.build_and_save():
        logger.info("Successfully built and saved impact data")
        return 0
    else:
        logger.error("Failed to build and save impact data")
        return 1


if __name__ == "__main__":
    import sys
    sys.exit(main())