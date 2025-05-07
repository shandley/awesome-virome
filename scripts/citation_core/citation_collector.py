#!/usr/bin/env python3
"""
Citation Collector for Awesome-Virome

This module implements the main citation collection process for the Awesome-Virome
project. It:

1. Loads tool data from data.json
2. Extracts DOIs from tool metadata
3. Fetches citation data from iCite API
4. Processes and validates the citation data
5. Stores the results for integration into impact_data.json

The collector only uses real citation data from iCite and never generates
synthetic/mock data for visualization purposes.
"""

import os
import json
import logging
import time
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple

# Import local modules
from icite_api import ICiteClient
from doi_handler import extract_dois_from_tools, is_valid_doi

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler("citation_collection.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger("citation_collector")

# Constants
DEFAULT_DATA_PATH = "data.json"
DEFAULT_OUTPUT_DIR = "metadata/citations"
CITATION_RESULTS_FILE = "citation_data.json"
CITATION_SUMMARY_FILE = "citation_summary.json"


class CitationCollector:
    """Collects citation data for Awesome-Virome tools."""
    
    def __init__(self, data_path: str = DEFAULT_DATA_PATH, 
                 output_dir: str = DEFAULT_OUTPUT_DIR,
                 icite_rate_limit: float = 3.0):
        """
        Initialize the citation collector.
        
        Args:
            data_path: Path to the data.json file
            output_dir: Directory to store citation data
            icite_rate_limit: Rate limit for iCite API requests (seconds)
        """
        self.data_path = data_path
        self.output_dir = output_dir
        
        # Ensure output directory exists
        os.makedirs(output_dir, exist_ok=True)
        
        # Initialize API client
        self.icite_client = ICiteClient(rate_limit=icite_rate_limit)
        
        # Statistics for reporting
        self.stats = {
            "tools_processed": 0,
            "tools_with_doi": 0,
            "tools_with_citations": 0,
            "total_citations": 0,
            "start_time": datetime.now().isoformat()
        }
        
        logger.info(f"Initialized citation collector (data_path={data_path}, output_dir={output_dir})")
    
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
    
    def collect_citations(self) -> Dict[str, Dict[str, Any]]:
        """
        Collect citation data for all tools.
        
        Returns:
            Dictionary mapping tool names to their citation data
        """
        # Load tools
        tools = self.load_tools()
        if not tools:
            logger.error("No tools found, aborting collection")
            return {}
        
        # Extract DOIs
        tool_dois = extract_dois_from_tools(tools)
        self.stats["tools_with_doi"] = len(tool_dois)
        
        if not tool_dois:
            logger.warning("No DOIs found for any tools, aborting collection")
            return {}
        
        logger.info(f"Found DOIs for {len(tool_dois)} tools")
        
        # Create lookup from tool name to tool data
        tool_lookup = {tool.get("name", ""): tool for tool in tools if tool.get("name")}
        
        # Collect citation data
        results = {}
        
        for tool_name, doi in tool_dois.items():
            try:
                self.stats["tools_processed"] += 1
                
                logger.info(f"Processing tool {self.stats['tools_processed']}/{len(tool_dois)}: {tool_name}")
                
                # Skip invalid DOIs
                if not is_valid_doi(doi):
                    logger.warning(f"Invalid DOI for {tool_name}: {doi}")
                    continue
                
                # Get citation data from iCite
                citation_data = self.icite_client.get_citation_data_by_doi(doi)
                
                if not citation_data:
                    logger.warning(f"No citation data found for {tool_name} (DOI: {doi})")
                    continue
                
                # Extract tool metadata
                tool_data = tool_lookup.get(tool_name, {})
                
                # Create citation result
                result = {
                    "name": tool_name,
                    "url": tool_data.get("url", ""),
                    "doi": doi,
                    "total_citations": citation_data.get("citation_count", 0),
                    "citation_source": "icite",
                    "citations_by_year": citation_data.get("citations_by_year", {}),
                    "publication": {
                        "title": citation_data.get("title", ""),
                        "authors": citation_data.get("authors", ""),
                        "journal": citation_data.get("journal", ""),
                        "year": citation_data.get("year", 0),
                        "pmid": citation_data.get("pmid", "")
                    },
                    "category": tool_data.get("category", "Other Tools"),
                    "last_updated": datetime.now().isoformat()
                }
                
                # Update statistics
                citation_count = result["total_citations"]
                if citation_count > 0:
                    self.stats["tools_with_citations"] += 1
                    self.stats["total_citations"] += citation_count
                
                # Store result
                results[tool_name] = result
                
                # Log progress
                if citation_count > 0:
                    logger.info(f"Found {citation_count} citations for {tool_name}")
                else:
                    logger.info(f"No citations found for {tool_name}")
                
            except Exception as e:
                logger.error(f"Error processing {tool_name}: {e}")
                continue
        
        # Save results
        self._save_results(results)
        
        return results
    
    def _save_results(self, results: Dict[str, Dict[str, Any]]) -> bool:
        """
        Save citation results to files.
        
        Args:
            results: Dictionary mapping tool names to citation data
            
        Returns:
            True if saved successfully, False otherwise
        """
        try:
            # Save individual tool data
            for tool_name, data in results.items():
                safe_name = tool_name.replace("/", "_").replace(" ", "_")
                file_path = os.path.join(self.output_dir, f"{safe_name}.json")
                
                with open(file_path, 'w', encoding='utf-8') as f:
                    json.dump(data, f, indent=2)
            
            # Save consolidated results
            results_path = os.path.join(self.output_dir, CITATION_RESULTS_FILE)
            with open(results_path, 'w', encoding='utf-8') as f:
                json.dump(results, f, indent=2)
            
            # Generate and save summary
            self._save_summary(results)
            
            logger.info(f"Saved citation data for {len(results)} tools to {self.output_dir}")
            return True
            
        except IOError as e:
            logger.error(f"Error saving results: {e}")
            return False
    
    def _save_summary(self, results: Dict[str, Dict[str, Any]]) -> bool:
        """
        Generate and save a summary of citation data.
        
        Args:
            results: Dictionary mapping tool names to citation data
            
        Returns:
            True if saved successfully, False otherwise
        """
        try:
            # Update stats
            self.stats["end_time"] = datetime.now().isoformat()
            self.stats["elapsed_seconds"] = (datetime.fromisoformat(self.stats["end_time"]) - 
                                           datetime.fromisoformat(self.stats["start_time"])).total_seconds()
            
            # Citations by year across all tools
            citations_by_year = {}
            
            for tool_data in results.values():
                for year, count in tool_data.get("citations_by_year", {}).items():
                    if year in citations_by_year:
                        citations_by_year[year] += count
                    else:
                        citations_by_year[year] = count
            
            # Sort years
            citations_by_year = dict(sorted(citations_by_year.items()))
            
            # Most cited tools
            tools_by_citations = sorted(
                results.values(), 
                key=lambda x: x.get("total_citations", 0), 
                reverse=True
            )
            
            most_cited = []
            for tool in tools_by_citations[:20]:  # Top 20
                if tool.get("total_citations", 0) > 0:
                    most_cited.append({
                        "name": tool.get("name", ""),
                        "citations": tool.get("total_citations", 0),
                        "doi": tool.get("doi", "")
                    })
            
            # Create summary
            summary = {
                "total_tools": self.stats["tools_processed"],
                "tools_with_doi": self.stats["tools_with_doi"],
                "tools_with_citations": self.stats["tools_with_citations"],
                "total_citations": self.stats["total_citations"],
                "most_cited": most_cited,
                "citations_by_year": citations_by_year,
                "collection_stats": {
                    "start_time": self.stats["start_time"],
                    "end_time": self.stats["end_time"],
                    "elapsed_seconds": self.stats["elapsed_seconds"]
                },
                "citation_source": "icite"
            }
            
            # Save summary
            summary_path = os.path.join(self.output_dir, CITATION_SUMMARY_FILE)
            with open(summary_path, 'w', encoding='utf-8') as f:
                json.dump(summary, f, indent=2)
            
            logger.info(f"Saved citation summary to {summary_path}")
            
            # Log summary
            logger.info(f"Citation collection summary:")
            logger.info(f"- Total tools processed: {summary['total_tools']}")
            logger.info(f"- Tools with DOI: {summary['tools_with_doi']}")
            logger.info(f"- Tools with citations: {summary['tools_with_citations']}")
            logger.info(f"- Total citations: {summary['total_citations']}")
            
            return True
            
        except (IOError, ValueError) as e:
            logger.error(f"Error saving summary: {e}")
            return False


def main():
    """Main function to run the citation collector."""
    import argparse
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Collect citation data for Awesome-Virome tools")
    parser.add_argument("--data", default=DEFAULT_DATA_PATH, help="Path to data.json")
    parser.add_argument("--output", default=DEFAULT_OUTPUT_DIR, help="Output directory for citation data")
    parser.add_argument("--rate-limit", type=float, default=3.0, help="Rate limit for iCite API (seconds)")
    args = parser.parse_args()
    
    # Run collector
    collector = CitationCollector(
        data_path=args.data,
        output_dir=args.output,
        icite_rate_limit=args.rate_limit
    )
    
    start_time = time.time()
    results = collector.collect_citations()
    elapsed_time = time.time() - start_time
    
    # Report results
    logger.info(f"Citation collection completed in {elapsed_time:.2f} seconds")
    logger.info(f"Collected citation data for {len(results)} tools")
    
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())