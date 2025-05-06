#!/usr/bin/env python3
"""
Enhanced citation collector for retrieving and aggregating citation data.
"""

import datetime
import json
import logging
import os
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple, Union

from ..api.citation_registry import get_citation_source, get_prioritized_sources
from ..config import (
    BATCH_SIZE, DATA_JSON_PATH, IMPACT_DATA_PATH, METADATA_DIR, PARALLEL_REQUESTS
)
from ..utils.logging_utils import log_section, log_summary, setup_logging
from ..validators.doi_validator import DOIValidator
from .citation_aggregator import CitationAggregator

logger = logging.getLogger(__name__)


class CitationCollector:
    """Collects and aggregates citation data from various sources."""
    
    def __init__(self):
        """Initialize the citation collector."""
        self.doi_validator = DOIValidator()
        self.citation_aggregator = CitationAggregator()
        self.collection_timestamp = datetime.datetime.now().isoformat()
    
    def collect_tool_dois(self) -> Dict[str, str]:
        """
        Collect DOIs for all tools in the repository.
        
        Returns:
            Dictionary mapping tool IDs to DOIs
        """
        logger.info("Collecting DOIs for tools")
        
        # Dictionary to store tool ID to DOI mapping
        tool_dois = {}
        
        # First, check data.json for DOIs
        try:
            if not DATA_JSON_PATH.exists():
                logger.warning(f"data.json not found at {DATA_JSON_PATH}")
            else:
                with open(DATA_JSON_PATH, 'r') as f:
                    data = json.load(f)
                
                # Check the structure of data.json
                if 'tools' in data:
                    tools = data['tools']
                elif 'nodes' in data:
                    # Filter for nodes of type 'tool'
                    tools = [node for node in data.get('nodes', []) 
                            if node.get('type') == 'tool']
                else:
                    tools = []
                
                for tool in tools:
                    tool_id = tool.get('id')
                    doi = tool.get('doi')
                    
                    if tool_id and doi:
                        # Normalize and validate DOI
                        normalized_doi = self.doi_validator.normalize_doi(doi)
                        if normalized_doi:
                            tool_dois[tool_id] = normalized_doi
        
        except Exception as e:
            logger.error(f"Error reading data.json: {e}")
        
        # Next, check individual metadata files for DOIs
        try:
            for metadata_dir in [METADATA_DIR, METADATA_DIR / "academic_impact", METADATA_DIR / "bioinformatics"]:
                if not metadata_dir.exists():
                    continue
                
                for file_path in metadata_dir.glob("*.json"):
                    try:
                        tool_id = f"tool-{file_path.stem}"
                        
                        with open(file_path, 'r') as f:
                            metadata = json.load(f)
                        
                        doi = metadata.get('doi')
                        if doi:
                            normalized_doi = self.doi_validator.normalize_doi(doi)
                            if normalized_doi:
                                tool_dois[tool_id] = normalized_doi
                    
                    except Exception as e:
                        logger.warning(f"Error reading metadata file {file_path}: {e}")
        
        except Exception as e:
            logger.error(f"Error reading metadata files: {e}")
        
        logger.info(f"Found {len(tool_dois)} tools with DOIs")
        return tool_dois
    
    def collect_citation_data(
        self,
        dois: Dict[str, str],
        use_cache: bool = True,
        force_refresh: bool = False
    ) -> Dict[str, Dict[str, Any]]:
        """
        Collect citation data for DOIs.
        
        Args:
            dois: Dictionary mapping tool IDs to DOIs
            use_cache: Whether to use cached responses
            force_refresh: Whether to force refresh all citation data
        
        Returns:
            Dictionary mapping DOIs to citation data
        """
        log_section(logger, "Collecting Citation Data")
        
        # Deduplicate DOIs (multiple tools can have the same DOI)
        unique_dois = list(set(dois.values()))
        total_dois = len(unique_dois)
        
        logger.info(f"Collecting citation data for {total_dois} unique DOIs")
        
        # Split DOIs into batches
        doi_batches = [
            unique_dois[i:i + BATCH_SIZE]
            for i in range(0, total_dois, BATCH_SIZE)
        ]
        
        # Initialize results dictionary
        citation_data = {}
        
        # Process batches of DOIs
        for batch_index, doi_batch in enumerate(doi_batches):
            batch_start_time = time.time()
            logger.info(f"Processing batch {batch_index + 1}/{len(doi_batches)} ({len(doi_batch)} DOIs)")
            
            # Process DOIs in parallel
            with ThreadPoolExecutor(max_workers=PARALLEL_REQUESTS) as executor:
                # Submit tasks for each DOI
                future_to_doi = {
                    executor.submit(
                        self.citation_aggregator.collect_citations, 
                        doi, 
                        use_cache and not force_refresh
                    ): doi
                    for doi in doi_batch
                }
                
                # Process completed tasks
                for future in as_completed(future_to_doi):
                    doi = future_to_doi[future]
                    try:
                        data = future.result()
                        citation_data[doi] = data
                    except Exception as e:
                        logger.error(f"Error processing DOI {doi}: {e}")
                        citation_data[doi] = {
                            "doi": doi,
                            "error": str(e),
                            "timestamp": self.collection_timestamp
                        }
            
            batch_duration = time.time() - batch_start_time
            logger.info(f"Batch {batch_index + 1} completed in {batch_duration:.2f}s")
        
        # Citation data integrity check
        citation_counts = [data.get('total_citations', 0) for data in citation_data.values() if "error" not in data]
        nonzero_citations = [c for c in citation_counts if c > 0]
        logger.info(f"Citation integrity check: {len(nonzero_citations)}/{len(citation_counts)} DOIs have non-zero citations")
        if nonzero_citations:
            logger.info(f"Sample citation counts: {nonzero_citations[:5]}")
            logger.info(f"Total citations summed from all DOIs: {sum(citation_counts)}")
            
            # Identify top cited DOIs for verification
            top_cited = []
            for doi, data in citation_data.items():
                if "error" not in data and data.get('total_citations', 0) > 0:
                    top_cited.append((doi, data.get('total_citations', 0)))
            top_cited.sort(key=lambda x: x[1], reverse=True)
            logger.info(f"Top 5 cited DOIs: {top_cited[:5]}")
        
        # Log summary
        successful_dois = sum(1 for d in citation_data.values() if "error" not in d)
        logger.info(f"Citation data collection complete: {successful_dois}/{total_dois} DOIs processed successfully")
        
        return citation_data
    
    def _get_tool_metadata(self) -> Dict[str, Dict[str, Any]]:
        """
        Get metadata for all tools.
        
        Returns:
            Dictionary mapping tool IDs to metadata
        """
        tool_metadata = {}
        
        # Try to load from data.json first
        try:
            if DATA_JSON_PATH.exists():
                with open(DATA_JSON_PATH, 'r') as f:
                    data = json.load(f)
                
                # Check the structure of data.json
                tools = []
                if 'tools' in data:
                    tools = data['tools']
                elif 'nodes' in data:
                    # Filter for nodes of type 'tool'
                    tools = [node for node in data.get('nodes', []) 
                            if node.get('type') == 'tool']
                
                for tool in tools:
                    tool_id = tool.get('id')
                    if tool_id:
                        tool_metadata[tool_id] = {
                            'name': tool.get('name', tool_id),
                            'url': tool.get('url', ''),
                            'category': tool.get('category', 'Other Tools')
                        }
        
        except Exception as e:
            logger.error(f"Error reading data.json: {e}")
        
        # Check individual metadata files
        try:
            for metadata_dir in [METADATA_DIR, METADATA_DIR / "academic_impact", METADATA_DIR / "bioinformatics"]:
                if not metadata_dir.exists():
                    continue
                
                for file_path in metadata_dir.glob("*.json"):
                    try:
                        tool_id = f"tool-{file_path.stem}"
                        
                        with open(file_path, 'r') as f:
                            metadata = json.load(f)
                        
                        tool_metadata[tool_id] = {
                            'name': metadata.get('name', file_path.stem),
                            'url': metadata.get('url', ''),
                            'category': metadata.get('category', 'Other Tools')
                        }
                    
                    except Exception as e:
                        logger.warning(f"Error reading metadata file {file_path}: {e}")
        
        except Exception as e:
            logger.error(f"Error reading metadata files: {e}")
        
        return tool_metadata
    
    def _create_tool_citation_data(
        self,
        tool_id: str,
        tool_name: str,
        tool_url: str,
        citation_data: Dict[str, Any],
        category: str = "Other Tools"
    ) -> Dict[str, Any]:
        """
        Create tool citation data from citation data.
        
        Args:
            tool_id: Tool ID
            tool_name: Tool name
            tool_url: Tool URL
            citation_data: Citation data for the tool's DOI
            category: Tool category
        
        Returns:
            Tool citation data ready for impact_data.json
        """
        # Add debug logging for citation processing
        if citation_data.get('total_citations', 0) > 0:
            logger.debug(f"Processing tool {tool_name} with {citation_data.get('total_citations', 0)} citations, source: {citation_data.get('source', 'unknown')}")
        
        # Extract citation information
        total_citations = citation_data.get('total_citations', 0)
        
        # More detailed logging for non-zero citation counts
        if total_citations > 0:
            logger.info(f"Tool {tool_name} has {total_citations} citations from {citation_data.get('source', 'unknown')}")
        
        # Extract influential citations if available
        influential_citations = citation_data.get('influential_citations', 0)
        
        # Extract citations by year
        citations_by_year = citation_data.get('citations_by_year', {})
        
        # Create the tool citation data
        tool_entry = {
            'name': tool_name,
            'url': tool_url,
            'doi': citation_data.get('doi', ''),
            'total_citations': total_citations,
            'influential_citations': influential_citations,
            'category': category or 'Other Tools',
            'citations_by_year': citations_by_year
        }
        
        # Add specialized metrics if available
        if 'rcr' in citation_data:
            tool_entry['rcr'] = citation_data['rcr']
        if 'expected_citations' in citation_data:
            tool_entry['expected_citations'] = citation_data['expected_citations']
        if 'field_citation_rate' in citation_data:
            tool_entry['field_citation_rate'] = citation_data['field_citation_rate']
        
        return tool_entry
    
    def generate_impact_data(
        self, 
        citation_data: Dict[str, Dict[str, Any]],
        write_to_file: bool = True
    ) -> Dict[str, Any]:
        """
        Generate impact_data.json file from collected citation data.
        
        Args:
            citation_data: Dictionary mapping DOIs to citation data
            write_to_file: Whether to write the impact data to a file
        
        Returns:
            Generated impact data
        """
        log_section(logger, "Generating Impact Data")
        
        # Collect tool metadata
        tool_metadata = self._get_tool_metadata()
        
        # Map DOIs to tool IDs
        tool_dois = self.collect_tool_dois()
        doi_to_tools = defaultdict(list)
        
        for tool_id, doi in tool_dois.items():
            doi_to_tools[doi].append(tool_id)
        
        # Generate tool citation data
        tool_citation_data = []
        
        for doi, citation in citation_data.items():
            # Skip entries with errors
            if "error" in citation:
                continue
            
            # Get tools with this DOI
            tool_ids = doi_to_tools.get(doi, [])
            
            for tool_id in tool_ids:
                # Get tool metadata
                metadata = tool_metadata.get(tool_id, {})
                
                tool_name = metadata.get('name', tool_id)
                tool_url = metadata.get('url', '')
                category = metadata.get('category', 'Other Tools')
                
                # Create the tool citation entry
                tool_entry = self._create_tool_citation_data(
                    tool_id, 
                    tool_name,
                    tool_url,
                    citation,
                    category
                )
                
                tool_citation_data.append(tool_entry)
        
        # Sort tools by total citations (descending)
        tool_citation_data.sort(key=lambda x: x['total_citations'], reverse=True)
        
        # Enhanced logging for citation data
        logger.info(f"Raw citation count sample: {[tool['total_citations'] for tool in tool_citation_data[:10]]}")
        logger.info(f"First 5 tools with citations: {[{'name': t['name'], 'citations': t['total_citations']} for t in tool_citation_data[:5]]}")
        
        # Calculate summary statistics
        total_tools = len(tool_metadata)
        tools_with_citations = len(tool_citation_data)
        total_citations = sum(tool['total_citations'] for tool in tool_citation_data)
        logger.info(f"Citation sum calculation: sum of {len(tool_citation_data)} tools = {total_citations}")
        average_citations = total_citations / tools_with_citations if tools_with_citations > 0 else 0
        
        # Create impact data structure
        impact_data = {
            'last_updated': self.collection_timestamp,
            'total_tools': total_tools,
            'tools_with_citations': tools_with_citations,
            'total_citations': total_citations,
            'average_citations': average_citations,
            'tools': tool_citation_data
        }
        
        # Log summary
        logger.info(f"Generated impact data with {tools_with_citations}/{total_tools} tools having citation data")
        logger.info(f"Total citations: {total_citations}, Average: {average_citations:.2f}")
        
        # Write to file if requested
        if write_to_file:
            try:
                with open(IMPACT_DATA_PATH, 'w') as f:
                    json.dump(impact_data, f, indent=2)
                logger.info(f"Impact data written to {IMPACT_DATA_PATH}")
                
                # Verify the file was written correctly
                try:
                    with open(IMPACT_DATA_PATH, 'r') as f:
                        written_data = json.load(f)
                    
                    # Verify total_citations was written correctly
                    if written_data.get('total_citations', 0) != total_citations:
                        logger.error(f"ERROR: File validation failed - calculated {total_citations} citations but file contains {written_data.get('total_citations', 0)}")
                        # Fix the file by writing it again
                        with open(IMPACT_DATA_PATH, 'w') as f:
                            json.dump(impact_data, f, indent=2)
                        logger.info("File was corrected after validation failure")
                    else:
                        logger.info(f"File validation successful - {total_citations} citations written correctly")
                except Exception as e:
                    logger.error(f"Error during file validation: {e}")
            except Exception as e:
                logger.error(f"Error writing impact data to file: {e}")
        
        return impact_data
    
    def run_full_collection(
        self,
        use_cache: bool = True,
        force_refresh: bool = False
    ) -> Dict[str, Any]:
        """
        Run the full citation collection process.
        
        Args:
            use_cache: Whether to use cached responses
            force_refresh: Whether to force refresh all citation data
        
        Returns:
            Generated impact data
        """
        # Step 1: Collect tool DOIs
        tool_dois = self.collect_tool_dois()
        
        # Step 2: Collect citation data
        citation_data = self.collect_citation_data(
            tool_dois,
            use_cache=use_cache,
            force_refresh=force_refresh
        )
        
        # Step 3: Generate impact data
        impact_data = self.generate_impact_data(citation_data)
        
        return impact_data


def main():
    """Main entry point for the citation collector."""
    import argparse
    import sys
    
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Citation Collector")
    parser.add_argument(
        "--no-cache",
        action="store_true",
        help="Don't use cached API responses"
    )
    parser.add_argument(
        "--force-refresh",
        action="store_true",
        help="Force refresh all citation data"
    )
    parser.add_argument(
        "--output",
        default=str(IMPACT_DATA_PATH),
        help="Output file path (defaults to impact_data.json)"
    )
    
    args = parser.parse_args()
    
    # Set up logging
    logger = setup_logging("citation_collector")
    
    # Create collector and run
    collector = CitationCollector()
    
    try:
        impact_data = collector.run_full_collection(
            use_cache=not args.no_cache,
            force_refresh=args.force_refresh
        )
        
        # Log summary
        logger.info(f"Citation collection complete: {impact_data['tools_with_citations']} tools with citation data")
        return 0
    
    except Exception as e:
        logger.error(f"Citation collection failed: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())