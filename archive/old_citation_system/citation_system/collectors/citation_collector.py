#!/usr/bin/env python3
"""
Enhanced citation collector for retrieving and aggregating citation data.
"""

import datetime
import json
import logging
import os
import sys
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple, Union

from ..api.citation_registry import get_citation_source, get_prioritized_sources
from ..config import (
    BATCH_SIZE, DATA_JSON_PATH, IMPACT_DATA_PATH, METADATA_DIR, PARALLEL_REQUESTS
)
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
        logger.info("================================================================================")
        logger.info("========================= Collecting Citation Data =============================")
        logger.info("================================================================================")
        
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
        
        # Add source attribution information
        if 'citation_source' in citation_data:
            tool_entry['citation_source'] = citation_data['citation_source']
        if 'yearly_citation_source' in citation_data:
            tool_entry['yearly_citation_source'] = citation_data['yearly_citation_source']
        
        # Maintain backward compatibility with existing source fields
        if 'primary_source' in citation_data:
            tool_entry['primary_source'] = citation_data['primary_source']
        if 'yearly_data_source' in citation_data:
            tool_entry['yearly_data_source'] = citation_data['yearly_data_source']
        
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
        logger.info("================================================================================")
        logger.info("========================= Generating Impact Data ===============================")
        logger.info("================================================================================")
        
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
        
        # Calculate citation source statistics
        citation_sources = {}
        yearly_citation_sources = {}
        
        for tool in tool_citation_data:
            citation_source = tool.get('citation_source')
            if citation_source:
                citation_sources[citation_source] = citation_sources.get(citation_source, 0) + 1
                
            yearly_source = tool.get('yearly_citation_source')
            if yearly_source:
                yearly_citation_sources[yearly_source] = yearly_citation_sources.get(yearly_source, 0) + 1
        
        # Create impact data structure
        impact_data = {
            'last_updated': self.collection_timestamp,
            'total_tools': total_tools,
            'tools_with_citations': tools_with_citations,
            'total_citations': total_citations,
            'average_citations': average_citations,
            'citation_sources': {
                'total_citations': citation_sources,
                'yearly_citations': yearly_citation_sources
            },
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


def collect_citations(use_cache: bool = True, test_mode: bool = False, limit: Optional[int] = None) -> bool:
    """
    Run citation collection process with specified parameters.
    
    Args:
        use_cache: Whether to use cached responses
        test_mode: Whether to run in test mode
        limit: Maximum number of DOIs to process
    
    Returns:
        True if collection was successful, False otherwise
    """
    try:
        # Create collector
        collector = CitationCollector()
        
        # Step 1: Collect tool DOIs
        tool_dois = collector.collect_tool_dois()
        
        # If in test mode, use a limited subset
        if test_mode or limit:
            # Sort DOIs by value for consistent results
            sorted_items = sorted(tool_dois.items(), key=lambda x: x[1])
            
            # Use first 'limit' items or 10 if in test mode
            max_items = limit or 10
            limited_dois = dict(sorted_items[:max_items])
            
            logger.info(f"Test mode: Using {len(limited_dois)}/{len(tool_dois)} tools")
            tool_dois = limited_dois
        
        # Step 2: Collect citation data
        citation_data = collector.collect_citation_data(
            tool_dois,
            use_cache=use_cache
        )
        
        # Step 3: Generate impact data
        collector.generate_impact_data(citation_data)
        
        return True
    
    except Exception as e:
        logger.exception(f"Error during citation collection: {e}")
        return False


def main():
    """Main entry point for the citation collector."""
    import argparse
    
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
        "--test",
        action="store_true",
        help="Run in test mode (limited subset)"
    )
    parser.add_argument(
        "--limit",
        type=int,
        help="Maximum number of DOIs to process"
    )
    parser.add_argument(
        "--output",
        default=str(IMPACT_DATA_PATH),
        help="Output file path (defaults to impact_data.json)"
    )
    
    args = parser.parse_args()
    
    # Configure logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(Path(IMPACT_DATA_PATH).parent / "citation_collector.log"),
            logging.StreamHandler()
        ]
    )
    
    # Run collection
    success = collect_citations(
        use_cache=not args.no_cache,
        test_mode=args.test,
        limit=args.limit
    )
    
    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())