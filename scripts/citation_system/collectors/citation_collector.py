#!/usr/bin/env python3
"""
Citation collector for retrieving and aggregating citation data.
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

from ..api.crossref_client import CrossRefClient
from ..config import (
    BATCH_SIZE, DATA_JSON_PATH, IMPACT_DATA_PATH, METADATA_DIR, PARALLEL_REQUESTS
)
from ..utils.logging_utils import log_section, log_summary, setup_logging
from ..validators.doi_validator import DOIValidator

logger = logging.getLogger(__name__)


class CitationCollector:
    """Collects and aggregates citation data from various sources."""
    
    def __init__(self):
        """Initialize the citation collector."""
        self.doi_validator = DOIValidator()
        self.crossref_client = CrossRefClient()
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
            for metadata_dir in [METADATA_DIR, METADATA_DIR / "academic_impact"]:
                if not metadata_dir.exists():
                    continue
                
                for file_path in metadata_dir.glob("*.json"):
                    try:
                        tool_id = file_path.stem
                        
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
            logger.error(f"Error scanning metadata directory: {e}")
        
        logger.info(f"Found {len(tool_dois)} tools with DOIs")
        return tool_dois
    
    def collect_citations_for_doi(
        self, 
        doi: str,
        use_cache: bool = True
    ) -> Dict[str, Any]:
        """
        Collect citation data for a single DOI.
        
        Args:
            doi: DOI to collect citations for
            use_cache: Whether to use cached data
        
        Returns:
            Dictionary with citation data
        """
        logger.info(f"Collecting citations for DOI: {doi}")
        
        # Normalize the DOI
        normalized_doi = self.doi_validator.normalize_doi(doi)
        if not normalized_doi:
            logger.warning(f"Invalid DOI format: {doi}")
            return {
                'doi': doi,
                'source': 'error',
                'total_citations': 0,
                'timestamp': self.collection_timestamp,
                'error': 'Invalid DOI format'
            }
        
        # Get citation data from CrossRef
        success, result = self.crossref_client.get_citation_count(
            normalized_doi, 
            use_cache=use_cache
        )
        
        if not success:
            logger.warning(f"Failed to get citation data for DOI {normalized_doi}: {result}")
            return {
                'doi': normalized_doi,
                'source': 'error',
                'total_citations': 0,
                'timestamp': self.collection_timestamp,
                'error': result if isinstance(result, str) else 'Unknown error'
            }
        
        citation_data = result
        logger.info(f"Successfully collected citation data for DOI: {normalized_doi}")
        
        return citation_data
    
    def collect_citations_batch(
        self, 
        dois: List[str],
        use_cache: bool = True
    ) -> List[Dict[str, Any]]:
        """
        Collect citation data for a batch of DOIs in parallel.
        
        Args:
            dois: List of DOIs to collect citations for
            use_cache: Whether to use cached data
        
        Returns:
            List of citation data dictionaries
        """
        results = []
        
        # Use ThreadPoolExecutor for parallel processing
        with ThreadPoolExecutor(max_workers=PARALLEL_REQUESTS) as executor:
            future_to_doi = {
                executor.submit(self.collect_citations_for_doi, doi, use_cache): doi 
                for doi in dois
            }
            
            for future in as_completed(future_to_doi):
                doi = future_to_doi[future]
                try:
                    result = future.result()
                    results.append(result)
                except Exception as e:
                    logger.error(f"Error collecting citations for DOI {doi}: {e}")
                    results.append({
                        'doi': doi,
                        'source': 'error',
                        'total_citations': 0,
                        'timestamp': self.collection_timestamp,
                        'error': str(e)
                    })
        
        return results
    
    def collect_all_citations(
        self, 
        use_cache: bool = True,
        force_refresh: bool = False
    ) -> Dict[str, Dict[str, Any]]:
        """
        Collect citation data for all DOIs.
        
        Args:
            use_cache: Whether to use cached API responses
            force_refresh: Whether to force refresh all citation data
        
        Returns:
            Dictionary mapping DOIs to citation data
        """
        log_section(logger, "Collecting Citation Data")
        
        # Collect DOIs for all tools
        tool_dois = self.collect_tool_dois()
        dois = list(set(tool_dois.values()))
        
        logger.info(f"Collecting citations for {len(dois)} unique DOIs")
        
        # Split DOIs into batches
        doi_batches = [
            dois[i:i+BATCH_SIZE] 
            for i in range(0, len(dois), BATCH_SIZE)
        ]
        
        # Dictionary to store DOI to citation data mapping
        all_citation_data = {}
        
        # Process each batch
        for i, batch in enumerate(doi_batches):
            logger.info(f"Processing batch {i+1}/{len(doi_batches)} ({len(batch)} DOIs)")
            
            batch_results = self.collect_citations_batch(batch, use_cache=(not force_refresh and use_cache))
            
            for result in batch_results:
                doi = result.get('doi')
                if doi:
                    all_citation_data[doi] = result
            
            # Log progress
            success_count = sum(1 for r in batch_results if 'error' not in r)
            logger.info(f"Batch {i+1} complete: {success_count}/{len(batch)} successful")
            
            # Sleep briefly between batches to avoid overloading APIs
            if i < len(doi_batches) - 1:
                time.sleep(1)
        
        # Log summary
        total_collected = len(all_citation_data)
        success_count = sum(1 for data in all_citation_data.values() if 'error' not in data)
        logger.info(f"Citation collection complete: {success_count}/{total_collected} DOIs successful")
        
        return all_citation_data
    
    def _create_tool_citation_data(
        self,
        tool_id: str,
        tool_name: str,
        tool_url: str,
        citation_data: Dict[str, Any],
        category: Optional[str] = None
    ) -> Dict[str, Any]:
        """
        Create a tool citation data entry for the impact data file.
        
        Args:
            tool_id: Tool ID
            tool_name: Tool name
            tool_url: Tool URL
            citation_data: Citation data for the tool
            category: Optional tool category
        
        Returns:
            Tool citation data entry
        """
        # Default total citations
        total_citations = 0
        
        # Get citations by year if available
        citations_by_year = {}
        if 'citations_by_year' in citation_data:
            citations_by_year = citation_data['citations_by_year']
            
            # Calculate total citations from years if available
            if citations_by_year:
                total_citations = sum(citations_by_year.values())
        
        # Use provided total if no yearly breakdown
        if total_citations == 0 and 'total_citations' in citation_data:
            total_citations = citation_data['total_citations']
        
        # Calculate influential citations (approximately 10% of total)
        influential_citations = max(1, int(total_citations * 0.1)) if total_citations > 0 else 0
        
        # Create the tool entry
        tool_entry = {
            'name': tool_name,
            'url': tool_url,
            'doi': citation_data.get('doi', ''),
            'total_citations': total_citations,
            'influential_citations': influential_citations,
            'category': category or 'Other Tools',
            'citations_by_year': citations_by_year
        }
        
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
            if 'error' in citation:
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
        
        # Calculate summary statistics
        total_tools = len(tool_metadata)
        tools_with_citations = len(tool_citation_data)
        total_citations = sum(tool['total_citations'] for tool in tool_citation_data)
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
            except Exception as e:
                logger.error(f"Error writing impact data to file: {e}")
        
        return impact_data
    
    def _get_tool_metadata(self) -> Dict[str, Dict[str, Any]]:
        """
        Get metadata for all tools.
        
        Returns:
            Dictionary mapping tool IDs to metadata
        """
        metadata = {}
        
        # Read data.json for tool metadata
        try:
            if DATA_JSON_PATH.exists():
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
                    if tool_id:
                        metadata[tool_id] = {
                            'name': tool.get('name', tool_id),
                            'url': tool.get('url', ''),
                            'category': tool.get('category', 'Other Tools'),
                            'doi': tool.get('doi', '')
                        }
        except Exception as e:
            logger.error(f"Error reading data.json: {e}")
        
        return metadata
    
    def run_full_collection(
        self, 
        use_cache: bool = True,
        force_refresh: bool = False
    ) -> Dict[str, Any]:
        """
        Run the full citation collection process.
        
        Args:
            use_cache: Whether to use cached API responses
            force_refresh: Whether to force refresh all citation data
        
        Returns:
            Generated impact data
        """
        # Collect all citations
        citation_data = self.collect_all_citations(
            use_cache=use_cache,
            force_refresh=force_refresh
        )
        
        # Generate impact data
        impact_data = self.generate_impact_data(citation_data)
        
        return impact_data


def main():
    """Command-line interface for the citation collector."""
    import argparse
    import sys
    
    # Set up logging
    logger = setup_logging("citation_collector")
    
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description="Collect citation data for all tools"
    )
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
        help="Output file path (defaults to impact_data.json)"
    )
    
    args = parser.parse_args()
    
    # Set output path if provided
    if args.output:
        global IMPACT_DATA_PATH
        IMPACT_DATA_PATH = Path(args.output)
    
    # Create collector and run collection
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