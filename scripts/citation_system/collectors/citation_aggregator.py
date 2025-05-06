#!/usr/bin/env python3
"""
Citation aggregator for combining data from multiple sources.
"""

import datetime
import logging
from typing import Any, Dict, List, Optional, Set, Tuple

from ..api.citation_registry import get_citation_source, get_prioritized_sources
from ..config import CITATION_PRIORITY

logger = logging.getLogger(__name__)


class CitationAggregator:
    """Aggregates citation data from multiple sources."""
    
    def __init__(self):
        """Initialize the citation aggregator."""
        self.timestamp = datetime.datetime.now().isoformat()
    
    def collect_citations(self, doi: str, use_cache: bool = True) -> Dict[str, Any]:
        """
        Collect citation data from all available sources for a DOI.
        
        Args:
            doi: DOI to look up
            use_cache: Whether to use cached data
        
        Returns:
            Aggregated citation data
        """
        # Get available sources in priority order
        source_names = get_prioritized_sources()
        if not source_names:
            logger.warning("No citation sources available")
            return {
                "doi": doi,
                "error": "No citation sources available",
                "timestamp": self.timestamp
            }
        
        # Collect data from each source
        source_data = {}
        for source_name in source_names:
            source = get_citation_source(source_name)
            if source:
                try:
                    if hasattr(source, "get_citation_data"):
                        success, data = source.get_citation_data(doi, use_cache)
                        if success:
                            source_data[source_name] = data
                    elif hasattr(source, "get_citation_count"):  # For backward compatibility
                        success, data = source.get_citation_count(doi, use_cache)
                        if success:
                            source_data[source_name] = data
                except Exception as e:
                    logger.error(f"Error getting citation data from {source_name} for {doi}: {e}")
        
        # Aggregate the data
        return self._aggregate_citation_data(source_data, doi)
    
    def _aggregate_citation_data(self, source_data: Dict[str, Dict[str, Any]], doi: str) -> Dict[str, Any]:
        """
        Aggregate citation data from multiple sources.
        
        Args:
            source_data: Dictionary mapping source names to citation data
            doi: DOI for the citation data
        
        Returns:
            Aggregated citation data
        """
        if not source_data:
            return {
                "doi": doi,
                "error": "No citation data available",
                "timestamp": self.timestamp
            }
        
        # Start with a blank aggregated data structure
        aggregated = {
            "doi": doi,
            "source": "aggregated",
            "sources_used": list(source_data.keys()),
            "total_citations": 0,
            "citations_by_year": {},
            "timestamp": self.timestamp,
            "metadata": {}
        }
        
        # Find the highest priority source for metadata
        prioritized_sources = sorted(
            source_data.keys(),
            key=lambda s: CITATION_PRIORITY.get(s, 999)
        )
        
        # Use metadata from the highest priority source
        if prioritized_sources:
            highest_priority = prioritized_sources[0]
            if "metadata" in source_data[highest_priority]:
                aggregated["metadata"] = source_data[highest_priority]["metadata"]
        
        # Determine citation count based on source priority
        for source_name in prioritized_sources:
            data = source_data[source_name]
            if "total_citations" in data and data["total_citations"] > 0:
                aggregated["total_citations"] = data["total_citations"]
                aggregated["primary_source"] = source_name
                
                # Use citations by year if available
                if "citations_by_year" in data and data["citations_by_year"]:
                    aggregated["citations_by_year"] = data["citations_by_year"]
                
                # Special metrics from different sources
                if source_name == "icite":
                    if "rcr" in data:
                        aggregated["rcr"] = data["rcr"]
                    if "expected_citations" in data:
                        aggregated["expected_citations"] = data["expected_citations"]
                    if "field_citation_rate" in data:
                        aggregated["field_citation_rate"] = data["field_citation_rate"]
                
                # Once we have citation data from highest priority source, we can stop
                break
        
        return aggregated
    
    def collect_citations_batch(self, dois: List[str], use_cache: bool = True) -> Dict[str, Dict[str, Any]]:
        """
        Collect citation data for multiple DOIs.
        
        Args:
            dois: List of DOIs to look up
            use_cache: Whether to use cached data
        
        Returns:
            Dictionary mapping DOIs to aggregated citation data
        """
        results = {}
        
        for doi in dois:
            results[doi] = self.collect_citations(doi, use_cache)
        
        return results