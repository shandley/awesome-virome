#!/usr/bin/env python3
"""
Command-line tool for collecting and processing citation data.
"""

import argparse
import json
import logging
import sys
import os
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from scripts.citation_system.collectors.citation_collector import collect_citations
from scripts.citation_system.validators.doi_validator import DOIValidator
from scripts.citation_system.api.citation_registry import get_available_sources, get_prioritized_sources
from scripts.citation_system.api.base_client import BaseAPIClient
from scripts.citation_system.config import LOG_DIR, LOG_FORMAT, DATA_JSON_PATH, IMPACT_DATA_PATH

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format=LOG_FORMAT,
    handlers=[
        logging.FileHandler(LOG_DIR / "citation_system.log"),
        logging.StreamHandler()
    ]
)

logger = logging.getLogger(__name__)


def setup_logger():
    """Set up the logger with detailed formatting."""
    # Create a formatter that includes timestamp, logger name, and level
    formatter = logging.Formatter(LOG_FORMAT)
    
    # Set up file handler
    file_handler = logging.FileHandler(LOG_DIR / "citation_system.log")
    file_handler.setFormatter(formatter)
    
    # Set up console handler
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)
    
    # Configure root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.INFO)
    
    # Remove any existing handlers and add our own
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)
    
    root_logger.addHandler(file_handler)
    root_logger.addHandler(console_handler)


def list_available_sources():
    """List available citation sources."""
    sources = get_available_sources()
    prioritized = get_prioritized_sources()
    
    logger.info("================================================================================")
    logger.info("=========================== Available Citation Sources =========================")
    logger.info("================================================================================")
    
    if not sources:
        logger.info("No citation sources available.")
        return
    
    logger.info(f"Available citation sources: {', '.join(sources)}")
    logger.info(f"Sources in priority order: {', '.join(prioritized)}")
    
    # List each source with more details
    for i, source_name in enumerate(prioritized):
        logger.info(f"{i+1}. {source_name.upper()}")


def test_citation_source(source_name: str, doi: str, use_cache: bool = True):
    """
    Test a specific citation source with a given DOI.
    
    Args:
        source_name: Name of the citation source to test
        doi: DOI to look up
        use_cache: Whether to use cached data
    """
    from scripts.citation_system.api.citation_registry import get_citation_source
    
    logger.info("================================================================================")
    logger.info(f"========================= Testing {source_name.upper()} ==========================")
    logger.info("================================================================================")
    
    # Get the citation source
    source = get_citation_source(source_name)
    if not source:
        logger.error(f"Citation source '{source_name}' not found or not enabled.")
        return
    
    # Validate the DOI
    validator = DOIValidator()
    if not validator.is_valid_format(doi):
        logger.error(f"Invalid DOI format: {doi}")
        return
    
    # Get citation data
    logger.info(f"Getting citation data for DOI {doi} from {source_name}...")
    
    try:
        if hasattr(source, "get_citation_data"):
            success, data = source.get_citation_data(doi, use_cache)
        else:
            success, data = source.get_citation_count(doi, use_cache)
        
        if success:
            # Pretty print the data
            logger.info(f"Successfully retrieved citation data from {source_name}:")
            logger.info(json.dumps(data, indent=2))
            
            # Log citation count
            citation_count = data.get("total_citations", 0)
            logger.info(f"Citation count: {citation_count}")
            
            # Log yearly breakdown if available
            yearly_data = data.get("citations_by_year", {})
            if yearly_data:
                logger.info(f"Yearly citation data available: {len(yearly_data)} years")
                logger.info(f"Years with data: {', '.join(sorted(yearly_data.keys()))}")
            else:
                logger.info("No yearly citation data available")
        else:
            logger.error(f"Failed to get citation data: {data}")
    except Exception as e:
        logger.exception(f"Error testing citation source: {e}")


def validate_dois():
    """Validate DOIs in the repository data."""
    logger.info("================================================================================")
    logger.info("=========================== Validating DOIs ==================================")
    logger.info("================================================================================")
    
    try:
        # Load data.json
        with open(DATA_JSON_PATH, "r") as f:
            data = json.load(f)
        
        # Extract DOIs from tool entries
        dois = []
        for tool in data.get("tools", []):
            tool_dois = tool.get("dois", [])
            if tool_dois:
                dois.extend(tool_dois)
        
        # Validate DOIs
        validator = DOIValidator()
        valid_count = 0
        invalid_count = 0
        invalid_dois = []
        
        for doi in dois:
            if validator.is_valid_format(doi):
                valid_count += 1
            else:
                invalid_count += 1
                invalid_dois.append(doi)
        
        # Log results
        logger.info(f"Total DOIs: {len(dois)}")
        logger.info(f"Valid DOIs: {valid_count}")
        logger.info(f"Invalid DOIs: {invalid_count}")
        
        if invalid_count > 0:
            logger.warning("Invalid DOIs found:")
            for doi in invalid_dois:
                logger.warning(f"  - {doi}")
        else:
            logger.info("All DOIs are valid.")
        
    except Exception as e:
        logger.exception(f"Error validating DOIs: {e}")


def run_collection(no_cache: bool = False, test_mode: bool = False, limit: Optional[int] = None):
    """
    Run citation collection process.
    
    Args:
        no_cache: If True, bypass cache for all API requests
        test_mode: If True, run in test mode (limited subset)
        limit: Maximum number of DOIs to process
    """
    logger.info("================================================================================")
    logger.info("=========================== Collecting Citation Data ===========================")
    logger.info("================================================================================")
    
    try:
        # Get available sources
        sources = get_available_sources()
        if not sources:
            logger.error("No citation sources available. Cannot collect citations.")
            return
        
        # Log collection parameters
        logger.info(f"Collection mode: {'Test mode' if test_mode else 'Full collection'}")
        logger.info(f"Using cache: {not no_cache}")
        if limit is not None:
            logger.info(f"Processing limit: {limit} DOIs")
        logger.info(f"Available sources: {', '.join(sources)}")
        logger.info(f"Sources in priority order: {', '.join(get_prioritized_sources())}")
        
        # Run collection
        success = collect_citations(
            use_cache=not no_cache,
            test_mode=test_mode,
            limit=limit
        )
        
        if success:
            # Verify impact data was created
            if os.path.exists(IMPACT_DATA_PATH):
                # Load and log summary stats
                with open(IMPACT_DATA_PATH, "r") as f:
                    impact_data = json.load(f)
                
                logger.info("Collection completed successfully")
                logger.info(f"Total tools: {impact_data.get('total_tools', 0)}")
                logger.info(f"Tools with citations: {impact_data.get('tools_with_citations', 0)}")
                logger.info(f"Total citations: {impact_data.get('total_citations', 0)}")
                
                # Log data sources used
                if "tools" in impact_data:
                    citation_sources = set()
                    yearly_sources = set()
                    tools_with_yearly = 0
                    
                    for tool in impact_data["tools"]:
                        if "primary_source" in tool:
                            citation_sources.add(tool["primary_source"])
                        if "yearly_data_source" in tool:
                            yearly_sources.add(tool["yearly_data_source"])
                            tools_with_yearly += 1
                    
                    if citation_sources:
                        logger.info(f"Citation sources used: {', '.join(citation_sources)}")
                    if yearly_sources:
                        logger.info(f"Yearly data sources used: {', '.join(yearly_sources)}")
                        logger.info(f"Tools with yearly citation data: {tools_with_yearly}")
            else:
                logger.error("Collection process did not create impact_data.json")
        else:
            logger.error("Collection process failed")
    except Exception as e:
        logger.exception(f"Error running collection: {e}")


def main():
    """Main entry point for the script."""
    parser = argparse.ArgumentParser(description="Collect and process citation data")
    subparsers = parser.add_subparsers(dest="command", help="Command to run")
    
    # Parser for 'sources' command
    sources_parser = subparsers.add_parser("sources", help="List available citation sources")
    
    # Parser for 'test' command
    test_parser = subparsers.add_parser("test", help="Test a specific citation source")
    test_parser.add_argument("source", help="Name of the citation source to test")
    test_parser.add_argument("doi", help="DOI to look up")
    test_parser.add_argument("--no-cache", action="store_true", help="Bypass cache for API requests")
    
    # Parser for 'validate' command
    validate_parser = subparsers.add_parser("validate", help="Validate DOIs in the repository data")
    
    # Parser for 'collect' command
    collect_parser = subparsers.add_parser("collect", help="Run citation collection process")
    collect_parser.add_argument("--no-cache", action="store_true", help="Bypass cache for API requests")
    collect_parser.add_argument("--test", action="store_true", help="Run in test mode (limited subset)")
    collect_parser.add_argument("--limit", type=int, help="Maximum number of DOIs to process")
    
    args = parser.parse_args()
    
    # Set up logging
    setup_logger()
    
    # Handle commands
    if args.command == "sources":
        list_available_sources()
    elif args.command == "test":
        test_citation_source(args.source, args.doi, not args.no_cache)
    elif args.command == "validate":
        validate_dois()
    elif args.command == "collect":
        run_collection(args.no_cache, args.test, args.limit)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()