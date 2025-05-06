#!/usr/bin/env python3
"""
Main entry point for the citation collection system.
"""

import argparse
import json
import logging
import sys
import time
from pathlib import Path
from typing import Dict, List

from .collectors.citation_collector import CitationCollector
from .collectors.doi_scanner import DOIScanner
from .config import IMPACT_DATA_PATH, ROOT_DIR
from .utils.logging_utils import log_section, log_summary, setup_logging
from .validators.doi_validator import DOIValidator

logger = logging.getLogger(__name__)


def validate_dois(dois: List[str], check_resolvable: bool = False) -> Dict[str, int]:
    """
    Validate a list of DOIs.
    
    Args:
        dois: List of DOIs to validate
        check_resolvable: Whether to check if DOIs are resolvable
    
    Returns:
        Summary statistics
    """
    validator = DOIValidator()
    
    total = len(dois)
    valid_format = 0
    resolvable = 0
    
    for doi in dois:
        # Check if the DOI has a valid format
        if validator.is_valid_format(doi):
            valid_format += 1
            
            # Check if the DOI is resolvable
            if check_resolvable and validator.is_resolvable(doi):
                resolvable += 1
    
    logger.info(f"DOI validation: {valid_format}/{total} have valid format")
    
    if check_resolvable:
        logger.info(f"DOI resolution: {resolvable}/{valid_format} are resolvable")
    
    return {
        'total': total,
        'valid_format': valid_format,
        'resolvable': resolvable if check_resolvable else None
    }


def scan_for_missing_dois(output_path: Path) -> Dict[str, int]:
    """
    Scan for missing DOIs and output potential matches.
    
    Args:
        output_path: Path to write the output file
    
    Returns:
        Summary statistics
    """
    scanner = DOIScanner()
    scanner.export_potential_dois(output_path)
    
    # Load the output file to get statistics
    with open(output_path, 'r') as f:
        data = json.load(f)
    
    return {
        'total_scanned': data.get('total', 0),
        'potentials_found': len(data.get('tools', []))
    }


def collect_citations(
    use_cache: bool = True,
    force_refresh: bool = False,
    output_path: Path = IMPACT_DATA_PATH
) -> Dict[str, int]:
    """
    Collect citation data and generate impact_data.json.
    
    Args:
        use_cache: Whether to use cached API responses
        force_refresh: Whether to force refresh all citation data
        output_path: Path to write the output file
    
    Returns:
        Summary statistics
    """
    collector = CitationCollector()
    
    # Override the default output path if specified
    global IMPACT_DATA_PATH
    IMPACT_DATA_PATH = output_path
    
    # Run the full collection process
    impact_data = collector.run_full_collection(
        use_cache=use_cache,
        force_refresh=force_refresh
    )
    
    return {
        'total_tools': impact_data.get('total_tools', 0),
        'tools_with_citations': impact_data.get('tools_with_citations', 0),
        'total_citations': impact_data.get('total_citations', 0)
    }


def main():
    """Command-line interface for the citation system."""
    # Set up logging
    logger = setup_logging("citation_system")
    
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description="Awesome Virome Citation Collection System"
    )
    
    subparsers = parser.add_subparsers(dest="command", help="Command to run")
    
    # Debug command to extract DOIs
    debug_parser = subparsers.add_parser(
        "debug",
        help="Debug tools and DOIs"
    )
    
    # Validate DOIs command
    validate_parser = subparsers.add_parser(
        "validate", 
        help="Validate DOIs"
    )
    validate_parser.add_argument(
        "dois", 
        nargs="*", 
        help="DOIs to validate"
    )
    validate_parser.add_argument(
        "-f", "--file", 
        help="File containing DOIs (one per line)"
    )
    validate_parser.add_argument(
        "-r", "--resolve", 
        action="store_true",
        help="Check if DOIs are resolvable"
    )
    
    # Scan for missing DOIs command
    scan_parser = subparsers.add_parser(
        "scan", 
        help="Scan for missing DOIs"
    )
    scan_parser.add_argument(
        "-o", "--output", 
        default=str(ROOT_DIR / "potential_dois.json"),
        help="Output file for potential DOIs (JSON format)"
    )
    
    # Collect citations command
    collect_parser = subparsers.add_parser(
        "collect", 
        help="Collect citation data"
    )
    collect_parser.add_argument(
        "--no-cache", 
        action="store_true",
        help="Don't use cached API responses"
    )
    collect_parser.add_argument(
        "--force-refresh", 
        action="store_true",
        help="Force refresh all citation data"
    )
    collect_parser.add_argument(
        "-o", "--output", 
        default=str(IMPACT_DATA_PATH),
        help="Output file path (defaults to impact_data.json)"
    )
    
    # Full workflow command
    full_parser = subparsers.add_parser(
        "full", 
        help="Run the full citation workflow"
    )
    full_parser.add_argument(
        "--no-cache", 
        action="store_true",
        help="Don't use cached API responses"
    )
    full_parser.add_argument(
        "--force-refresh", 
        action="store_true",
        help="Force refresh all citation data"
    )
    
    args = parser.parse_args()
    
    # Handle the selected command
    try:
        if args.command == "validate":
            # Collect DOIs from arguments and/or file
            dois = args.dois or []
            
            if args.file:
                try:
                    with open(args.file, 'r') as f:
                        file_dois = [line.strip() for line in f if line.strip()]
                        dois.extend(file_dois)
                except Exception as e:
                    logger.error(f"Error reading DOI file: {e}")
                    return 1
            
            # If no DOIs provided, collect from tools
            if not dois:
                logger.info("No DOIs provided, collecting from tools...")
                collector = CitationCollector()
                tool_dois = collector.collect_tool_dois()
                logger.info(f"Found {len(tool_dois)} tools with DOIs")
                dois = list(tool_dois.values())
                logger.info(f"First 5 DOIs: {dois[:5]}")
                
            if not dois:
                logger.error("No DOIs found in tools")
                return 1
            
            # Validate DOIs
            stats = validate_dois(dois, args.resolve)
            
            logger.info("DOI validation complete")
            log_summary(logger, stats)
        
        elif args.command == "scan":
            log_section(logger, "Scanning for Missing DOIs")
            
            output_path = Path(args.output)
            stats = scan_for_missing_dois(output_path)
            
            logger.info(f"DOI scanning complete, results saved to {output_path}")
            log_summary(logger, stats)
        
        elif args.command == "collect":
            log_section(logger, "Collecting Citation Data")
            
            output_path = Path(args.output)
            stats = collect_citations(
                use_cache=not args.no_cache,
                force_refresh=args.force_refresh,
                output_path=output_path
            )
            
            logger.info(f"Citation collection complete, results saved to {output_path}")
            log_summary(logger, stats)
        
        elif args.command == "debug":
            log_section(logger, "Debugging Tool DOIs")
            
            collector = CitationCollector()
            tool_dois = collector.collect_tool_dois()
            
            logger.info(f"Found {len(tool_dois)} tools with DOIs")
            logger.info(f"DOI list: {list(tool_dois.values())[:20]}")
            
            # Find empty or bad DOIs
            empty_dois = [name for name, doi in tool_dois.items() if not doi or doi.strip() == ""]
            logger.info(f"Tools with empty DOIs: {len(empty_dois)}")
            if empty_dois:
                logger.info(f"Examples: {empty_dois[:5]}")
            
            # Log a sample of tools and their DOIs
            logger.info("Sample of tools and their DOIs:")
            sample = list(tool_dois.items())[:10]
            for name, doi in sample:
                logger.info(f"  {name}: {doi}")
            
            return 0
            
        elif args.command == "full":
            log_section(logger, "Running Full Citation Workflow")
            
            # Step 1: Scan for missing DOIs
            logger.info("Step 1/3: Scanning for missing DOIs")
            scan_output_path = ROOT_DIR / "potential_dois.json"
            scan_stats = scan_for_missing_dois(scan_output_path)
            log_summary(logger, scan_stats)
            
            # Step 2: Validate known DOIs
            logger.info("Step 2/3: Validating existing DOIs")
            collector = CitationCollector()
            tool_dois = collector.collect_tool_dois()
            validate_stats = validate_dois(list(tool_dois.values()))
            log_summary(logger, validate_stats)
            
            # Step 3: Collect citation data
            logger.info("Step 3/3: Collecting citation data")
            collect_stats = collect_citations(
                use_cache=not args.no_cache,
                force_refresh=args.force_refresh
            )
            log_summary(logger, collect_stats)
            
            logger.info("Full citation workflow complete")
        
        else:
            parser.print_help()
            return 1
        
        return 0
    
    except Exception as e:
        logger.error(f"Error running command: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())