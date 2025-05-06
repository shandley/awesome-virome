#!/usr/bin/env python3
"""
Example script demonstrating how to use the iCite client.
"""

import argparse
import json
import logging
import sys
from pathlib import Path

# Add parent directory to path for imports
parent_dir = Path(__file__).resolve().parent.parent
sys.path.append(str(parent_dir.parent))

from scripts.citation_system.api.sources.icite_client import ICiteClient
from scripts.citation_system.config import ICITE_API_URL, ICITE_RATE_LIMIT
from scripts.citation_system.utils.logging_utils import setup_logging


def main():
    """Main entry point for the example script."""
    # Set up logging
    logger = setup_logging("icite_example")
    
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description="iCite Client Example"
    )
    parser.add_argument(
        "doi",
        help="DOI to look up"
    )
    parser.add_argument(
        "--output",
        help="Output file for the citation data (JSON)"
    )
    parser.add_argument(
        "--pmid",
        action="store_true",
        help="Interpret the input as a PubMed ID instead of a DOI"
    )
    parser.add_argument(
        "--no-cache",
        action="store_true",
        help="Don't use cached API responses"
    )
    
    args = parser.parse_args()
    
    try:
        # Create iCite client
        client = ICiteClient(
            api_url=ICITE_API_URL,
            rate_limit=ICITE_RATE_LIMIT
        )
        
        logger.info(f"Looking up {'PMID' if args.pmid else 'DOI'}: {args.doi}")
        
        # Get citation data
        if args.pmid:
            success, data = client.search_by_pmid(args.doi, use_cache=not args.no_cache)
        else:
            success, data = client.get_citation_data(args.doi, use_cache=not args.no_cache)
        
        if success:
            # Format the data for display
            pretty_data = json.dumps(data, indent=2)
            
            # Print citation count and other key metrics
            logger.info(f"Citation count: {data.get('total_citations', 0)}")
            
            if 'rcr' in data:
                logger.info(f"Relative Citation Ratio (RCR): {data.get('rcr')}")
            
            if 'expected_citations' in data:
                logger.info(f"Expected citations per year: {data.get('expected_citations')}")
            
            logger.info(f"Full citation data:\n{pretty_data}")
            
            # Write to output file if specified
            if args.output:
                with open(args.output, 'w') as f:
                    json.dump(data, f, indent=2)
                logger.info(f"Citation data written to {args.output}")
        else:
            logger.error(f"Failed to get citation data: {data}")
            return 1
        
        return 0
    
    except Exception as e:
        logger.error(f"Error: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())