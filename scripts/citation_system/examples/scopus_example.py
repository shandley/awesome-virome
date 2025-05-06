#!/usr/bin/env python3
"""
Example usage of the Scopus citation client.

This script demonstrates how to use the Scopus API client to retrieve citation data
for scientific publications using DOIs.

To run this example, you need a valid Scopus API key. Set it as an environment variable:
export SCOPUS_API_KEY="your-api-key-here"

Optional: If you have an institutional token, you can set that as well:
export SCOPUS_INSTITUTIONAL_TOKEN="your-token-here"
"""

import os
import sys
import json
import logging
from pathlib import Path

# Add parent directory to path to allow imports
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent.parent))

from scripts.citation_system.api.sources.scopus_client import ScopusClient
from scripts.citation_system.config import SCOPUS_API_URL, SCOPUS_API_KEY, SCOPUS_INSTITUTIONAL_TOKEN

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger("scopus_example")

def main():
    """Run the Scopus client example."""
    # Check if API key is available
    api_key = os.environ.get("SCOPUS_API_KEY", SCOPUS_API_KEY)
    if not api_key:
        logger.error("No Scopus API key found. Please set SCOPUS_API_KEY environment variable.")
        return 1
    
    # Initialize the Scopus client
    client = ScopusClient(
        api_url=SCOPUS_API_URL,
        api_key=api_key,
        institutional_token=os.environ.get("SCOPUS_INSTITUTIONAL_TOKEN", SCOPUS_INSTITUTIONAL_TOKEN),
        rate_limit=5.0
    )
    
    # Example DOIs to look up
    # These are well-known papers that should be in Scopus
    dois = [
        "10.1038/s41586-020-2008-3",  # COVID-19 paper
        "10.1093/bioinformatics/btab213",  # Bioinformatics paper
        "10.1126/science.1132563"  # Highly cited Science paper
    ]
    
    print(f"Looking up citation data for {len(dois)} DOIs using Scopus API...")
    
    for doi in dois:
        print(f"\nLooking up DOI: {doi}")
        success, data = client.get_citation_data(doi)
        
        if success:
            print(f"✅ Found in Scopus")
            print(f"  Total citations: {data['total_citations']}")
            
            if "metadata" in data:
                metadata = data["metadata"]
                print(f"  Title: {metadata.get('title', 'Unknown')}")
                print(f"  Journal: {metadata.get('journal', 'Unknown')}")
                print(f"  Year: {metadata.get('year', 'Unknown')}")
                print(f"  Authors: {', '.join(metadata.get('authors', []))[:100]}")
            
            if "citations_by_year" in data and data["citations_by_year"]:
                print("  Citations by year:")
                for year, count in sorted(data["citations_by_year"].items()):
                    print(f"    {year}: {count}")
        else:
            print(f"❌ Error: {data}")
    
    print("\nBatch citation lookup example:")
    results = client.batch_get_citation_data(dois)
    
    print(f"Batch results: {len(results)} DOIs processed")
    total_citations = sum(data.get("total_citations", 0) for data in results.values() if "error" not in data)
    print(f"Total citations across all papers: {total_citations}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())