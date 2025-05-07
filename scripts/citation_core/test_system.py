#!/usr/bin/env python3
"""
Test Script for Citation Core System

This script performs a basic test of the citation core system by:
1. Testing the iCite API client with a known DOI
2. Testing DOI extraction and validation
3. Testing the citation collector and impact data builder with sample data

Run this script to verify that the citation system is working correctly.
"""

import os
import sys
import json
import logging
from typing import Dict, Any

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger("test_system")

# Import local modules
from icite_api import ICiteClient, test_icite_client
from doi_handler import test_doi_handler


def create_sample_tools() -> str:
    """Create a sample data.json file for testing."""
    sample_file = "sample_data.json"
    
    sample_data = {
        "nodes": [
            {
                "name": "CheckV",
                "type": "tool",
                "url": "https://github.com/nayfach/checkv",
                "doi": "10.1038/s41587-020-00774-7",
                "description": "CheckV is a pipeline for assessing the quality of metagenome-assembled viral genomes."
            },
            {
                "name": "VirSorter2",
                "type": "tool",
                "url": "https://github.com/jiarong/VirSorter2",
                "doi": "10.1186/s40168-020-00990-y",
                "description": "VirSorter2 applies a multi-classifier, expert-guided approach to detect diverse DNA and RNA virus genomes."
            },
            {
                "name": "Pharokka",
                "type": "tool",
                "url": "https://github.com/gbouras13/pharokka",
                "doi": "10.1093/bioinformatics/btad311",
                "description": "Pharokka is a fast, scalable, and comprehensive pipeline to annotate phage genomes."
            },
            {
                "name": "Test Tool Without DOI",
                "type": "tool",
                "url": "https://example.com/tool",
                "description": "This is a test tool without a DOI."
            }
        ]
    }
    
    with open(sample_file, 'w') as f:
        json.dump(sample_data, f, indent=2)
    
    logger.info(f"Created sample data file: {sample_file}")
    return sample_file


def test_citation_collector(data_file: str) -> bool:
    """Test the citation collector."""
    from citation_collector import CitationCollector
    
    output_dir = "test_citations"
    os.makedirs(output_dir, exist_ok=True)
    
    logger.info("Testing CitationCollector...")
    
    # Initialize collector with sample data
    collector = CitationCollector(
        data_path=data_file,
        output_dir=output_dir,
        icite_rate_limit=1.0  # Fast for testing
    )
    
    # Collect citations
    results = collector.collect_citations()
    
    # Check if we got results
    if not results:
        logger.error("Citation collector returned no results")
        return False
    
    logger.info(f"Citation collector returned data for {len(results)} tools")
    
    # Print results summary
    for tool_name, data in results.items():
        citations = data.get("total_citations", 0)
        logger.info(f"Tool: {tool_name}, Citations: {citations}")
    
    return True


def test_impact_data_builder(data_file: str) -> bool:
    """Test the impact data builder."""
    from impact_data_builder import ImpactDataBuilder
    
    citation_dir = "test_citations"
    output_path = "test_impact_data.json"
    
    logger.info("Testing ImpactDataBuilder...")
    
    # Initialize builder
    builder = ImpactDataBuilder(
        data_path=data_file,
        citation_dir=citation_dir,
        output_path=output_path
    )
    
    # Build and save impact data
    success = builder.build_and_save()
    
    if not success:
        logger.error("Impact data builder failed")
        return False
    
    # Verify the output file exists
    if not os.path.exists(output_path):
        logger.error(f"Output file not created: {output_path}")
        return False
    
    # Load and check the impact data
    try:
        with open(output_path, 'r') as f:
            impact_data = json.load(f)
        
        # Check required fields
        required_fields = ['last_updated', 'tools', 'total_tools', 'tools_with_citations', 'total_citations']
        for field in required_fields:
            if field not in impact_data:
                logger.error(f"Missing required field in impact data: {field}")
                return False
        
        logger.info(f"Impact data successfully built with {len(impact_data['tools'])} tools")
        return True
        
    except (json.JSONDecodeError, IOError) as e:
        logger.error(f"Error loading impact data: {e}")
        return False


def run_full_test():
    """Run a full test of the citation system."""
    logger.info("Starting full system test")
    
    # Step 1: Test iCite API client
    logger.info("\n=== Testing iCite API Client ===")
    if not test_icite_client():
        logger.error("iCite API client test failed")
        return False
    
    # Step 2: Test DOI handler
    logger.info("\n=== Testing DOI Handler ===")
    if not test_doi_handler():
        logger.error("DOI handler test failed")
        return False
    
    # Step 3: Create sample data
    logger.info("\n=== Creating Sample Data ===")
    data_file = create_sample_tools()
    
    # Step 4: Test citation collector
    logger.info("\n=== Testing Citation Collector ===")
    if not test_citation_collector(data_file):
        logger.error("Citation collector test failed")
        return False
    
    # Step 5: Test impact data builder
    logger.info("\n=== Testing Impact Data Builder ===")
    if not test_impact_data_builder(data_file):
        logger.error("Impact data builder test failed")
        return False
    
    logger.info("\n=== All tests passed! ===")
    logger.info("The citation core system is working correctly.")
    return True


if __name__ == "__main__":
    success = run_full_test()
    sys.exit(0 if success else 1)