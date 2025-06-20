#!/usr/bin/env python3
"""
Debug script for citation source attribution.

This script checks if citation source attribution is working correctly
and diagnoses any issues with the citation collection process.
"""

import os
import sys
import json
import logging
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler("debug_attribution.log")
    ]
)
logger = logging.getLogger(__name__)

# Constants
BASE_DIR = Path(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
IMPACT_DATA_PATH = BASE_DIR / "impact_data.json"


def check_citation_code():
    """Check if the citation source attribution code is correctly implemented."""
    # Import our key modules to check for attribution fields
    logger.info("Checking citation aggregator implementation...")
    
    try:
        # Add the scripts directory to the path
        sys.path.append(str(BASE_DIR))
        
        # Import the citation aggregator
        from scripts.citation_system.collectors.citation_aggregator import CitationAggregator
        
        # Check if the new fields are in the initialization
        aggregator = CitationAggregator()
        
        # Create a basic aggregation to check field presence
        test_data = {
            "test_source": {
                "total_citations": 100,
                "citations_by_year": {"2020": 50, "2021": 50}
            }
        }
        
        # Perform test aggregation
        result = aggregator._aggregate_citation_data(test_data, "10.1000/test.123")
        
        # Check for attribution fields
        missing_fields = []
        for field in ["citation_source", "yearly_citation_source"]:
            if field not in result:
                missing_fields.append(field)
        
        if missing_fields:
            logger.error(f"Missing attribution fields in aggregator output: {missing_fields}")
            return False
        else:
            logger.info("Citation aggregator implementation looks correct!")
            logger.info(f"Test result contains attribution fields: {result['citation_source']} and {result['yearly_citation_source']}")
            return True
        
    except Exception as e:
        logger.error(f"Error checking citation code: {e}")
        return False


def check_impact_data():
    """Check if the impact_data.json file has citation source attribution."""
    logger.info(f"Checking impact data file at {IMPACT_DATA_PATH}...")
    
    if not IMPACT_DATA_PATH.exists():
        logger.error(f"Impact data file not found: {IMPACT_DATA_PATH}")
        return False
    
    try:
        with open(IMPACT_DATA_PATH, 'r') as f:
            impact_data = json.load(f)
        
        # Check if we have tools
        tools = impact_data.get('tools', [])
        if not tools:
            logger.error("No tools found in impact_data.json")
            return False
        
        logger.info(f"Found {len(tools)} tools in impact_data.json")
        
        # Check for citation sources in the summary
        if 'citation_sources' in impact_data:
            logger.info(f"Found citation_sources in summary: {impact_data['citation_sources']}")
        else:
            logger.warning("No citation_sources found in summary!")
        
        # Check for attribution fields in tools
        tools_with_attribution = 0
        for tool in tools:
            if 'citation_source' in tool or 'yearly_citation_source' in tool:
                tools_with_attribution += 1
        
        if tools_with_attribution > 0:
            logger.info(f"Found {tools_with_attribution}/{len(tools)} tools with source attribution")
            return True
        else:
            logger.error("No tools have source attribution!")
            
            # Show sample of tools for debugging
            logger.info("First 3 tools structure:")
            for i, tool in enumerate(tools[:3]):
                logger.info(f"Tool {i+1}: {json.dumps(tool, indent=2)}")
            
            return False
        
    except Exception as e:
        logger.error(f"Error reading impact data: {e}")
        return False


def check_collection_script():
    """Check if the collect_citations.py script has the correct command structure."""
    logger.info("Checking citation collection script...")
    
    collect_script_path = BASE_DIR / "scripts" / "citation_system" / "collect_citations.py"
    
    if not collect_script_path.exists():
        logger.error(f"Collection script not found: {collect_script_path}")
        return False
    
    try:
        with open(collect_script_path, 'r') as f:
            script_content = f.read()
        
        # Check for the command parsing
        if "ArgumentParser" in script_content:
            logger.info("Found ArgumentParser in collection script")
            
            # Check for 'full' command
            if "full" in script_content:
                logger.info("Found 'full' command in script")
            else:
                logger.warning("'full' command not found in script - workflow may fail!")
                
                # Try to determine what commands are available
                import re
                commands = re.findall(r'parser_(\w+)\s*=\s*subparsers\.add_parser', script_content)
                logger.info(f"Available commands appear to be: {commands}")
        
        return True
        
    except Exception as e:
        logger.error(f"Error checking collection script: {e}")
        return False


def main():
    """Main function to run diagnostics."""
    logger.info("=========== Citation Attribution Diagnostics ===========")
    
    # Check if the code is correctly implemented
    code_check = check_citation_code()
    
    # Check if the impact data has attribution
    data_check = check_impact_data()
    
    # Check the collection script
    script_check = check_collection_script()
    
    # Summarize findings
    logger.info("\n=========== Diagnostic Summary ===========")
    logger.info(f"Citation attribution code: {'OK' if code_check else 'ISSUES FOUND'}")
    logger.info(f"Impact data attribution: {'OK' if data_check else 'ISSUES FOUND'}")
    logger.info(f"Collection script: {'OK' if script_check else 'ISSUES FOUND'}")
    
    # Provide recommendation
    if not code_check or not data_check:
        logger.info("\nRECOMMENDATION:")
        logger.info("1. Ensure the citation_aggregator.py changes are properly deployed")
        logger.info("2. Check that the collect_citations.py command structure matches the workflow")
        logger.info("3. Try running a local test collection and monitor attribution fields")
    
    return 0 if (code_check and data_check and script_check) else 1


if __name__ == "__main__":
    sys.exit(main())