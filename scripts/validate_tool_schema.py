#!/usr/bin/env python3
"""
Tool Schema Validator for Awesome-Virome Repository

This script validates the data.json file against the defined JSON schema
for tool entries, ensuring data consistency and quality.
"""

import os
import sys
import json
import logging
import argparse
from pathlib import Path
from typing import Dict, List, Any, Optional
import jsonschema

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Constants
ROOT_DIR = Path(__file__).parent.parent
DATA_JSON_PATH = ROOT_DIR / "data.json"
SCHEMA_PATH = ROOT_DIR / "scripts" / "schema.json"


def load_schema(schema_path: Path) -> Dict[str, Any]:
    """Load and parse JSON schema."""
    try:
        with open(schema_path, 'r') as f:
            schema = json.load(f)
        return schema
    except (FileNotFoundError, json.JSONDecodeError) as e:
        logger.error(f"Error loading schema: {e}")
        sys.exit(1)


def load_data_json(data_json_path: Path) -> Dict[str, Any]:
    """Load and parse data.json file."""
    try:
        with open(data_json_path, 'r') as f:
            data = json.load(f)
        return data
    except (FileNotFoundError, json.JSONDecodeError) as e:
        logger.error(f"Error loading data.json: {e}")
        sys.exit(1)


def extract_tools_from_data_json(data: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Extract tool nodes from data.json."""
    tools = []
    for node in data.get('nodes', []):
        if node.get('type') == 'tool':
            tools.append(node)
    return tools


def validate_tool_against_schema(tool: Dict[str, Any], schema: Dict[str, Any]) -> List[str]:
    """
    Validate a single tool against the JSON schema.
    
    Args:
        tool: Tool data to validate
        schema: JSON schema
        
    Returns:
        List of validation errors, empty if valid
    """
    validator = jsonschema.Draft7Validator(schema)
    errors = list(validator.iter_errors(tool))
    
    return [
        f"{error.path}: {error.message}" if error.path else error.message
        for error in errors
    ]


def validate_all_tools(tools: List[Dict[str, Any]], schema: Dict[str, Any]) -> Dict[str, List[str]]:
    """
    Validate all tools against schema.
    
    Args:
        tools: List of tool data dictionaries
        schema: JSON schema
        
    Returns:
        Dictionary mapping tool names to validation errors
    """
    results = {}
    
    for tool in tools:
        name = tool.get('name', 'Unknown Tool')
        errors = validate_tool_against_schema(tool, schema)
        
        if errors:
            results[name] = errors
    
    return results


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Validate tool entries against JSON schema')
    parser.add_argument('--data', help='Path to data.json file', default=str(DATA_JSON_PATH))
    parser.add_argument('--schema', help='Path to JSON schema file', default=str(SCHEMA_PATH))
    parser.add_argument('--output', help='Path to write validation results')
    parser.add_argument('--verbose', '-v', action='store_true', help='Enable verbose output')
    
    return parser.parse_args()


def main():
    """Main function."""
    args = parse_arguments()
    
    # Configure logging level
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    
    # Load schema and data
    schema = load_schema(Path(args.schema))
    data = load_data_json(Path(args.data))
    
    # Extract and validate tools
    tools = extract_tools_from_data_json(data)
    logger.info(f"Found {len(tools)} tools in data.json")
    
    validation_results = validate_all_tools(tools, schema)
    
    # Report validation results
    if validation_results:
        logger.error(f"Validation failed for {len(validation_results)} tools:")
        for tool_name, errors in validation_results.items():
            logger.error(f"  {tool_name}:")
            for error in errors:
                logger.error(f"    - {error}")
        
        # Write results to output file if specified
        if args.output:
            try:
                with open(args.output, 'w') as f:
                    json.dump(validation_results, f, indent=2)
                logger.info(f"Validation results written to {args.output}")
            except IOError as e:
                logger.error(f"Error writing output file: {e}")
        
        sys.exit(1)
    else:
        logger.info("All tools passed schema validation")
        sys.exit(0)


if __name__ == "__main__":
    main()