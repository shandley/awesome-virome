#!/usr/bin/env python3
"""
Utility functions for schema validation.
"""

import json
import logging
from pathlib import Path
from typing import Any, Dict, Optional, Union

import jsonschema

logger = logging.getLogger(__name__)


def load_schema(schema_path: Path) -> Dict[str, Any]:
    """
    Load a JSON schema from a file.
    
    Args:
        schema_path: Path to the schema file
    
    Returns:
        Loaded schema as a dictionary
    
    Raises:
        FileNotFoundError: If the schema file doesn't exist
        json.JSONDecodeError: If the schema file contains invalid JSON
    """
    try:
        with open(schema_path, 'r') as f:
            return json.load(f)
    except FileNotFoundError:
        logger.error(f"Schema file not found: {schema_path}")
        raise
    except json.JSONDecodeError:
        logger.error(f"Invalid JSON in schema file: {schema_path}")
        raise


def validate_against_schema(
    data: Dict[str, Any], 
    schema: Union[Dict[str, Any], Path],
    error_title: Optional[str] = None
) -> bool:
    """
    Validate data against a JSON schema.
    
    Args:
        data: Data to validate
        schema: JSON schema as a dictionary or path to schema file
        error_title: Optional title for error messages
    
    Returns:
        True if validation succeeds, False otherwise
    """
    if isinstance(schema, Path):
        schema = load_schema(schema)
    
    try:
        jsonschema.validate(instance=data, schema=schema)
        return True
    except jsonschema.exceptions.ValidationError as e:
        title = error_title or "Schema validation failed"
        logger.error(f"{title}: {e.message}")
        # Log the path where validation failed
        if e.path:
            path_str = ".".join(str(p) for p in e.path)
            logger.error(f"Failed at path: {path_str}")
        return False


def check_schema(file_path: Path, schema_path: Path) -> bool:
    """
    Check if a JSON file conforms to a schema.
    
    Args:
        file_path: Path to the JSON file to validate
        schema_path: Path to the schema file
    
    Returns:
        True if validation succeeds, False otherwise
    """
    try:
        with open(file_path, 'r') as f:
            data = json.load(f)
        schema = load_schema(schema_path)
        return validate_against_schema(
            data, 
            schema, 
            f"File {file_path} failed schema validation"
        )
    except (FileNotFoundError, json.JSONDecodeError) as e:
        logger.error(f"Error validating {file_path}: {e}")
        return False