#!/usr/bin/env python3
"""
DOI Handler for Awesome-Virome

This module provides tools for extracting and validating Digital Object Identifiers
(DOIs) from tool metadata. It includes functionality for:

- Validating DOI syntax
- Extracting DOIs from various text formats
- Normalizing DOIs to a consistent format
"""

import re
import logging
from typing import List, Optional, Dict, Any, Set

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger("doi_handler")

# Regular expression for DOI validation
# DOIs typically start with "10." followed by a number (registrant code),
# then a slash, and then an alphanumeric string (suffix)
DOI_PATTERN = re.compile(r"^10\.\d{4,9}/[-._;()/:a-zA-Z0-9]+$")

# Regular expression for finding DOIs in text
DOI_EXTRACT_PATTERN = re.compile(r"(?:doi:|https?://doi\.org/)?(10\.\d{4,9}/[-._;()/:a-zA-Z0-9]+)")


def is_valid_doi(doi: str) -> bool:
    """
    Check if a string is a valid DOI.
    
    Args:
        doi: The DOI string to validate
        
    Returns:
        True if the DOI is valid, False otherwise
    """
    if not doi or not isinstance(doi, str):
        return False
    
    # Clean the DOI before validation
    cleaned_doi = clean_doi(doi)
    
    # Check if the cleaned DOI matches the pattern
    return bool(DOI_PATTERN.match(cleaned_doi))


def clean_doi(doi: str) -> str:
    """
    Clean and normalize a DOI string.
    
    Args:
        doi: The DOI string to clean
        
    Returns:
        Cleaned DOI string
    """
    if not doi or not isinstance(doi, str):
        return ""
    
    # Remove whitespace
    doi = doi.strip()
    
    # Remove common prefixes
    if doi.lower().startswith("doi:"):
        doi = doi[4:].strip()
    
    if doi.lower().startswith("https://doi.org/"):
        doi = doi[16:].strip()
    
    if doi.lower().startswith("http://doi.org/"):
        doi = doi[15:].strip()
    
    return doi


def extract_doi_from_text(text: str) -> Optional[str]:
    """
    Extract a DOI from text content.
    
    Args:
        text: The text to search for a DOI
        
    Returns:
        Extracted DOI if found, None otherwise
    """
    if not text or not isinstance(text, str):
        return None
    
    # Search for DOI pattern in the text
    matches = DOI_EXTRACT_PATTERN.findall(text)
    
    if not matches:
        return None
    
    # Return the first valid DOI found
    for match in matches:
        cleaned_doi = clean_doi(match)
        if is_valid_doi(cleaned_doi):
            return cleaned_doi
    
    return None


def extract_doi_from_tool(tool: Dict[str, Any]) -> Optional[str]:
    """
    Extract DOI from tool metadata.
    
    Args:
        tool: Dictionary containing tool metadata
        
    Returns:
        Extracted DOI if found, None otherwise
    """
    # Check explicit DOI field
    if "doi" in tool and tool["doi"]:
        doi = clean_doi(tool["doi"])
        if is_valid_doi(doi):
            return doi
    
    # Check citation fields if present
    if "citation" in tool and tool["citation"]:
        citation = tool["citation"]
        
        # If citation is a string, try to extract DOI
        if isinstance(citation, str):
            doi = extract_doi_from_text(citation)
            if doi:
                return doi
        
        # If citation is a dictionary, check for DOI field
        elif isinstance(citation, dict) and "doi" in citation:
            doi = clean_doi(citation["doi"])
            if is_valid_doi(doi):
                return doi
    
    # Check description for DOI
    if "description" in tool and tool["description"]:
        doi = extract_doi_from_text(tool["description"])
        if doi:
            return doi
    
    return None


def extract_dois_from_tools(tools: List[Dict[str, Any]]) -> Dict[str, str]:
    """
    Extract DOIs from a list of tools.
    
    Args:
        tools: List of tool metadata dictionaries
        
    Returns:
        Dictionary mapping tool names to their DOIs
    """
    results = {}
    tools_with_doi = 0
    
    for tool in tools:
        tool_name = tool.get("name", "")
        if not tool_name:
            continue
        
        doi = extract_doi_from_tool(tool)
        if doi:
            results[tool_name] = doi
            tools_with_doi += 1
    
    logger.info(f"Extracted DOIs for {tools_with_doi} out of {len(tools)} tools")
    return results


def test_doi_handler():
    """Test the DOI handler functionality."""
    # Test valid DOIs
    valid_dois = [
        "10.1093/nar/gkaa939",
        "10.1038/s41586-020-2008-3",
        "10.1093/bioinformatics/btab213"
    ]
    
    for doi in valid_dois:
        assert is_valid_doi(doi), f"Should be valid: {doi}"
    
    # Test DOI cleaning
    doi_variants = [
        "doi:10.1093/nar/gkaa939",
        "https://doi.org/10.1093/nar/gkaa939",
        "  10.1093/nar/gkaa939  "
    ]
    
    for variant in doi_variants:
        cleaned = clean_doi(variant)
        assert cleaned == "10.1093/nar/gkaa939", f"Failed to clean: {variant} -> {cleaned}"
    
    # Test DOI extraction from text
    text_samples = [
        "The tool was published in Nature (doi:10.1038/s41586-020-2008-3) and has been widely cited.",
        "For more information, visit https://doi.org/10.1093/nar/gkaa939",
        "Reference: Smith et al., 2020. Journal of Virology. 10.1093/bioinformatics/btab213"
    ]
    
    expected_dois = [
        "10.1038/s41586-020-2008-3",
        "10.1093/nar/gkaa939",
        "10.1093/bioinformatics/btab213"
    ]
    
    for i, text in enumerate(text_samples):
        extracted = extract_doi_from_text(text)
        assert extracted == expected_dois[i], f"Failed to extract from: {text} -> {extracted}"
    
    logger.info("All DOI handler tests passed!")
    return True


if __name__ == "__main__":
    # Run tests if executed directly
    test_doi_handler()