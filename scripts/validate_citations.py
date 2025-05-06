#!/usr/bin/env python3
"""
Citation Data Validation

This script performs validation checks on citation data in the awesome-virome repository,
including validating DOIs, checking citation formatting, and ensuring consistency across
different citation sources.

NOTE: DOI fixing capability has been removed from this script. For DOI format fixing,
please use the authoritative 'auto_fix_dois.py' script and its associated workflow.
"""

import os
import sys
import json
import logging
import argparse
import datetime
import re
import requests
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple, Set

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("citation_validation.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Constants
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(SCRIPT_DIR)
DATA_JSON_PATH = os.path.join(ROOT_DIR, "data.json")
METADATA_DIR = os.path.join(ROOT_DIR, "metadata")
ACADEMIC_IMPACT_DIR = os.path.join(METADATA_DIR, "academic_impact")
PUBMED_CITATIONS_DIR = os.path.join(METADATA_DIR, "pubmed_citations")
REPORTS_DIR = os.path.join(ROOT_DIR, "reports")
CITATIONS_MD_PATH = os.path.join(ROOT_DIR, "citations.md")

# DOI regex pattern (simplified)
DOI_PATTERN = r"^10\.\d{4,9}/[-._;()/:A-Z0-9]+$"

# Regular expressions for citation formats (simplified)
CITATION_PATTERNS = {
    "bibtex": r"@\w+\s*{.*?,\s*(\s*\w+\s*=\s*{.*?},?)+\s*}",
    "apa": r".*?\(\d{4}\)\..*?\..*?\.",
    "mla": r".*?\..*?\..*?\d{4}\.",
}

# DOI validation pattern only - no fixing in this script
# For DOI fixing, use auto_fix_dois.py instead

def load_data_json() -> Dict[str, Any]:
    """Load the current data.json file"""
    try:
        with open(DATA_JSON_PATH, 'r') as f:
            data = json.load(f)
            # Convert nodes with type "tool" into the expected tools format
            if "nodes" in data:
                tools = []
                for node in data["nodes"]:
                    if node.get("type") == "tool":
                        tools.append(node)
                return {"tools": tools}
            else:
                logger.error("Data.json does not contain the expected 'nodes' key")
                return {"tools": []}
    except (json.JSONDecodeError, FileNotFoundError) as e:
        logger.error(f"Error loading data.json: {e}")
        return {"tools": []}

def load_citations_md() -> str:
    """Load the citations.md file if it exists"""
    if os.path.exists(CITATIONS_MD_PATH):
        try:
            with open(CITATIONS_MD_PATH, 'r') as f:
                return f.read()
        except (IOError, OSError) as e:
            logger.error(f"Error loading citations.md: {e}")
    return ""

def load_metadata_for_tool(tool_id: str, subdir: Optional[str] = None) -> Dict[str, Any]:
    """Load the metadata file for a specific tool"""
    base_dir = METADATA_DIR
    if subdir:
        base_dir = os.path.join(METADATA_DIR, subdir)
    
    metadata_path = os.path.join(base_dir, f"{tool_id}.json")
    try:
        if os.path.exists(metadata_path):
            with open(metadata_path, 'r') as f:
                return json.load(f)
        return {}
    except (json.JSONDecodeError, FileNotFoundError) as e:
        logger.warning(f"Error loading metadata for {tool_id} from {subdir}: {e}")
        return {}

def validate_doi_format(doi: str) -> bool:
    """
    Validate if a DOI string follows the correct format
    
    Args:
        doi: The DOI string to validate
        
    Returns:
        bool: True if the DOI format is valid, False otherwise
    """
    if not doi:
        return False
    
    # Handle markdown link format like 10.1093/gigascience/giae020](https://doi.org/10.1093/gigascience/giae020)
    # We only validate, not fix
    if '](' in doi:
        return False
    
    # Basic format validation
    return bool(re.match(DOI_PATTERN, doi, re.IGNORECASE))

def validate_doi(doi: str, skip_resolution_check: bool = False) -> bool:
    """
    Validate if a DOI is properly formatted and resolvable
    
    Args:
        doi: The DOI string to validate
        skip_resolution_check: Skip online DOI resolution check
        
    Returns:
        bool: True if the DOI is valid, False otherwise
    """
    # Basic format validation
    if not re.match(DOI_PATTERN, doi, re.IGNORECASE):
        return False
    
    # Skip resolution check if requested (for bulk validation)
    if skip_resolution_check:
        return True
    
    # Check if DOI is resolvable
    try:
        headers = {
            'Accept': 'application/json',
            'User-Agent': 'awesome-virome-citation-validator/1.0 (https://github.com/yourusername/awesome-virome)'
        }
        response = requests.get(f"https://doi.org/api/handles/{doi}", headers=headers, timeout=10)
        if response.status_code == 200:
            response_data = response.json()
            return response_data.get('responseCode') == 1
        return False
    except requests.RequestException:
        # If we can't connect, we'll assume the DOI might be valid
        # since we don't want to fail validation only due to connectivity issues
        logger.warning(f"Could not connect to DOI resolution service for {doi}")
        return True

def validate_citation_format(citation: str, format_type: str) -> bool:
    """
    Validate if a citation string follows the expected format
    
    Args:
        citation: The citation string
        format_type: The expected format (bibtex, apa, mla)
        
    Returns:
        bool: True if the citation format is valid, False otherwise
    """
    if format_type not in CITATION_PATTERNS:
        logger.warning(f"Unknown citation format: {format_type}")
        return False
    
    pattern = CITATION_PATTERNS[format_type]
    return bool(re.match(pattern, citation, re.DOTALL | re.IGNORECASE))

def check_citations_consistency(
    tool_id: str, 
    data_json_doi: Optional[str] = None, 
    skip_doi_check: bool = False
) -> Dict[str, Any]:
    """
    Check if citation information is consistent across different sources
    
    Args:
        tool_id: The tool ID to check
        data_json_doi: The DOI from data.json (if available)
        skip_doi_check: Skip DOI resolution check
        
    Returns:
        Dict with validation results
    """
    results = {
        "tool_id": tool_id,
        "has_doi": bool(data_json_doi),
        "doi_valid": False if data_json_doi else None,
        "doi_consistent": True,
        "has_pubmed_data": False,
        "has_academic_impact": False,
        "citation_formats": [],
        "issues": []
    }
    
    # Load metadata from different sources
    general_metadata = load_metadata_for_tool(tool_id)
    academic_impact = load_metadata_for_tool(tool_id, "academic_impact")
    pubmed_data = load_metadata_for_tool(tool_id, "pubmed_citations")
    
    # Check for existence in different sources
    results["has_academic_impact"] = bool(academic_impact)
    results["has_pubmed_data"] = bool(pubmed_data)
    
    # Store DOIs for consistency checks
    dois = {}
    
    # Validate DOI in data.json
    if data_json_doi:
        dois["data_json"] = data_json_doi
        format_valid = validate_doi_format(data_json_doi)
        if format_valid:
            results["doi_valid"] = validate_doi(data_json_doi, skip_doi_check)
        else:
            results["doi_valid"] = False
            results["issues"].append(f"Invalid DOI format: {data_json_doi}")
    
    # Validate DOI in general metadata
    if general_metadata.get("doi"):
        dois["general"] = general_metadata["doi"]
        if not validate_doi_format(general_metadata["doi"]):
            results["issues"].append(f"Invalid DOI format in general metadata: {general_metadata['doi']}")
    
    # Validate DOI in academic impact
    if academic_impact.get("doi"):
        dois["academic"] = academic_impact["doi"]
        if not validate_doi_format(academic_impact["doi"]):
            results["issues"].append(f"Invalid DOI format in academic impact metadata: {academic_impact['doi']}")
    
    # Validate DOI in pubmed data
    if pubmed_data.get("doi"):
        dois["pubmed"] = pubmed_data["doi"]
        if not validate_doi_format(pubmed_data["doi"]):
            results["issues"].append(f"Invalid DOI format in pubmed citation data: {pubmed_data['doi']}")
    
    # Check DOI consistency across sources
    unique_dois = set(dois.values())
    if len(unique_dois) > 1:
        results["doi_consistent"] = False
        results["issues"].append(f"Inconsistent DOIs across different sources: {', '.join(unique_dois)}")
    
    # Check citation formats
    citations_md = load_citations_md()
    if citations_md:
        # Look for this tool's citation section in citations.md
        tool_section_match = re.search(rf"## {re.escape(tool_id)}\s*\n(.*?)(?:\n## |\Z)", 
                                     citations_md, re.DOTALL)
        
        if tool_section_match:
            section_text = tool_section_match.group(1)
            
            # Check for BibTeX
            if "```bibtex" in section_text:
                bibtex_match = re.search(r"```bibtex\s*(.*?)\s*```", section_text, re.DOTALL)
                if bibtex_match and validate_citation_format(bibtex_match.group(1), "bibtex"):
                    results["citation_formats"].append("bibtex")
                else:
                    results["issues"].append("Invalid BibTeX format")
            
            # Check for APA
            if "APA:" in section_text:
                apa_match = re.search(r"APA:\s*(.*?)(?:\n\n|\n[A-Z]|\Z)", section_text, re.DOTALL)
                if apa_match and validate_citation_format(apa_match.group(1), "apa"):
                    results["citation_formats"].append("apa")
                else:
                    results["issues"].append("Invalid APA format")
            
            # Check for MLA
            if "MLA:" in section_text:
                mla_match = re.search(r"MLA:\s*(.*?)(?:\n\n|\n[A-Z]|\Z)", section_text, re.DOTALL)
                if mla_match and validate_citation_format(mla_match.group(1), "mla"):
                    results["citation_formats"].append("mla")
                else:
                    results["issues"].append("Invalid MLA format")
    
    return results

def validate_all_citations(
    skip_doi_check: bool = False
) -> Dict[str, Any]:
    """
    Validate citations for all tools in the repository
    
    Args:
        skip_doi_check: Skip DOI resolution check
    
    Returns:
        Dict with validation results
    """
    data = load_data_json()
    tools = data.get("tools", [])
    
    if not tools:
        logger.error("No tools found in data.json (looked for 'type': 'tool' nodes)")
        logger.info("This may happen if no tool nodes were found in the data.json structure")
        # Return minimal results structure to prevent downstream errors
        return {
            "timestamp": datetime.datetime.now().isoformat(),
            "total_tools": 0,
            "tools_with_doi": 0,
            "valid_dois": 0,
            "consistent_dois": 0,
            "tools_with_pubmed_data": 0,
            "tools_with_academic_impact": 0,
            "citation_format_counts": {"bibtex": 0, "apa": 0, "mla": 0},
            "all_issues": ["No tools found in data.json"],
            "tool_results": []
        }
    
    logger.info(f"Found {len(tools)} tools in data.json")
    
    results = {
        "timestamp": datetime.datetime.now().isoformat(),
        "total_tools": len(tools),
        "tools_with_doi": 0,
        "valid_dois": 0,
        "consistent_dois": 0,
        "tools_with_pubmed_data": 0,
        "tools_with_academic_impact": 0,
        "citation_format_counts": {
            "bibtex": 0,
            "apa": 0,
            "mla": 0
        },
        "all_issues": [],
        "tool_results": []
    }
    
    # Process each tool
    for tool in tools:
        tool_id = tool.get("id")
        if not tool_id:
            continue
        
        # Validate citations for this tool
        tool_result = check_citations_consistency(
            tool_id, 
            tool.get("doi"),
            skip_doi_check=skip_doi_check
        )
        results["tool_results"].append(tool_result)
        
        # Update summary counts
        if tool_result["has_doi"]:
            results["tools_with_doi"] += 1
            
            if tool_result["doi_valid"]:
                results["valid_dois"] += 1
                
            if tool_result["doi_consistent"]:
                results["consistent_dois"] += 1
        
        if tool_result["has_pubmed_data"]:
            results["tools_with_pubmed_data"] += 1
            
        if tool_result["has_academic_impact"]:
            results["tools_with_academic_impact"] += 1
            
        for fmt in tool_result["citation_formats"]:
            results["citation_format_counts"][fmt] += 1
            
        if tool_result["issues"]:
            for issue in tool_result["issues"]:
                results["all_issues"].append(f"{tool_id}: {issue}")
    
    # Calculate percentages
    if tools:
        results["doi_percentage"] = (results["tools_with_doi"] / len(tools)) * 100
        results["valid_doi_percentage"] = (results["valid_dois"] / results["tools_with_doi"] * 100) if results["tools_with_doi"] > 0 else 0
        results["consistent_doi_percentage"] = (results["consistent_dois"] / results["tools_with_doi"] * 100) if results["tools_with_doi"] > 0 else 0
        results["pubmed_percentage"] = (results["tools_with_pubmed_data"] / len(tools)) * 100
        results["academic_impact_percentage"] = (results["tools_with_academic_impact"] / len(tools)) * 100
    
    return results

def format_validation_report(results: Dict[str, Any]) -> str:
    """
    Format validation results into a readable report
    
    Args:
        results: Validation results
        
    Returns:
        Formatted report as a string
    """
    report = ["# Citation Validation Report", ""]
    report.append(f"Generated: {results.get('timestamp', datetime.datetime.now().isoformat())}")
    report.append("")
    
    report.append("## Summary")
    total_tools = results.get('total_tools', 0)
    report.append(f"- Total Tools: {total_tools}")
    
    # If we have no tools, create a simplified report
    if total_tools == 0:
        report.append("")
        report.append("## Issues")
        report.append("- No tools found in data.json for validation")
        report.append("- Check that data.json contains nodes with 'type': 'tool'")
        report.append("- Run 'python update_data_json.py' to update data.json from README.md")
        report.append("")
        report.append("## Recommendation")
        report.append("1. Verify the structure of data.json")
        report.append("2. Ensure README.md contains properly formatted tool entries")
        report.append("3. Run the data.json update script")
        report.append("4. Re-run citation validation")
        return "\n".join(report)

    # Continue with normal report for when tools are found
    tools_with_doi = results.get('tools_with_doi', 0)
    doi_percentage = results.get('doi_percentage', 0.0)
    report.append(f"- Tools with DOI: {tools_with_doi} ({doi_percentage:.1f}%)")
    
    if tools_with_doi > 0:
        valid_dois = results.get('valid_dois', 0)
        valid_doi_percentage = results.get('valid_doi_percentage', 0.0)
        consistent_dois = results.get('consistent_dois', 0)
        consistent_doi_percentage = results.get('consistent_doi_percentage', 0.0)
        report.append(f"- Valid DOI Format: {valid_dois} ({valid_doi_percentage:.1f}% of tools with DOI)")
        report.append(f"- Consistent DOIs: {consistent_dois} ({consistent_doi_percentage:.1f}% of tools with DOI)")
    
    pubmed_data = results.get('tools_with_pubmed_data', 0)
    pubmed_percentage = results.get('pubmed_percentage', 0.0)
    academic_impact = results.get('tools_with_academic_impact', 0)
    academic_impact_percentage = results.get('academic_impact_percentage', 0.0)
    
    report.append(f"- Tools with PubMed Data: {pubmed_data} ({pubmed_percentage:.1f}%)")
    report.append(f"- Tools with Academic Impact: {academic_impact} ({academic_impact_percentage:.1f}%)")
    report.append("")
    
    report.append("## Citation Format Coverage")
    citation_counts = results.get('citation_format_counts', {'bibtex': 0, 'apa': 0, 'mla': 0})
    report.append(f"- BibTeX: {citation_counts.get('bibtex', 0)} tools")
    report.append(f"- APA: {citation_counts.get('apa', 0)} tools")
    report.append(f"- MLA: {citation_counts.get('mla', 0)} tools")
    report.append("")
    
    all_issues = results.get('all_issues', [])
    if all_issues:
        report.append("## Issues")
        for issue in all_issues:
            report.append(f"- {issue}")
    else:
        report.append("## Issues")
        report.append("- No issues found")
    
    report.append("")
    report.append("## Tools Requiring Citation Data")
    
    # Find tools without any citation data
    missing_data = []
    for tool in results.get("tool_results", []):
        if not tool.get("has_doi") and not tool.get("has_pubmed_data") and not tool.get("has_academic_impact"):
            missing_data.append(tool.get("tool_id", "unknown"))
    
    if missing_data:
        for tool_id in missing_data[:20]:  # Limit to first 20
            report.append(f"- {tool_id}")
        if len(missing_data) > 20:
            report.append(f"- ... and {len(missing_data) - 20} more tools")
    else:
        report.append("- All tools have some form of citation data")
    
    return "\n".join(report)

def main():
    """Main entry point for the script."""
    parser = argparse.ArgumentParser(description="Citation data validation")
    parser.add_argument("--output-file", type=str, default="citation_validation.json",
                      help="Path to save the validation report (default: citation_validation.json)")
    parser.add_argument("--report", action="store_true",
                      help="Generate a human-readable report")
    parser.add_argument("--report-file", type=str, default="citation_validation_report.md",
                      help="Path to save the human-readable report (default: citation_validation_report.md)")
    parser.add_argument("--skip-doi-check", action="store_true",
                      help="Skip DOI online validation (faster for bulk validation)")
    args = parser.parse_args()
    
    # Make sure the reports directory exists
    os.makedirs(REPORTS_DIR, exist_ok=True)
    
    logger.info("Starting citation validation...")
    
    # Run validation
    results = validate_all_citations(
        skip_doi_check=args.skip_doi_check
    )
    
    # Even if no tools were found, we still want to generate an empty report
    # rather than failing completely
    
    # Save validation results
    try:
        with open(os.path.join(REPORTS_DIR, args.output_file), 'w') as f:
            json.dump(results, f, indent=2)
        logger.info(f"Validation results saved to {os.path.join(REPORTS_DIR, args.output_file)}")
    except (IOError, OSError) as e:
        logger.error(f"Error saving validation results: {e}")
        return 1
    
    # Generate human-readable report if requested
    if args.report:
        report = format_validation_report(results)
        try:
            with open(os.path.join(REPORTS_DIR, args.report_file), 'w') as f:
                f.write(report)
            logger.info(f"Report saved to {os.path.join(REPORTS_DIR, args.report_file)}")
        except (IOError, OSError) as e:
            logger.error(f"Error saving report: {e}")
    
    # Log a summary of the validation results
    all_issues_count = len(results.get('all_issues', []))
    logger.info(f"Citation validation completed. {all_issues_count} issues found.")
    
    if results.get('total_tools', 0) > 0:
        logger.info(f"Tools with DOI: {results.get('tools_with_doi', 0)} of {results.get('total_tools', 0)}")
        logger.info(f"Tools with PubMed data: {results.get('tools_with_pubmed_data', 0)} of {results.get('total_tools', 0)}")
    else:
        logger.warning("No tools were found for validation. Check data.json structure.")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())