#!/usr/bin/env python3
"""
DOI scanner for finding DOIs for tools without them.
"""

import json
import logging
import os
import re
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from ..api.crossref_client import CrossRefClient
from ..config import METADATA_DIR, ROOT_DIR
from ..utils.logging_utils import setup_logging
from ..validators.doi_validator import DOIValidator

logger = logging.getLogger(__name__)


class DOIScanner:
    """Scanner for finding DOIs for tools without them."""
    
    def __init__(self):
        """Initialize the DOI scanner."""
        self.doi_validator = DOIValidator()
        self.crossref_client = CrossRefClient()
        self.doi_pattern = re.compile(r'10\.\d{4,9}/[-._;()/:A-Za-z0-9]+')
    
    def scan_readme_for_dois(self) -> Dict[str, str]:
        """
        Scan the README.md file for potential DOIs associated with tools.
        
        Returns:
            Dictionary mapping tool names to potential DOIs
        """
        readme_path = ROOT_DIR / "README.md"
        if not readme_path.exists():
            logger.error("README.md not found")
            return {}
        
        tool_dois = {}
        current_tool = None
        
        with open(readme_path, 'r', encoding='utf-8') as f:
            for line in f:
                # Check for tool headers (# Tool Name)
                if line.startswith('##') and not line.startswith('###'):
                    current_tool = line.strip('#').strip()
                
                # If we're in a tool section, look for DOIs
                if current_tool and 'doi' in line.lower():
                    dois = self.doi_pattern.findall(line)
                    if dois:
                        # Use the first DOI found
                        doi = dois[0]
                        if self.doi_validator.is_valid_format(doi):
                            tool_dois[current_tool] = doi
        
        logger.info(f"Found {len(tool_dois)} DOIs in README")
        return tool_dois
    
    def scan_metadata_for_dois(self) -> Dict[str, str]:
        """
        Scan metadata files for DOIs.
        
        Returns:
            Dictionary mapping tool IDs to DOIs
        """
        tool_dois = {}
        
        # Scan files in the metadata directory
        metadata_files = list(METADATA_DIR.glob("*.json"))
        metadata_files.extend(METADATA_DIR.glob("*/*.json"))
        
        for file_path in metadata_files:
            try:
                with open(file_path, 'r', encoding='utf-8') as f:
                    data = json.load(f)
                
                # Get tool ID from filename
                tool_id = file_path.stem
                
                # Check if the file has a DOI
                if isinstance(data, dict) and 'doi' in data:
                    doi = data['doi']
                    if doi and self.doi_validator.is_valid_format(doi):
                        tool_dois[tool_id] = doi
            except Exception as e:
                logger.warning(f"Error reading {file_path}: {e}")
        
        logger.info(f"Found {len(tool_dois)} DOIs in metadata files")
        return tool_dois
    
    def search_for_doi_by_name(self, tool_name: str) -> Optional[str]:
        """
        Search for a DOI based on tool name.
        
        Args:
            tool_name: Name of the tool
        
        Returns:
            DOI if found, None otherwise
        """
        logger.info(f"Searching for DOI for tool: {tool_name}")
        
        # Try to find the DOI using CrossRef
        success, results = self.crossref_client.search_by_title(tool_name)
        
        if not success or not results:
            logger.info(f"No results found for tool: {tool_name}")
            return None
        
        # Check the results for the most relevant match
        best_match = None
        best_score = 0
        
        for result in results:
            score = result.get('score', 0)
            doi = result.get('doi')
            
            # Check if the DOI is valid and score is better than previous best
            if doi and self.doi_validator.is_valid_format(doi) and score > best_score:
                best_match = doi
                best_score = score
        
        if best_match:
            logger.info(f"Found potential DOI for {tool_name}: {best_match}")
            return best_match
        
        return None
    
    def find_missing_dois(self) -> Dict[str, Dict[str, str]]:
        """
        Find tools without DOIs and search for potential DOIs.
        
        Returns:
            Dictionary mapping tool IDs to potential DOI info
        """
        # Get all tools with DOIs
        known_dois = self.scan_metadata_for_dois()
        readme_dois = self.scan_readme_for_dois()
        
        # Merge DOIs from different sources
        all_tools_with_dois = {**known_dois, **readme_dois}
        
        # Get all tools from the data.json file
        data_path = ROOT_DIR / "data.json"
        all_tools = set()
        
        try:
            with open(data_path, 'r', encoding='utf-8') as f:
                data = json.load(f)
                
                # Check if data has a tools key or nodes key
                if 'tools' in data:
                    tools = data['tools']
                elif 'nodes' in data:
                    # Filter for nodes of type 'tool'
                    tools = [node for node in data.get('nodes', []) 
                            if node.get('type') == 'tool']
                else:
                    tools = []
                
                for tool in tools:
                    tool_id = tool.get('id')
                    if tool_id:
                        all_tools.add(tool_id)
        except Exception as e:
            logger.error(f"Error reading data.json: {e}")
        
        logger.info(f"Found {len(all_tools)} total tools")
        logger.info(f"Found {len(all_tools_with_dois)} tools with DOIs")
        
        # Find tools without DOIs
        tools_without_dois = all_tools - set(all_tools_with_dois.keys())
        logger.info(f"Found {len(tools_without_dois)} tools without DOIs")
        
        # Search for DOIs for each tool without one
        potential_dois = {}
        
        for tool_id in tools_without_dois:
            # Clean up tool_id to get a searchable name
            tool_name = tool_id
            if tool_name.startswith("tool-"):
                tool_name = tool_name[5:]
            
            # Replace underscores and dashes with spaces
            tool_name = tool_name.replace('_', ' ').replace('-', ' ')
            
            doi = self.search_for_doi_by_name(tool_name)
            
            if doi:
                potential_dois[tool_id] = {
                    'tool_id': tool_id,
                    'tool_name': tool_name,
                    'potential_doi': doi,
                    'confidence': 'medium'  # Default confidence
                }
        
        logger.info(f"Found {len(potential_dois)} potential DOIs for tools without DOIs")
        return potential_dois
    
    def export_potential_dois(self, output_path: Optional[Path] = None) -> Path:
        """
        Export potential DOIs to a file.
        
        Args:
            output_path: Path to write the output file (defaults to ROOT_DIR/potential_dois.json)
        
        Returns:
            Path to the output file
        """
        if output_path is None:
            output_path = ROOT_DIR / "potential_dois.json"
        
        potential_dois = self.find_missing_dois()
        
        # Write the potential DOIs to a file
        with open(output_path, 'w', encoding='utf-8') as f:
            json.dump({
                'total': len(potential_dois),
                'tools': list(potential_dois.values())
            }, f, indent=2)
        
        logger.info(f"Exported {len(potential_dois)} potential DOIs to {output_path}")
        return output_path


def main():
    """Command-line interface for the DOI scanner."""
    import argparse
    import sys
    
    # Set up logging
    logger = setup_logging("doi_scanner")
    
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Scan for and find missing DOIs")
    parser.add_argument(
        "-o", "--output", 
        help="Output file for potential DOIs (JSON format)"
    )
    
    args = parser.parse_args()
    
    # Create the scanner and find missing DOIs
    scanner = DOIScanner()
    
    output_path = Path(args.output) if args.output else None
    output_file = scanner.export_potential_dois(output_path)
    
    logger.info(f"DOI scanning complete. Results saved to {output_file}")
    return 0


if __name__ == "__main__":
    sys.exit(main())