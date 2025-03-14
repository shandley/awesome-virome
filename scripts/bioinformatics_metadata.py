#!/usr/bin/env python3
"""
Bioinformatics Metadata Collection for Awesome-Virome

This script collects specialized bioinformatics metadata for tools in the
Awesome-Virome repository, including:
- Input/output formats from Bio.tools
- Bioinformatics operations and topics from Bio.tools
- Package information and dependencies from Bioconda
- Academic impact and citation information

Usage:
    python bioinformatics_metadata.py [--output OUTPUT] [--token TOKEN]
"""

import os
import sys
import json
import argparse
import logging
from datetime import datetime
from typing import Dict, List, Optional, Any, Tuple
from concurrent.futures import ThreadPoolExecutor
import re

# Import our API modules
from apis.biotools_api import BioToolsAPI
from apis.bioconda_api import BiocondaAPI

# Import academic impact module
from scripts.academic_impact import AcademicImpactCollector

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("metadata_collection.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Directory constants
METADATA_DIR = os.path.join("metadata", "bioinformatics")
OUTPUT_FILE = os.path.join(METADATA_DIR, "summary.json")


class BioinformaticsMetadataCollector:
    """Collects specialized bioinformatics metadata for Awesome-Virome tools."""
    
    def __init__(self, github_token: Optional[str] = None, 
                 semantic_scholar_key: Optional[str] = None,
                 contact_email: Optional[str] = None):
        """Initialize the bioinformatics metadata collector."""
        self.biotools_api = BioToolsAPI()
        self.bioconda_api = BiocondaAPI()
        
        # Initialize academic impact collector
        self.academic_impact_collector = AcademicImpactCollector(
            github_token=github_token,
            semantic_scholar_key=semantic_scholar_key,
            contact_email=contact_email
        )
        
        # Create metadata directory if it doesn't exist
        os.makedirs(METADATA_DIR, exist_ok=True)
        
        # Load existing metadata if available
        self.existing_metadata = {}
        if os.path.exists(OUTPUT_FILE):
            try:
                with open(OUTPUT_FILE, 'r') as f:
                    self.existing_metadata = json.load(f)
            except json.JSONDecodeError:
                logger.warning(f"Could not parse existing metadata from {OUTPUT_FILE}")
    
    def search_biotools(self, tool_name: str) -> Dict[str, Any]:
        """Search Bio.tools for a tool by name and return its metadata."""
        # Try different search terms
        search_terms = [
            tool_name,
            tool_name.replace('-', ' '),
            tool_name.lower(),
            tool_name.replace('_', ' ')
        ]
        
        for term in search_terms:
            results = self.biotools_api.search_tools(term)
            
            if results:
                # Try to find an exact match first
                for result in results:
                    name = result.get('name', '').lower()
                    if name == tool_name.lower() or name == tool_name.lower().replace('-', ''):
                        tool_id = result.get('biotoolsID')
                        return self.biotools_api.get_tool_details(tool_id)
                
                # If no exact match, return the first result
                tool_id = results[0].get('biotoolsID')
                return self.biotools_api.get_tool_details(tool_id)
        
        return {}
    
    def search_bioconda(self, tool_name: str) -> Dict[str, Any]:
        """Search Bioconda for a package by name and return its metadata."""
        # Try different search terms
        search_terms = [
            tool_name,
            tool_name.replace('-', ''),
            tool_name.lower(),
            f"bioconda/{tool_name}",
            f"bioconda-{tool_name}"
        ]
        
        for term in search_terms:
            results = self.bioconda_api.search_package(term)
            
            if results:
                # Try to find an exact match first
                for result in results:
                    if result.lower() == tool_name.lower() or \
                       result.lower() == f"bioconda-{tool_name.lower()}" or \
                       result.lower() == tool_name.lower().replace('-', ''):
                        return self.bioconda_api.get_package_info(result)
                
                # If no exact match, return the first result
                return self.bioconda_api.get_package_info(results[0])
        
        return {}
    
    def collect_tool_metadata(self, tool: Dict[str, Any]) -> Dict[str, Any]:
        """Collect bioinformatics metadata for a single tool."""
        tool_name = tool.get('name', '')
        tool_url = tool.get('url', '')
        
        logger.info(f"Collecting metadata for tool: {tool_name}")
        
        # Check if we already have metadata for this tool
        if tool_name in self.existing_metadata:
            last_updated = self.existing_metadata[tool_name].get('last_updated', '')
            # If updated within the last 30 days, use cached data
            if last_updated and (datetime.now() - datetime.fromisoformat(last_updated)).days < 30:
                logger.info(f"Using cached metadata for {tool_name}")
                return self.existing_metadata[tool_name]
        
        # Collect data from Bio.tools
        biotools_data = self.search_biotools(tool_name)
        
        # Collect data from Bioconda
        bioconda_data = self.search_bioconda(tool_name)
        
        # Collect academic impact data
        academic_impact_data = self.academic_impact_collector.process_tool(tool)
        
        # Compile all metadata
        metadata = {
            'name': tool_name,
            'url': tool_url,
            'biotools': {
                'id': biotools_data.get('biotoolsID', ''),
                'name': biotools_data.get('name', ''),
                'description': biotools_data.get('description', ''),
                'homepage': biotools_data.get('homepage', ''),
                'input_formats': self._extract_data_formats(biotools_data, 'input'),
                'output_formats': self._extract_data_formats(biotools_data, 'output'),
                'operations': self._extract_operations(biotools_data),
                'topics': self._extract_topics(biotools_data)
            },
            'bioconda': {
                'package': bioconda_data.get('package', ''),
                'version': bioconda_data.get('version', ''),
                'dependencies': bioconda_data.get('dependencies', []),
                'installation': self._format_installation_command(bioconda_data.get('package', ''))
            },
            'academic_impact': {
                'doi': academic_impact_data.get('doi', ''),
                'citation_file': academic_impact_data.get('citation_info', {}).get('citation_file', ''),
                'citation_format': academic_impact_data.get('citation_info', {}).get('citation_format', ''),
                'total_citations': academic_impact_data.get('citation_metrics', {}).get('metrics', {}).get('total_citations', 0),
                'influential_citations': academic_impact_data.get('citation_metrics', {}).get('metrics', {}).get('influential_citations', 0),
                'citations_by_year': academic_impact_data.get('citation_metrics', {}).get('metrics', {}).get('citations_by_year', {}),
                'formatted_citations': academic_impact_data.get('citation_metrics', {}).get('formatted_citations', {})
            },
            'last_updated': datetime.now().isoformat()
        }
        
        # Save individual tool metadata
        tool_file = os.path.join(METADATA_DIR, f"{tool_name.replace('/', '_')}.json")
        with open(tool_file, 'w') as f:
            json.dump(metadata, f, indent=2)
        
        return metadata
    
    def _extract_data_formats(self, tool_data: Dict[str, Any], io_type: str) -> List[str]:
        """Extract input or output data formats from Bio.tools data."""
        formats = []
        
        for function in tool_data.get('function', []):
            for io in function.get('operation', []):
                for data in io.get(f'{io_type}', []):
                    data_format = data.get('data', {}).get('format', [])
                    formats.extend([fmt.get('term', '') for fmt in data_format if fmt.get('term')])
        
        return list(set(formats))
    
    def _extract_operations(self, tool_data: Dict[str, Any]) -> List[str]:
        """Extract bioinformatics operations from Bio.tools data."""
        operations = []
        
        for function in tool_data.get('function', []):
            for operation in function.get('operation', []):
                term = operation.get('term', '')
                if term:
                    operations.append(term)
        
        return list(set(operations))
    
    def _extract_topics(self, tool_data: Dict[str, Any]) -> List[str]:
        """Extract bioinformatics topics from Bio.tools data."""
        topics = []
        
        for topic in tool_data.get('topic', []):
            term = topic.get('term', '')
            if term:
                topics.append(term)
        
        return list(set(topics))
    
    def _format_installation_command(self, package_name: str) -> str:
        """Format a conda installation command for the package."""
        if not package_name:
            return ""
        
        return f"conda install -c bioconda {package_name}"
    
    def collect_all_metadata(self, tools: List[Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
        """Collect bioinformatics metadata for all tools."""
        results = {}
        
        # Process tools in parallel using ThreadPoolExecutor
        with ThreadPoolExecutor(max_workers=5) as executor:
            # Map each tool to its metadata
            for tool, metadata in zip(tools, executor.map(self.collect_tool_metadata, tools)):
                tool_name = tool.get('name', '')
                results[tool_name] = metadata
        
        # Save all results
        with open(OUTPUT_FILE, 'w') as f:
            json.dump(results, f, indent=2)
        
        return results
    
    def generate_summary(self, metadata: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
        """Generate a summary of the collected metadata."""
        total_tools = len(metadata)
        
        # Count tools with different types of metadata
        tools_with_biotools = sum(1 for data in metadata.values() if data.get('biotools', {}).get('id'))
        tools_with_bioconda = sum(1 for data in metadata.values() if data.get('bioconda', {}).get('package'))
        tools_with_doi = sum(1 for data in metadata.values() if data.get('academic_impact', {}).get('doi'))
        tools_with_citations = sum(1 for data in metadata.values() 
                              if data.get('academic_impact', {}).get('total_citations', 0) > 0)
        
        # Collect input/output formats
        all_input_formats = []
        all_output_formats = []
        
        for data in metadata.values():
            all_input_formats.extend(data.get('biotools', {}).get('input_formats', []))
            all_output_formats.extend(data.get('biotools', {}).get('output_formats', []))
        
        # Count format frequencies
        input_format_counts = {}
        for fmt in all_input_formats:
            input_format_counts[fmt] = input_format_counts.get(fmt, 0) + 1
        
        output_format_counts = {}
        for fmt in all_output_formats:
            output_format_counts[fmt] = output_format_counts.get(fmt, 0) + 1
        
        # Get most common formats
        most_common_input = sorted(input_format_counts.items(), key=lambda x: x[1], reverse=True)[:10]
        most_common_output = sorted(output_format_counts.items(), key=lambda x: x[1], reverse=True)[:10]
        
        # Collect citation statistics
        total_citations = sum(data.get('academic_impact', {}).get('total_citations', 0) for data in metadata.values())
        most_cited = sorted(
            [(name, data.get('academic_impact', {}).get('total_citations', 0)) 
             for name, data in metadata.items() if data.get('academic_impact', {}).get('total_citations', 0) > 0],
            key=lambda x: x[1],
            reverse=True
        )[:10]
        
        summary = {
            'total_tools': total_tools,
            'tools_with_biotools': tools_with_biotools,
            'tools_with_bioconda': tools_with_bioconda,
            'tools_with_doi': tools_with_doi,
            'tools_with_citations': tools_with_citations,
            'total_citations': total_citations,
            'most_cited_tools': most_cited,
            'most_common_input_formats': most_common_input,
            'most_common_output_formats': most_common_output,
            'generated': datetime.now().isoformat()
        }
        
        summary_file = os.path.join(METADATA_DIR, "statistics.json")
        with open(summary_file, 'w') as f:
            json.dump(summary, f, indent=2)
        
        return summary


def load_tools_data() -> List[Dict[str, Any]]:
    """Load tools data from data.json."""
    try:
        with open('data.json', 'r') as f:
            data = json.load(f)
            return data.get('nodes', [])
    except (json.JSONDecodeError, FileNotFoundError) as e:
        logger.error(f"Error loading data.json: {e}")
        return []


def main():
    """Main function to run the bioinformatics metadata collection."""
    parser = argparse.ArgumentParser(description='Collect bioinformatics metadata for Awesome-Virome tools')
    parser.add_argument('--output', help='Output directory for metadata', default=METADATA_DIR)
    parser.add_argument('--github-token', help='GitHub API token')
    parser.add_argument('--semantic-scholar-key', help='Semantic Scholar API key')
    parser.add_argument('--contact-email', help='Contact email for API rate limiting')
    args = parser.parse_args()
    
    # Update output directory if specified
    global METADATA_DIR, OUTPUT_FILE
    if args.output != METADATA_DIR:
        METADATA_DIR = args.output
        OUTPUT_FILE = os.path.join(METADATA_DIR, "summary.json")
    
    # Ensure metadata directory exists
    os.makedirs(METADATA_DIR, exist_ok=True)
    
    # Load tools data
    tools = load_tools_data()
    if not tools:
        logger.error("No tools found in data.json")
        return
    
    logger.info(f"Loaded {len(tools)} tools from data.json")
    
    # Initialize metadata collector
    collector = BioinformaticsMetadataCollector(
        github_token=args.github_token,
        semantic_scholar_key=args.semantic_scholar_key,
        contact_email=args.contact_email
    )
    
    # Collect metadata for all tools
    metadata = collector.collect_all_metadata(tools)
    
    # Generate summary statistics
    summary = collector.generate_summary(metadata)
    
    logger.info(f"Metadata collection complete. Processed {len(metadata)} tools.")
    logger.info(f"Found Bio.tools entries for {summary['tools_with_biotools']} tools.")
    logger.info(f"Found Bioconda packages for {summary['tools_with_bioconda']} tools.")
    logger.info(f"Found DOIs for {summary['tools_with_doi']} tools.")
    logger.info(f"Found citations for {summary['tools_with_citations']} tools.")
    logger.info(f"Total citations: {summary['total_citations']}")


if __name__ == "__main__":
    main()