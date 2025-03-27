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
import apis.biotools_api as biotools_api
import apis.bioconda_api as bioconda_api

# Import academic impact module
from academic_impact import AcademicImpactCollector

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

# Directory constants (declare at module level, no global needed)
METADATA_DIR = os.path.join("metadata", "bioinformatics")
OUTPUT_FILE = os.path.join(METADATA_DIR, "summary.json")


class BioinformaticsMetadataCollector:
    """Collects specialized bioinformatics metadata for Awesome-Virome tools."""
    
    def __init__(self, github_token: Optional[str] = None, 
                 semantic_scholar_key: Optional[str] = None,
                 contact_email: Optional[str] = None,
                 metadata_dir: str = METADATA_DIR,
                 output_file: str = OUTPUT_FILE):
        """Initialize the bioinformatics metadata collector."""
        # Store directory paths
        self.metadata_dir = metadata_dir
        self.output_file = output_file
        
        # Initialize academic impact collector
        self.academic_impact_collector = AcademicImpactCollector(
            github_token=github_token,
            semantic_scholar_key=semantic_scholar_key,
            contact_email=contact_email
        )
        
        # Create metadata directory if it doesn't exist
        os.makedirs(self.metadata_dir, exist_ok=True)
        
        # Load existing metadata if available
        self.existing_metadata = {}
        if os.path.exists(self.output_file):
            try:
                with open(self.output_file, 'r') as f:
                    self.existing_metadata = json.load(f)
            except json.JSONDecodeError:
                logger.warning(f"Could not parse existing metadata from {self.output_file}")
    
    def search_biotools(self, tool_name: str) -> Dict[str, Any]:
        """Search Bio.tools for a tool by name and return its metadata."""
        # Try different search terms
        search_terms = [
            tool_name,
            tool_name.replace('-', ' '),
            tool_name.lower(),
            tool_name.replace('_', ' ')
        ]
        
        try:
            for term in search_terms:
                results = biotools_api.search_tool(term)
                
                if results and 'list' in results:
                    result_list = results.get('list', [])
                    # Try to find an exact match first
                    for result in result_list:
                        name = result.get('name', '').lower()
                        if name == tool_name.lower() or name == tool_name.lower().replace('-', ''):
                            tool_id = result.get('biotoolsID')
                            tool_details = biotools_api.get_tool_details(tool_id)
                            return tool_details if tool_details else {}
                    
                    # If no exact match, return the first result
                    if result_list:
                        tool_id = result_list[0].get('biotoolsID')
                        tool_details = biotools_api.get_tool_details(tool_id)
                        return tool_details if tool_details else {}
        except Exception as e:
            logger.error(f"Error searching Bio.tools for {tool_name}: {e}")
        
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
        
        try:
            for term in search_terms:
                package_data = bioconda_api.search_package(term)
                
                if package_data:
                    # Get additional information
                    package_name = package_data.get('name', '')
                    files_data = bioconda_api.get_package_files(package_name)
                    recipe_data = bioconda_api.get_package_recipe(package_name)
                    
                    # Extract structured metadata
                    metadata = bioconda_api.extract_package_metadata(package_data, files_data, recipe_data)
                    return metadata if metadata else {}
        except Exception as e:
            logger.error(f"Error searching Bioconda for {tool_name}: {e}")
        
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
        if bioconda_data is None:
            bioconda_data = {}
        
        # Collect academic impact data
        academic_impact_data = self.academic_impact_collector.process_tool(tool)
        if academic_impact_data is None:
            academic_impact_data = {}
        
        # Ensure biotools_data is not None
        if biotools_data is None:
            biotools_data = {}
            
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
        
        # Note: File saving is handled by the subclass
        
        return metadata
    
    def _extract_data_formats(self, tool_data: Dict[str, Any], io_type: str) -> List[str]:
        """Extract input or output data formats from Bio.tools data."""
        formats = []
        
        # Ensure tool_data is a dict and not None
        if not tool_data or not isinstance(tool_data, dict):
            return formats
            
        for function in tool_data.get('function', []):
            if not function or not isinstance(function, dict):
                continue
                
            for io in function.get('operation', []):
                if not io or not isinstance(io, dict):
                    continue
                    
                for data in io.get(f'{io_type}', []):
                    if not data or not isinstance(data, dict):
                        continue
                        
                    data_format = data.get('data', {}).get('format', [])
                    if data_format and isinstance(data_format, list):
                        formats.extend([fmt.get('term', '') for fmt in data_format if isinstance(fmt, dict) and fmt.get('term')])
        
        return list(set(formats))
    
    def _extract_operations(self, tool_data: Dict[str, Any]) -> List[str]:
        """Extract bioinformatics operations from Bio.tools data."""
        operations = []
        
        # Ensure tool_data is a dict and not None
        if not tool_data or not isinstance(tool_data, dict):
            return operations
            
        for function in tool_data.get('function', []):
            if not function or not isinstance(function, dict):
                continue
                
            for operation in function.get('operation', []):
                if not operation or not isinstance(operation, dict):
                    continue
                    
                term = operation.get('term', '')
                if term:
                    operations.append(term)
        
        return list(set(operations))
    
    def _extract_topics(self, tool_data: Dict[str, Any]) -> List[str]:
        """Extract bioinformatics topics from Bio.tools data."""
        topics = []
        
        # Ensure tool_data is a dict and not None
        if not tool_data or not isinstance(tool_data, dict):
            return topics
            
        for topic in tool_data.get('topic', []):
            if not topic or not isinstance(topic, dict):
                continue
                
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
                
                # Save individual tool metadata
                tool_file = os.path.join(self.metadata_dir, f"{tool_name.replace('/', '_')}.json")
                with open(tool_file, 'w') as f:
                    json.dump(metadata, f, indent=2)
        
        # Save all results
        with open(self.output_file, 'w') as f:
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
        
        # Save summary statistics
        summary_file = os.path.join(self.metadata_dir, "statistics.json")
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
    
    # Set up directory paths
    metadata_dir = args.output
    output_file = os.path.join(metadata_dir, "summary.json")
    
    # Ensure metadata directory exists
    os.makedirs(metadata_dir, exist_ok=True)
    
    # Load tools data
    tools = load_tools_data()
    if not tools:
        logger.error("No tools found in data.json")
        return
    
    logger.info(f"Loaded {len(tools)} tools from data.json")
    
    # Initialize metadata collector with our directory paths
    collector = BioinformaticsMetadataCollector(
        github_token=args.github_token,
        semantic_scholar_key=args.semantic_scholar_key,
        contact_email=args.contact_email,
        metadata_dir=metadata_dir,
        output_file=output_file
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