#!/usr/bin/env python3
"""
Generate API endpoints from the main data.json file

This script creates a set of API endpoints for the awesome-virome repository,
allowing programmatic access to the database of virome analysis tools.
"""

import os
import sys
import json
import time
import shutil
import datetime
import logging
from pathlib import Path
from collections import defaultdict

# Configuration
DATA_FILE = "data.json"
API_DIR = "api/v1"
API_VERSION = "1.0"

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Category slug mapping
CATEGORY_MAP = {
    "Virus and Phage Identification": "virus-identification",
    "Host Prediction": "host-prediction",
    "Genome Analysis": "genome-analysis",
    "Genome Annotation": "genome-annotation",
    "Genome Assembly": "genome-assembly",
    "Genome Completeness": "genome-completeness",
    "Genome Comparison": "genome-comparison",
    "Gene Finding": "gene-finding",
    "Taxonomy": "taxonomy",
    "Quality Control": "quality-control",
    "Databases": "databases",
    "Sequence Databases": "sequence-databases",
    "Functional Analysis": "functional-analysis",
    "Evolutionary Analysis": "evolutionary-analysis",
    "Lifestyle Classification": "lifestyle-classification",
    "Phage-specific Analysis": "phage-specific-analysis",
    "Viral Orthologous Groups": "viral-orthologous-groups",
    "CRISPR Analysis": "crispr-analysis",
    "Sequence Analysis": "sequence-analysis",
    "Multiple Sequence Alignment": "multiple-sequence-alignment",
    "Sequence Translation": "sequence-translation",
    "Viral Strain Reconstruction": "viral-strain-reconstruction",
    "Viral Quasispecies Analysis": "viral-quasispecies-analysis",
    "Visualization and Infrastructure": "visualization",
    "Cyberinfrastructure": "cyberinfrastructure",
    "Plaque Analysis": "plaque-analysis",
    "Other Tools": "other-tools",
    "Machine Learning Models": "machine-learning",
    "Structural Analysis Tools": "structural-analysis",
    "Antimicrobial Resistance Analysis": "amr-analysis",
    "Viral Metatranscriptomics": "metatranscriptomics",
    "Cloud-based Viral Analysis": "cloud-based-analysis",
    "Simulation": "simulation",
    "Amplicon Analysis": "amplicon-analysis",
    "Interaction Analysis": "interaction-analysis",
    "Viral Single-Cell Analysis": "single-cell-analysis",
    "Viral Glycoprotein Analysis": "glycoprotein-analysis",
    "Ancient Viral Sequence Analysis": "ancient-viral-analysis",
    "Viral Immune Epitope Prediction": "immune-epitope-prediction",
    "Viral Molecular Dynamics": "molecular-dynamics",
    "Dark Matter Viral Analysis": "dark-matter-analysis",
    "Transduction": "transduction",
    "Metagenome Analysis": "metagenome-analysis",
    "Integrated Viruses": "integrated-viruses",
    "RNA Virus Identification": "rna-virus-identification",
}

# Map slug back to pretty name
CATEGORY_NAMES = {v: k for k, v in CATEGORY_MAP.items()}

def slugify(name):
    """
    Convert a name to a URL-friendly slug.
    
    Args:
        name: The name to convert
        
    Returns:
        A URL-friendly slug
    """
    if name in CATEGORY_MAP:
        return CATEGORY_MAP[name]
    
    # Generic fallback slugify
    return name.lower().replace(" ", "-").replace("/", "-")

def create_dirs():
    """Create the API directory structure."""
    api_path = Path(API_DIR)
    api_path.mkdir(parents=True, exist_ok=True)
    
    cat_path = api_path / "categories"
    cat_path.mkdir(exist_ok=True)
    
    return api_path, cat_path

def load_data():
    """
    Load the main data file.
    
    Returns:
        The parsed data JSON
    """
    try:
        with open(DATA_FILE, 'r') as f:
            return json.load(f)
    except (json.JSONDecodeError, IOError) as e:
        logger.error(f"Failed to load data file {DATA_FILE}: {e}")
        sys.exit(1)

def extract_tools(data):
    """
    Extract tools from the data structure.
    
    Args:
        data: The parsed data JSON
        
    Returns:
        List of tool objects
    """
    tools = []
    nodes = data.get('nodes', [])
    
    for node in nodes:
        if node.get('type') == 'tool':
            # Clean up tool object for API
            tool = {
                "id": node.get('id', '').replace('tool-', ''),
                "name": node.get('name', ''),
                "description": node.get('description', ''),
                "url": node.get('url', ''),
                "category": node.get('category', ''),
                "subcategory": node.get('subcategory', ''),
                "language": node.get('language'),
                "stars": node.get('stars', 0),
                "forks": node.get('forks', 0),
                "license": node.get('license'),
                "topics": node.get('topics', []),
                "updated_at": node.get('lastUpdated'),
                "created_at": node.get('createdAt'),
                "languages": node.get('languages', {}),
                "doi": node.get('doi'),
                "citation_count": node.get('citation_count', 0),
                "package_manager": node.get('package_manager')
            }
            tools.append(tool)
    
    return tools

def create_tools_endpoint(api_path, tools):
    """
    Create the main tools.json endpoint.
    
    Args:
        api_path: Path to the API directory
        tools: List of tool objects
    """
    endpoint = {
        "api_version": API_VERSION,
        "generated_at": datetime.datetime.now().isoformat(),
        "count": len(tools),
        "tools": tools
    }
    
    with open(api_path / "tools.json", 'w') as f:
        json.dump(endpoint, f, indent=2)
    
    logger.info(f"Created tools endpoint with {len(tools)} tools")

def create_category_endpoints(api_path, cat_path, tools):
    """
    Create category-specific endpoints.
    
    Args:
        api_path: Path to the API directory
        cat_path: Path to the categories directory
        tools: List of tool objects
    """
    # Group tools by category
    categories = defaultdict(list)
    category_counts = defaultdict(int)
    
    for tool in tools:
        category = tool.get('category')
        if category:
            slug = slugify(category)
            categories[slug].append(tool)
            category_counts[slug] += 1
        
        # Also add to subcategory if present
        subcategory = tool.get('subcategory')
        if subcategory:
            subslug = slugify(subcategory)
            categories[subslug].append(tool)
            category_counts[subslug] += 1
    
    # Create individual category files
    for slug, cat_tools in categories.items():
        category_name = CATEGORY_NAMES.get(slug, slug.replace("-", " ").title())
        
        endpoint = {
            "api_version": API_VERSION,
            "generated_at": datetime.datetime.now().isoformat(),
            "category": category_name,
            "count": len(cat_tools),
            "tools": cat_tools
        }
        
        with open(cat_path / f"{slug}.json", 'w') as f:
            json.dump(endpoint, f, indent=2)
    
    # Create categories index
    categories_list = []
    for slug, count in category_counts.items():
        category_name = CATEGORY_NAMES.get(slug, slug.replace("-", " ").title())
        categories_list.append({
            "slug": slug,
            "name": category_name,
            "count": count,
            "endpoint": f"{API_DIR}/categories/{slug}.json"
        })
    
    # Sort by name
    categories_list.sort(key=lambda x: x["name"])
    
    endpoint = {
        "api_version": API_VERSION,
        "generated_at": datetime.datetime.now().isoformat(),
        "count": len(categories_list),
        "categories": categories_list
    }
    
    with open(api_path / "categories.json", 'w') as f:
        json.dump(endpoint, f, indent=2)
    
    logger.info(f"Created category endpoints for {len(categories)} categories")

def create_search_index(api_path, tools):
    """
    Create a search index for client-side filtering.
    
    Args:
        api_path: Path to the API directory
        tools: List of tool objects
    """
    # Create lightweight search index
    search_index = []
    for tool in tools:
        index_entry = {
            "id": tool.get('id'),
            "name": tool.get('name'),
            "description": tool.get('description'),
            "category": tool.get('category'),
            "subcategory": tool.get('subcategory'),
            "language": tool.get('language'),
            "topics": tool.get('topics', []),
            "stars": tool.get('stars', 0),
            "url": tool.get('url')
        }
        search_index.append(index_entry)
    
    endpoint = {
        "api_version": API_VERSION,
        "generated_at": datetime.datetime.now().isoformat(),
        "count": len(search_index),
        "tools": search_index
    }
    
    with open(api_path / "search.json", 'w') as f:
        json.dump(endpoint, f, indent=2)
    
    logger.info(f"Created search index with {len(search_index)} entries")

def create_stats_endpoint(api_path, tools):
    """
    Create statistics endpoint.
    
    Args:
        api_path: Path to the API directory
        tools: List of tool objects
    """
    # Calculate statistics
    total_tools = len(tools)
    
    # Stars
    total_stars = sum(tool.get('stars', 0) or 0 for tool in tools)
    avg_stars = total_stars / total_tools if total_tools > 0 else 0
    
    # Languages
    languages = defaultdict(int)
    for tool in tools:
        tool_langs = tool.get('languages', {}) or {}
        for lang in tool_langs:
            languages[lang] += 1
    
    # Package managers
    package_managers = defaultdict(int)
    for tool in tools:
        pm = tool.get('package_manager')
        if pm:
            package_managers[pm] += 1
    
    # Category counts
    categories = defaultdict(int)
    for tool in tools:
        category = tool.get('category')
        if category:
            categories[category] += 1
    
    # Create stats endpoint
    stats = {
        "api_version": API_VERSION,
        "generated_at": datetime.datetime.now().isoformat(),
        "total_tools": total_tools,
        "total_stars": total_stars,
        "average_stars": round(avg_stars, 2),
        "languages": dict(sorted(languages.items(), key=lambda x: x[1], reverse=True)),
        "package_managers": dict(sorted(package_managers.items(), key=lambda x: x[1], reverse=True)),
        "categories": dict(sorted(categories.items(), key=lambda x: x[1], reverse=True))
    }
    
    with open(api_path / "stats.json", 'w') as f:
        json.dump(stats, f, indent=2)
    
    logger.info("Created statistics endpoint")

def create_api_metadata(api_path):
    """
    Create API metadata endpoint.
    
    Args:
        api_path: Path to the API directory
    """
    metadata = {
        "api_name": "Awesome Virome API",
        "api_version": API_VERSION,
        "generated_at": datetime.datetime.now().isoformat(),
        "description": "API for accessing the awesome-virome database of virus and phage analysis tools",
        "documentation_url": "https://github.com/shandley/awesome-virome/blob/main/API.md",
        "endpoints": [
            {
                "path": f"/{API_DIR}/tools.json",
                "description": "Complete list of tools",
                "methods": ["GET"]
            },
            {
                "path": f"/{API_DIR}/categories.json",
                "description": "List of available categories",
                "methods": ["GET"]
            },
            {
                "path": f"/{API_DIR}/categories/{{category_name}}.json",
                "description": "Tools filtered by category",
                "methods": ["GET"],
                "parameters": {
                    "category_name": "The category slug (e.g., 'virus-identification')"
                }
            },
            {
                "path": f"/{API_DIR}/search.json",
                "description": "Search index for client-side filtering",
                "methods": ["GET"]
            },
            {
                "path": f"/{API_DIR}/stats.json",
                "description": "Statistics about the tools database",
                "methods": ["GET"]
            },
            {
                "path": f"/{API_DIR}/metadata.json",
                "description": "API metadata and documentation",
                "methods": ["GET"]
            }
        ],
        "license": "CC0",
        "repository": "https://github.com/shandley/awesome-virome"
    }
    
    with open(api_path / "metadata.json", 'w') as f:
        json.dump(metadata, f, indent=2)
    
    logger.info("Created API metadata endpoint")

def main():
    """Main entry point for the script."""
    start_time = time.time()
    logger.info("Starting API generation")
    
    # Create directories
    api_path, cat_path = create_dirs()
    
    # Load data
    data = load_data()
    
    # Extract tools
    tools = extract_tools(data)
    
    # Create endpoints
    create_tools_endpoint(api_path, tools)
    create_category_endpoints(api_path, cat_path, tools)
    create_search_index(api_path, tools)
    create_stats_endpoint(api_path, tools)
    create_api_metadata(api_path)
    
    elapsed_time = time.time() - start_time
    logger.info(f"API generation completed in {elapsed_time:.2f} seconds")

if __name__ == "__main__":
    main()