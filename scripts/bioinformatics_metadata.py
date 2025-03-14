#!/usr/bin/env python3
"""
Bioinformatics metadata collection script for the Awesome-Virome repository.

This script orchestrates the collection of specialized bioinformatics metadata
for tools in the repository, using various API integrations:
- Bio.tools: A comprehensive registry of bioinformatics tools
- Bioconda: A channel for the conda package manager specializing in bioinformatics
- GitHub: For repository-specific metadata

The collected data enhances the repository with:
- Input/output data formats supported by each tool
- Bioinformatics categories and functions
- Installation methods (conda/bioconda)
- Tool dependencies
- Compatible workflows
"""

import os
import re
import json
import time
import logging
import argparse
import requests
from pathlib import Path
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed

# Import API modules
from apis import biotools_api, bioconda_api

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("bioinformatics_metadata.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Set up constants
REPO_ROOT = Path(__file__).parent.parent
BIOINFO_METADATA_DIR = REPO_ROOT / "metadata" / "bioinformatics"
BIOINFO_METADATA_SUMMARY = BIOINFO_METADATA_DIR / "summary.json"
README_PATH = REPO_ROOT / "README.md"
REPO_UPDATES_JSON = REPO_ROOT / "repo_updates.json"

# Create metadata directories if they don't exist
BIOINFO_METADATA_DIR.mkdir(exist_ok=True, parents=True)

# GitHub API token
GITHUB_TOKEN = os.environ.get("GITHUB_TOKEN", "")

# Maximum number of concurrent requests
MAX_WORKERS = 5

def sanitize_name(name):
    """Convert tool name to a valid filename."""
    sanitized = re.sub(r'[^\w\-\.]', '_', name)
    return sanitized

def extract_repos_from_readme():
    """Extract repository URLs and names from the README file."""
    try:
        with open(README_PATH, 'r', encoding='utf-8') as f:
            content = f.read()
        
        # Regular expression to match markdown links with tool info
        tool_pattern = r'- \[([^\]]+)\]\((https?://[^\s)]+)\)(.+?)(?:$|\n)'
        
        tools = []
        lines = content.split('\n')
        current_section = None
        current_subsection = None
        
        section_pattern = r'^## ([\w\s&\-]+)$'
        subsection_pattern = r'^### ([\w\s&\-]+)$'
        
        for i, line in enumerate(lines):
            section_match = re.match(section_pattern, line)
            if section_match:
                current_section = section_match.group(1)
                current_subsection = None
                continue
                
            subsection_match = re.match(subsection_pattern, line)
            if subsection_match:
                current_subsection = subsection_match.group(1)
                continue
            
            tool_match = re.match(tool_pattern, line)
            if tool_match and current_section and current_section not in ["Contents", "Contributing", "License"]:
                name = tool_match.group(1)
                url = tool_match.group(2)
                remaining_text = tool_match.group(3).strip()
                
                # Extract description - often after a dash
                description = ""
                if " - " in remaining_text:
                    description = remaining_text.split(" - ", 1)[1].strip()
                else:
                    description = remaining_text
                
                # Extract metadata tags like [Python] or [bioconda]
                metadata_pattern = r'\[([\w\s]+)\]'
                metadata_tags = re.findall(metadata_pattern, remaining_text)
                
                tools.append({
                    "name": name,
                    "url": url.strip(),
                    "description": description,
                    "category": current_section,
                    "subcategory": current_subsection,
                    "metadata_tags": metadata_tags
                })
        
        logger.info(f"Extracted {len(tools)} tools from README.md")
        return tools
    except Exception as e:
        logger.error(f"Error extracting tools from README: {e}")
        return []

def load_repo_data():
    """Load repository data from repo_updates.json if available."""
    try:
        with open(REPO_UPDATES_JSON, 'r') as f:
            repo_data = json.load(f)
        
        # Convert to a dictionary for easier lookup
        repo_dict = {}
        for repo in repo_data:
            repo_dict[repo["url"]] = repo
        
        logger.info(f"Loaded data for {len(repo_dict)} repositories from repo_updates.json")
        return repo_dict
    except (FileNotFoundError, json.JSONDecodeError) as e:
        logger.warning(f"Error loading repo_updates.json: {e}")
        return {}

def load_existing_metadata():
    """Load existing bioinformatics metadata to avoid duplicate work."""
    try:
        if not BIOINFO_METADATA_SUMMARY.exists():
            return {}
            
        with open(BIOINFO_METADATA_SUMMARY, 'r', encoding='utf-8') as f:
            summary = json.load(f)
        
        metadata_dict = {}
        for tool in summary.get("tools", []):
            if "name" in tool and "url" in tool:
                metadata_dict[tool["url"]] = tool
        
        logger.info(f"Loaded existing bioinformatics metadata for {len(metadata_dict)} tools")
        return metadata_dict
    except Exception as e:
        logger.error(f"Error loading existing bioinformatics metadata: {e}")
        return {}

def search_biotools(tool_info):
    """Search for a tool in Bio.tools registry."""
    name = tool_info["name"]
    search_terms = [name]
    
    # Add potential alternative search terms
    if " " in name:
        # Try without spaces
        search_terms.append(name.replace(" ", ""))
    
    # Try with common prefixes/suffixes removed
    for prefix in ["viral", "vir", "bio", "meta"]:
        if name.lower().startswith(prefix):
            search_terms.append(name[len(prefix):].strip())
    
    # For tools with dashes or underscores, try both versions
    if "-" in name:
        search_terms.append(name.replace("-", "_"))
    if "_" in name:
        search_terms.append(name.replace("_", "-"))
    
    # Use description keywords
    description_words = tool_info["description"].lower().split()
    bioinformatics_keywords = [
        "assembly", "annotation", "alignment", "metagenome", "genome", 
        "sequence", "sequencing", "variant", "phylogeny", "taxonomy",
        "classification", "virus", "viral", "virome"
    ]
    
    for word in description_words:
        if word in bioinformatics_keywords and len(word) > 4:
            search_terms.append(word)
    
    # Remove duplicates and sort by length (shorter terms first)
    search_terms = sorted(set(search_terms), key=len)
    
    # Try each search term
    for term in search_terms:
        results = biotools_api.search_tool(term)
        if results and "list" in results and results["list"]:
            # Look for exact or close matches
            for item in results["list"]:
                item_name = item.get("name", "").lower()
                biotoolsID = item.get("biotoolsID", "").lower()
                
                # Check for matches in name, biotoolsID, or description
                name_lower = name.lower()
                if (
                    name_lower == item_name or 
                    name_lower == biotoolsID or
                    name_lower in item_name or
                    item_name in name_lower
                ):
                    logger.info(f"Found match in Bio.tools for {name}: {item.get('biotoolsID')}")
                    # Get detailed information
                    details = biotools_api.get_tool_details(item.get("biotoolsID"))
                    if details:
                        return biotools_api.extract_tool_metadata(details)
    
    logger.info(f"No match found in Bio.tools for {name}")
    return None

def search_bioconda(tool_info):
    """Search for a tool in Bioconda packages."""
    name = tool_info["name"]
    search_terms = [name.lower()]
    
    # Add variations of the name for search
    if " " in name:
        search_terms.append(name.lower().replace(" ", "-"))
        search_terms.append(name.lower().replace(" ", "_"))
    
    # Check if there are tags suggesting conda/bioconda
    has_conda_tag = any(tag.lower() in ["conda", "bioconda"] for tag in tool_info.get("metadata_tags", []))
    
    # Remove duplicates
    search_terms = list(set(search_terms))
    
    # Try each search term
    for term in search_terms:
        package_data = bioconda_api.search_package(term)
        if package_data:
            logger.info(f"Found match in Bioconda for {name}: {package_data.get('name')}")
            
            # Get additional data
            files_data = bioconda_api.get_package_files(package_data.get("name"))
            recipe_data = bioconda_api.get_package_recipe(package_data.get("name"))
            
            return bioconda_api.extract_package_metadata(package_data, files_data, recipe_data)
    
    # If no direct match but has conda tag, try harder with variations
    if has_conda_tag and not any(search_terms):
        # Try more aggressive normalization
        name_normalized = re.sub(r'[^\w]', '', name.lower())
        name_variants = [
            name_normalized,
            "bioconda-" + name_normalized,
            "conda-" + name_normalized
        ]
        
        for variant in name_variants:
            package_data = bioconda_api.search_package(variant)
            if package_data:
                logger.info(f"Found match in Bioconda for {name} (variant: {variant}): {package_data.get('name')}")
                
                # Get additional data
                files_data = bioconda_api.get_package_files(package_data.get("name"))
                recipe_data = bioconda_api.get_package_recipe(package_data.get("name"))
                
                return bioconda_api.extract_package_metadata(package_data, files_data, recipe_data)
    
    logger.info(f"No match found in Bioconda for {name}")
    return None

def extract_github_workflow_info(tool_info, repo_data):
    """Extract GitHub Actions workflow information for compatible workflows."""
    url = tool_info["url"]
    if not url.startswith("https://github.com/"):
        return None
    
    # Extract owner/repo from URL
    match = re.search(r'github\.com/([^/]+/[^/]+)', url)
    if not match:
        return None
    
    repo_path = match.group(1).rstrip('/')
    
    # Check if we have a GitHub token
    if not GITHUB_TOKEN:
        logger.warning(f"No GitHub token available to check workflow info for {repo_path}")
        return None
    
    try:
        # Try to get workflow files
        headers = {"Authorization": f"token {GITHUB_TOKEN}"}
        workflows_url = f"https://api.github.com/repos/{repo_path}/contents/.github/workflows"
        
        response = requests.get(workflows_url, headers=headers)
        if response.status_code != 200:
            return None
            
        workflow_files = response.json()
        
        # Extract workflow information
        workflows = []
        for workflow in workflow_files:
            if workflow["name"].endswith((".yml", ".yaml")):
                # Get the workflow file content
                workflow_content_response = requests.get(workflow["download_url"], headers=headers)
                if workflow_content_response.status_code == 200:
                    content = workflow_content_response.text
                    
                    # Look for common workflow patterns
                    is_nextflow = "nextflow" in content.lower()
                    is_snakemake = "snakemake" in content.lower()
                    is_cwl = "cwl" in content.lower() or "common workflow language" in content.lower()
                    is_wdl = "wdl" in content.lower() or "workflow description language" in content.lower()
                    
                    workflow_info = {
                        "name": workflow["name"],
                        "url": workflow["html_url"],
                        "type": None
                    }
                    
                    if is_nextflow:
                        workflow_info["type"] = "Nextflow"
                    elif is_snakemake:
                        workflow_info["type"] = "Snakemake"
                    elif is_cwl:
                        workflow_info["type"] = "CWL"
                    elif is_wdl:
                        workflow_info["type"] = "WDL"
                    else:
                        workflow_info["type"] = "GitHub Actions"
                    
                    workflows.append(workflow_info)
        
        if workflows:
            return {"workflows": workflows}
        
        return None
    except Exception as e:
        logger.error(f"Error fetching workflow info for {repo_path}: {e}")
        return None

def collect_tool_metadata(tool_info, repo_data, existing_metadata=None):
    """Collect comprehensive bioinformatics metadata for a tool."""
    name = tool_info["name"]
    url = tool_info["url"]
    
    logger.info(f"Collecting bioinformatics metadata for {name} ({url})")
    
    # Check if we already have metadata for this tool
    if existing_metadata and url in existing_metadata:
        existing = existing_metadata[url]
        # Check if the metadata is not too old (over 30 days)
        if "updated_at" in existing:
            try:
                updated_time = datetime.fromisoformat(existing["updated_at"])
                days_since_update = (datetime.now() - updated_time).days
                if days_since_update < 30:
                    logger.info(f"Using existing metadata for {name} (updated {days_since_update} days ago)")
                    return existing
            except (ValueError, TypeError):
                pass
    
    # Initialize metadata with basic info
    metadata = {
        "name": name,
        "url": url,
        "description": tool_info["description"],
        "category": tool_info["category"],
        "subcategory": tool_info["subcategory"],
        "metadata_tags": tool_info.get("metadata_tags", []),
        "updated_at": datetime.now().isoformat()
    }
    
    # Add repository data if available
    if url in repo_data:
        repo_info = repo_data[url]
        metadata["stars"] = repo_info.get("stars")
        metadata["last_updated"] = repo_info.get("last_updated")
    
    # Search Bio.tools for bioinformatics metadata
    biotools_metadata = search_biotools(tool_info)
    if biotools_metadata:
        metadata["biotools"] = biotools_metadata
    
    # Search Bioconda for package metadata
    bioconda_metadata = search_bioconda(tool_info)
    if bioconda_metadata:
        metadata["bioconda"] = bioconda_metadata
    
    # Extract GitHub workflow information
    workflow_info = extract_github_workflow_info(tool_info, repo_data)
    if workflow_info:
        metadata["workflow_info"] = workflow_info
    
    # Consolidate installation methods
    metadata["installation_methods"] = {
        "conda": False,
        "bioconda": False,
        "pip": False,
        "docker": False,
        "singularity": False,
        "source": True  # Default to source installation
    }
    
    # Update from metadata tags
    for tag in tool_info.get("metadata_tags", []):
        tag_lower = tag.lower()
        if tag_lower in ["conda", "bioconda", "pip", "docker", "singularity"]:
            metadata["installation_methods"][tag_lower] = True
    
    # Update from BioTools data
    if biotools_metadata and "installation" in biotools_metadata:
        for method, value in biotools_metadata["installation"].items():
            if method in metadata["installation_methods"] and value:
                metadata["installation_methods"][method] = True
    
    # Update from Bioconda data
    if bioconda_metadata and "installation" in bioconda_metadata:
        for method, value in bioconda_metadata["installation"].items():
            if method in metadata["installation_methods"] and value:
                metadata["installation_methods"][method] = True
    
    # Consolidate input/output formats
    metadata["input_formats"] = []
    metadata["output_formats"] = []
    
    if biotools_metadata:
        if "input_formats" in biotools_metadata:
            metadata["input_formats"].extend(biotools_metadata["input_formats"])
        if "output_formats" in biotools_metadata:
            metadata["output_formats"].extend(biotools_metadata["output_formats"])
    
    # Remove duplicates
    metadata["input_formats"] = list(set(metadata["input_formats"]))
    metadata["output_formats"] = list(set(metadata["output_formats"]))
    
    # Consolidate dependencies
    metadata["dependencies"] = []
    
    if bioconda_metadata and "dependencies" in bioconda_metadata:
        for dep_type, deps in bioconda_metadata["dependencies"].items():
            metadata["dependencies"].extend(deps)
    
    # Remove duplicates
    metadata["dependencies"] = list(set(metadata["dependencies"]))
    
    # Add bioinformatics categories
    metadata["bioinformatics_categories"] = []
    
    if biotools_metadata:
        if "topics" in biotools_metadata:
            metadata["bioinformatics_categories"].extend(biotools_metadata["topics"])
        if "operations" in biotools_metadata:
            metadata["bioinformatics_categories"].extend(biotools_metadata["operations"])
    
    # Add keywords from tags and categories
    inferred_keywords = []
    for tag in tool_info.get("metadata_tags", []):
        if tag.lower() not in ["conda", "bioconda", "pip", "docker", "singularity", "python", "r", "c++", "java"]:
            inferred_keywords.append(tag)
    
    inferred_keywords.append(tool_info["category"])
    if tool_info["subcategory"]:
        inferred_keywords.append(tool_info["subcategory"])
    
    metadata["keywords"] = list(set(inferred_keywords))
    
    # Save the metadata to a file
    save_tool_metadata(metadata)
    
    return metadata

def save_tool_metadata(metadata):
    """Save tool metadata to a JSON file."""
    name = metadata["name"]
    sanitized_name = sanitize_name(name)
    file_path = BIOINFO_METADATA_DIR / f"{sanitized_name}.json"
    
    with open(file_path, 'w', encoding='utf-8') as f:
        json.dump(metadata, f, indent=2)
    
    logger.info(f"Saved bioinformatics metadata for {name} to {file_path}")
    return file_path

def generate_summary_json(metadata_list):
    """Generate a summary JSON file with essential metadata."""
    summary = {
        "tools": [],
        "generated_at": datetime.now().isoformat(),
        "total_count": len(metadata_list),
        "stats": {
            "with_biotools": 0,
            "with_bioconda": 0,
            "installation_methods": {
                "conda": 0,
                "bioconda": 0,
                "pip": 0,
                "docker": 0,
                "singularity": 0,
                "source": 0
            },
            "input_formats": {},
            "output_formats": {},
            "bioinformatics_categories": {}
        }
    }
    
    for metadata in metadata_list:
        # Extract essential information for the summary
        tool_summary = {
            "name": metadata.get("name"),
            "url": metadata.get("url"),
            "description": metadata.get("description"),
            "category": metadata.get("category"),
            "subcategory": metadata.get("subcategory"),
            "stars": metadata.get("stars"),
            "updated_at": metadata.get("updated_at"),
            "installation_methods": metadata.get("installation_methods", {}),
            "input_formats": metadata.get("input_formats", []),
            "output_formats": metadata.get("output_formats", []),
            "dependencies": metadata.get("dependencies", []),
            "bioinformatics_categories": metadata.get("bioinformatics_categories", []),
            "keywords": metadata.get("keywords", [])
        }
        
        # Update statistics
        if "biotools" in metadata:
            summary["stats"]["with_biotools"] += 1
        if "bioconda" in metadata:
            summary["stats"]["with_bioconda"] += 1
        
        # Count installation methods
        for method, value in metadata.get("installation_methods", {}).items():
            if value:
                summary["stats"]["installation_methods"][method] += 1
        
        # Count input/output formats
        for fmt in metadata.get("input_formats", []):
            if fmt not in summary["stats"]["input_formats"]:
                summary["stats"]["input_formats"][fmt] = 0
            summary["stats"]["input_formats"][fmt] += 1
        
        for fmt in metadata.get("output_formats", []):
            if fmt not in summary["stats"]["output_formats"]:
                summary["stats"]["output_formats"][fmt] = 0
            summary["stats"]["output_formats"][fmt] += 1
        
        # Count bioinformatics categories
        for cat in metadata.get("bioinformatics_categories", []):
            if cat not in summary["stats"]["bioinformatics_categories"]:
                summary["stats"]["bioinformatics_categories"][cat] = 0
            summary["stats"]["bioinformatics_categories"][cat] += 1
        
        summary["tools"].append(tool_summary)
    
    # Save the summary JSON
    with open(BIOINFO_METADATA_SUMMARY, 'w', encoding='utf-8') as f:
        json.dump(summary, f, indent=2)
    
    logger.info(f"Generated bioinformatics metadata summary at {BIOINFO_METADATA_SUMMARY}")
    return summary

def batch_process_tools(tools, repo_data, existing_metadata, batch_size=10, batch_delay=5):
    """Process tools in batches to avoid overloading APIs."""
    results = []
    total_tools = len(tools)
    
    for i in range(0, total_tools, batch_size):
        batch = tools[i:i + batch_size]
        logger.info(f"Processing batch {i//batch_size + 1}/{(total_tools + batch_size - 1)//batch_size} ({len(batch)} tools)")
        
        batch_results = []
        with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
            futures = {executor.submit(collect_tool_metadata, tool, repo_data, existing_metadata): tool for tool in batch}
            
            for future in as_completed(futures):
                tool = futures[future]
                try:
                    metadata = future.result()
                    if metadata:
                        batch_results.append(metadata)
                    logger.info(f"Processed {tool['name']} ({tool['url']})")
                except Exception as e:
                    logger.error(f"Error processing {tool['name']} ({tool['url']}): {e}")
        
        results.extend(batch_results)
        
        # If there are more batches to process, wait to avoid rate limiting
        if i + batch_size < total_tools:
            logger.info(f"Waiting {batch_delay} seconds before processing next batch...")
            time.sleep(batch_delay)
    
    return results

def main():
    """Main function to collect bioinformatics metadata."""
    parser = argparse.ArgumentParser(description="Collect bioinformatics metadata for Awesome-Virome tools")
    parser.add_argument("--force-refresh", action="store_true", help="Force refresh all metadata even if recent")
    parser.add_argument("--limit", type=int, help="Limit processing to N tools (for testing)")
    args = parser.parse_args()
    
    # Create metadata directory if it doesn't exist
    BIOINFO_METADATA_DIR.mkdir(exist_ok=True, parents=True)
    
    # Load existing metadata (unless force refresh)
    existing_metadata = {} if args.force_refresh else load_existing_metadata()
    
    # Load repository data (stars, last update, etc.)
    repo_data = load_repo_data()
    
    # Extract tools from README
    tools = extract_repos_from_readme()
    logger.info(f"Found {len(tools)} tools in the README")
    
    # Apply limit if specified
    if args.limit and args.limit > 0:
        tools = tools[:args.limit]
        logger.info(f"Limited processing to {args.limit} tools")
    
    # Process tools in batches
    logger.info("Starting bioinformatics metadata collection...")
    metadata_results = batch_process_tools(tools, repo_data, existing_metadata)
    
    # Generate summary JSON
    summary = generate_summary_json(metadata_results)
    
    successful_count = len(metadata_results)
    logger.info(f"Bioinformatics metadata collection completed. Successfully processed {successful_count}/{len(tools)} tools.")
    
    # Print some statistics
    with_biotools = summary["stats"]["with_biotools"]
    with_bioconda = summary["stats"]["with_bioconda"]
    logger.info(f"Tools with Bio.tools metadata: {with_biotools} ({with_biotools/successful_count*100:.1f}%)")
    logger.info(f"Tools with Bioconda metadata: {with_bioconda} ({with_bioconda/successful_count*100:.1f}%)")
    
    # Installation methods
    for method, count in summary["stats"]["installation_methods"].items():
        if count > 0:
            logger.info(f"Tools with {method} installation: {count} ({count/successful_count*100:.1f}%)")
    
    return 0

if __name__ == "__main__":
    main()