#!/usr/bin/env python3
"""
Script to extract tool information from README.md and update data.json.
Can also incorporate enhanced metadata from the metadata directory.
"""

import re
import json
import logging
import argparse
from pathlib import Path
from datetime import datetime

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler()]
)
logger = logging.getLogger(__name__)

def load_data_json(data_json_path):
    """Load existing data.json or create a new structure if it doesn't exist."""
    try:
        with open(data_json_path, 'r') as f:
            data = json.load(f)
        return data
    except (FileNotFoundError, json.JSONDecodeError):
        # Create a new data structure if the file doesn't exist or is invalid
        return {
            "nodes": [],
            "links": [],
            "categories": [],
            "languages": [],
            "last_updated": datetime.now().isoformat()
        }

def extract_tools_from_readme(readme_path):
    """Extract tool information from README.md."""
    with open(readme_path, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Regular expressions to extract categories, subcategories, and tools
    section_pattern = r'^## ([\w\s&\-]+)$'
    subsection_pattern = r'^### ([\w\s&\-]+)$'
    # Modified pattern to be more permissive with names and capture the full line
    tool_pattern = r'- \[([^\]]+)\]\((https?://[^\s)]+)\)(.+?)(?:$|\n)'
    
    # Additional pattern to extract metadata tags
    metadata_pattern = r'\[([\w\s]+)\]'
    
    tools_data = []
    
    # Process the README line by line
    current_section = None
    current_subsection = None
    
    lines = content.split('\n')
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
            
            # Extract language and package manager from metadata
            metadata_tags = re.findall(metadata_pattern, remaining_text)
            language = None
            package_manager = None
            
            for tag in metadata_tags:
                if tag in ["Python", "R", "Java", "C++", "Perl", "JavaScript", "C", "Ruby", "Go"]:
                    language = tag
                elif tag in ["pip", "conda", "bioconda", "npm", "cargo", "source"]:
                    package_manager = tag
            
            logger.debug(f"Extracted tool: {name} ({url})")
            tools_data.append({
                "name": name,
                "url": url,
                "description": description,
                "category": current_section,
                "subcategory": current_subsection,
                "language": language,
                "package_manager": package_manager
            })
    
    return tools_data

def get_tool_id(name):
    """Generate a tool ID from the tool name."""
    return f"tool-{name.replace(' ', '')}"

def get_category_id(category):
    """Generate a category ID from the category name."""
    return f"category-{category.replace(' ', '-')}"

def get_subcategory_id(subcategory):
    """Generate a subcategory ID from the subcategory name."""
    return f"subcategory-{subcategory.replace(' ', '-')}"

def load_enhanced_metadata():
    """Load enhanced metadata from the metadata directory."""
    metadata_dir = Path(__file__).parent / "metadata"
    summary_path = metadata_dir / "summary.json"
    
    if not summary_path.exists():
        logger.warning(f"Metadata summary file not found at {summary_path}")
        return {}
    
    try:
        with open(summary_path, 'r', encoding='utf-8') as f:
            summary = json.load(f)
        
        # Create a dictionary mapping URLs to metadata
        metadata_by_url = {}
        for repo in summary.get("repositories", []):
            if "url" in repo:
                metadata_by_url[repo["url"]] = repo
                logger.info(f"Added metadata for URL: {repo['url']}")
        
        return metadata_by_url
    except Exception as e:
        logger.error(f"Error loading enhanced metadata: {e}")
        return {}

def sanitize_repo_name(repo_name):
    """Convert repo name to a valid filename."""
    sanitized = re.sub(r'[^\w\-\.]', '_', repo_name)
    return sanitized

def load_individual_metadata(repo_name, repo_url):
    """Load detailed metadata for a specific repository."""
    metadata_dir = Path(__file__).parent / "metadata"
    sanitized_name = sanitize_repo_name(repo_name)
    file_path = metadata_dir / f"{sanitized_name}.json"
    
    if not file_path.exists():
        return None
    
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            return json.load(f)
    except Exception as e:
        logger.warning(f"Error loading metadata for {repo_name}: {e}")
        return None

def update_data_json(data, tools_data, include_metadata=False):
    """Update data.json with extracted tool information and optional enhanced metadata."""
    # Track existing tools for preservation of data
    existing_nodes_by_id = {node["id"]: node for node in data.get("nodes", [])}
    updated_node_ids = set()
    
    # Reset nodes and links arrays
    new_nodes = []
    new_links = []
    
    # Process categories and subcategories
    all_categories = set()
    all_subcategories = {}  # Maps subcategory to parent category
    all_languages = set()
    
    for tool in tools_data:
        all_categories.add(tool["category"])
        if tool["subcategory"]:
            all_subcategories[tool["subcategory"]] = tool["category"]
        if tool["language"]:
            all_languages.add(tool["language"])
    
    # Load enhanced metadata if requested
    enhanced_metadata = {}
    if include_metadata:
        enhanced_metadata = load_enhanced_metadata()
        logger.info(f"Loaded enhanced metadata for {len(enhanced_metadata)} repositories")
    
    # Add category nodes
    for category in sorted(all_categories):
        category_id = get_category_id(category)
        updated_node_ids.add(category_id)
        
        if category_id in existing_nodes_by_id:
            # Preserve existing node but update name
            node = existing_nodes_by_id[category_id].copy()
            node["name"] = category
            new_nodes.append(node)
        else:
            # Create new category node
            new_nodes.append({
                "id": category_id,
                "name": category,
                "type": "category",
                "size": 15,
                "color": "#343a40"
            })
    
    # Add subcategory nodes
    for subcategory, parent_category in all_subcategories.items():
        subcategory_id = get_subcategory_id(subcategory)
        parent_category_id = get_category_id(parent_category)
        updated_node_ids.add(subcategory_id)
        
        if subcategory_id in existing_nodes_by_id:
            # Preserve existing node but update name and parent
            node = existing_nodes_by_id[subcategory_id].copy()
            node["name"] = subcategory
            node["parent"] = parent_category
            new_nodes.append(node)
        else:
            # Create new subcategory node
            new_nodes.append({
                "id": subcategory_id,
                "name": subcategory,
                "type": "subcategory",
                "parent": parent_category,
                "size": 10,
                "color": "#495057"
            })
        
        # Add link between subcategory and parent category
        new_links.append({
            "source": subcategory_id,
            "target": parent_category_id,
            "value": 1
        })
    
    # Count of metadata matches
    metadata_matches = 0
    
    # Add tool nodes
    for tool in tools_data:
        tool_id = get_tool_id(tool["name"])
        updated_node_ids.add(tool_id)
        
        # Create the base node data
        node_data = {
            "id": tool_id,
            "name": tool["name"],
            "type": "tool",
            "url": tool["url"],
            "description": tool["description"],
            "category": tool["category"],
            "subcategory": tool["subcategory"],
            "language": tool["language"],
            "package_manager": tool["package_manager"]
        }
        
        # If the tool already exists, preserve additional properties
        if tool_id in existing_nodes_by_id:
            existing_node = existing_nodes_by_id[tool_id]
            for key, value in existing_node.items():
                if key not in node_data:
                    node_data[key] = value
            
            # Preserve these properties if they exist
            for key in ["stars", "updateTime", "size", "color"]:
                if key in existing_node:
                    node_data[key] = existing_node[key]
        else:
            # Set defaults for new tools
            node_data["stars"] = 0
            node_data["updateTime"] = 0
            node_data["size"] = 10
            node_data["color"] = "#198754"  # Default color
        
        # Add enhanced metadata if available
        logger.info(f"Checking metadata for tool: {tool['name']} ({tool['url']})")
        if include_metadata and tool["url"] in enhanced_metadata:
            logger.info(f"Found enhanced metadata for {tool['name']} ({tool['url']})")
            metadata_matches += 1
            meta = enhanced_metadata[tool["url"]]
            
            # Add basic metadata from summary
            node_data["stars"] = meta.get("stars", node_data.get("stars", 0))
            node_data["forks"] = meta.get("forks", 0)
            node_data["open_issues"] = meta.get("open_issues", 0)
            node_data["license"] = meta.get("license")
            node_data["topics"] = meta.get("topics", [])
            node_data["lastUpdated"] = meta.get("updated_at")
            node_data["createdAt"] = meta.get("created_at")
            node_data["latestRelease"] = meta.get("latest_release")
            node_data["languages"] = meta.get("languages", {})
            
            # If language wasn't set in README, get it from metadata
            if not node_data.get("language") and meta.get("languages"):
                # For GitHub, languages is a dict with language name as key
                if isinstance(meta.get("languages"), dict):
                    languages = list(meta.get("languages").keys())
                    if languages:
                        primary_language = languages[0]
                        node_data["language"] = primary_language
                        all_languages.add(primary_language)
                    node_data["all_languages"] = languages
                # For others, it might be a simple string
                elif isinstance(meta.get("languages"), str):
                    node_data["language"] = meta.get("languages")
                    all_languages.add(meta.get("languages"))
            
            # Load full metadata for additional details
            detailed_meta = load_individual_metadata(tool["name"], tool["url"])
            if detailed_meta:
                # Add any additional fields from detailed metadata that might be useful
                node_data["provider"] = detailed_meta.get("provider")
                node_data["repo_path"] = detailed_meta.get("repo_path")
                node_data["is_archived"] = detailed_meta.get("is_archived", False)
                
                # Adjust node size based on popularity
                stars = node_data.get("stars", 0) or 0
                if stars:
                    if stars > 500:
                        node_data["size"] = 25
                    elif stars > 100:
                        node_data["size"] = 20
                    elif stars > 50:
                        node_data["size"] = 15
                    else:
                        node_data["size"] = 10
        
        new_nodes.append(node_data)
        
        # Add link between tool and category/subcategory
        if tool["subcategory"]:
            subcategory_id = get_subcategory_id(tool["subcategory"])
            new_links.append({
                "source": tool_id,
                "target": subcategory_id,
                "value": 1
            })
        else:
            category_id = get_category_id(tool["category"])
            new_links.append({
                "source": tool_id,
                "target": category_id,
                "value": 1
            })
    
    logger.info(f"Added enhanced metadata to {metadata_matches} tools")
    
    # Update the data object
    data["nodes"] = new_nodes
    data["links"] = new_links
    data["categories"] = sorted(all_categories)
    data["languages"] = sorted(all_languages)
    data["last_updated"] = datetime.now().isoformat()
    
    return data

def save_data_json(data, data_json_path):
    """Save data to data.json."""
    with open(data_json_path, 'w') as f:
        json.dump(data, f, indent=2)
    logger.info(f"Updated {data_json_path}")

def load_bioinformatics_metadata():
    """Load bioinformatics metadata from the metadata/bioinformatics directory."""
    metadata_dir = Path(__file__).parent / "metadata" / "bioinformatics"
    summary_path = metadata_dir / "summary.json"
    
    if not summary_path.exists():
        logger.warning(f"Bioinformatics metadata summary file not found at {summary_path}")
        return {}
    
    try:
        with open(summary_path, 'r', encoding='utf-8') as f:
            summary = json.load(f)
        
        # Create a dictionary mapping URLs to metadata
        metadata_by_url = {}
        for tool in summary.get("tools", []):
            if "url" in tool:
                metadata_by_url[tool["url"]] = tool
                logger.info(f"Added bioinformatics metadata for URL: {tool['url']}")
        
        return metadata_by_url
    except Exception as e:
        logger.error(f"Error loading bioinformatics metadata: {e}")
        return {}

def update_data_json_with_bioinformatics(data, bioinformatics_metadata):
    """Update data.json with bioinformatics metadata."""
    if not bioinformatics_metadata:
        return data
    
    # Get all nodes that are tools
    tool_nodes = [node for node in data.get("nodes", []) if node.get("type") == "tool"]
    
    # Count of metadata matches
    metadata_matches = 0
    
    # Update tool nodes with bioinformatics metadata
    for node in tool_nodes:
        if node.get("url") and node.get("url") in bioinformatics_metadata:
            metadata = bioinformatics_metadata[node.get("url")]
            metadata_matches += 1
            
            # Add bioinformatics fields
            node["input_formats"] = metadata.get("input_formats", [])
            node["output_formats"] = metadata.get("output_formats", [])
            node["bioinformatics_categories"] = metadata.get("bioinformatics_categories", [])
            node["dependencies"] = metadata.get("dependencies", [])
            
            # Update installation methods
            if "installation_methods" not in node:
                node["installation_methods"] = {}
            
            for method, value in metadata.get("installation_methods", {}).items():
                node["installation_methods"][method] = value
    
    logger.info(f"Added bioinformatics metadata to {metadata_matches} tools")
    return data

def main():
    """Main function to extract tools from README.md and update data.json."""
    parser = argparse.ArgumentParser(description="Extract tool information from README.md and update data.json")
    parser.add_argument("--include-metadata", action="store_true", help="Include enhanced metadata from metadata directory")
    parser.add_argument("--include-bioinformatics-metadata", action="store_true", help="Include bioinformatics metadata")
    args = parser.parse_args()
    
    repo_root = Path(__file__).parent
    readme_path = repo_root / "README.md"
    data_json_path = repo_root / "data.json"
    
    logger.info(f"Loading existing data from {data_json_path}")
    data = load_data_json(data_json_path)
    
    logger.info(f"Extracting tools from {readme_path}")
    tools_data = extract_tools_from_readme(readme_path)
    logger.info(f"Found {len(tools_data)} tools in README.md")
    
    logger.info(f"Updating {data_json_path}" + (" with enhanced metadata" if args.include_metadata else ""))
    updated_data = update_data_json(data, tools_data, include_metadata=args.include_metadata)
    
    # Load and add bioinformatics metadata if requested
    if args.include_bioinformatics_metadata:
        logger.info("Loading bioinformatics metadata...")
        bioinformatics_metadata = load_bioinformatics_metadata()
        if bioinformatics_metadata:
            logger.info(f"Loaded bioinformatics metadata for {len(bioinformatics_metadata)} tools")
            updated_data = update_data_json_with_bioinformatics(updated_data, bioinformatics_metadata)
    
    logger.info(f"Saving to {data_json_path}")
    save_data_json(updated_data, data_json_path)
    
    logger.info("Update completed successfully")

if __name__ == "__main__":
    main()