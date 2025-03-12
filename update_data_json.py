#!/usr/bin/env python3
"""
Script to extract tool information from README.md and update data.json.
"""

import re
import json
import logging
from pathlib import Path

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
            "languages": []
        }

def extract_tools_from_readme(readme_path):
    """Extract tool information from README.md."""
    with open(readme_path, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Regular expressions to extract categories, subcategories, and tools
    section_pattern = r'^## ([\w\s&\-]+)$'
    subsection_pattern = r'^### ([\w\s&\-]+)$'
    tool_pattern = r'- \[([\w\s\-\.]+)\]\((https?://[^\s)]+)\) - (.*?)(?:$|\n)'
    
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
            description = tool_match.group(3).strip()
            
            # Extract language and package manager from metadata
            metadata_tags = re.findall(metadata_pattern, line)
            language = None
            package_manager = None
            
            for tag in metadata_tags:
                if tag in ["Python", "R", "Java", "C++", "Perl", "JavaScript", "C", "Ruby", "Go"]:
                    language = tag
                elif tag in ["pip", "conda", "bioconda", "npm", "cargo", "source"]:
                    package_manager = tag
            
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

def update_data_json(data, tools_data):
    """Update data.json with extracted tool information."""
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
    
    # Update the data object
    data["nodes"] = new_nodes
    data["links"] = new_links
    data["categories"] = sorted(all_categories)
    data["languages"] = sorted(all_languages)
    
    return data

def save_data_json(data, data_json_path):
    """Save data to data.json."""
    with open(data_json_path, 'w') as f:
        json.dump(data, f, indent=2)
    logger.info(f"Updated {data_json_path}")

def main():
    """Main function to extract tools from README.md and update data.json."""
    repo_root = Path(__file__).parent
    readme_path = repo_root / "README.md"
    data_json_path = repo_root / "data.json"
    
    logger.info(f"Loading existing data from {data_json_path}")
    data = load_data_json(data_json_path)
    
    logger.info(f"Extracting tools from {readme_path}")
    tools_data = extract_tools_from_readme(readme_path)
    logger.info(f"Found {len(tools_data)} tools in README.md")
    
    logger.info(f"Updating {data_json_path}")
    updated_data = update_data_json(data, tools_data)
    
    logger.info(f"Saving to {data_json_path}")
    save_data_json(updated_data, data_json_path)
    
    logger.info("Update completed successfully")

if __name__ == "__main__":
    main()