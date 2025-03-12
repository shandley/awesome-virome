#!/usr/bin/env python3
"""
Script to verify that all tools in data.json appear in README.md and check for duplicates.
"""

import os
import re
import sys
import json
import argparse
import logging
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

def load_data_json(data_json_path):
    """Load and parse the data.json file."""
    try:
        with open(data_json_path, 'r') as f:
            data = json.load(f)
        return data
    except (FileNotFoundError, json.JSONDecodeError) as e:
        logger.error(f"Error loading data.json: {e}")
        sys.exit(1)

def load_readme(readme_path):
    """Load the README.md file."""
    try:
        with open(readme_path, 'r') as f:
            readme_content = f.read()
        return readme_content
    except FileNotFoundError as e:
        logger.error(f"Error loading README.md: {e}")
        sys.exit(1)

def extract_tools_from_data_json(data):
    """Extract tool information from data.json."""
    tools = []
    for node in data.get('nodes', []):
        if node.get('type') == 'tool':
            tools.append({
                'name': node.get('name'),
                'url': node.get('url'),
                'category': node.get('category')
            })
    return tools

def extract_links_from_readme(readme_content):
    """Extract markdown links from README.md."""
    link_pattern = r'\[([^\]]+)\]\(([^)]+)\)'
    links = re.findall(link_pattern, readme_content)
    
    # Convert to a more usable format
    return [{'name': name, 'url': url} for name, url in links]

def check_tools_in_readme(data_tools, readme_links):
    """Check if all tools from data.json appear in README.md."""
    readme_urls = [link['url'] for link in readme_links]
    missing_tools = []
    
    for tool in data_tools:
        if tool['url'] and tool['url'] not in readme_urls:
            missing_tools.append(tool)
    
    return missing_tools

def check_duplicates_in_readme(readme_content):
    """Check for duplicate tool entries in README.md."""
    # Regular expression to match markdown links
    link_pattern = r'\[([^\]]+)\]\(([^)]+)\)'
    links = re.findall(link_pattern, readme_content)
    
    # Check for duplicate URLs
    url_counts = {}
    duplicate_urls = {}
    
    for name, url in links:
        # Skip non-repository URLs
        if not (url.startswith('http') and 
                ('github.com' in url or 'gitlab.com' in url or 'bitbucket.org' in url)):
            continue
        
        url_counts[url] = url_counts.get(url, 0) + 1
        if url_counts[url] > 1:
            if url not in duplicate_urls:
                duplicate_urls[url] = []
            duplicate_urls[url].append(name)
    
    # Check for duplicate names
    name_counts = {}
    duplicate_names = {}
    
    for name, url in links:
        # Skip non-tool entries (e.g., section links, images)
        if not (url.startswith('http') and 
                ('github.com' in url or 'gitlab.com' in url or 'bitbucket.org' in url)):
            continue
        
        name_counts[name] = name_counts.get(name, 0) + 1
        if name_counts[name] > 1:
            if name not in duplicate_names:
                duplicate_names[name] = []
            duplicate_names[name].append(url)
    
    return duplicate_urls, duplicate_names

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Verify README.md content against data.json.')
    parser.add_argument('--data-json', default='data.json',
                        help='Path to data.json file (default: data.json)')
    parser.add_argument('--readme', default='README.md',
                        help='Path to README.md file (default: README.md)')
    parser.add_argument('--check-duplicates', action='store_true',
                        help='Check for duplicate entries in README.md')
    
    return parser.parse_args()

def main():
    """Main function to verify README content against data.json."""
    args = parse_arguments()
    
    repo_root = Path(__file__).parent.parent
    data_json_path = repo_root / args.data_json
    readme_path = repo_root / args.readme
    
    # Load data
    data = load_data_json(data_json_path)
    readme_content = load_readme(readme_path)
    
    # If only checking for duplicates
    if args.check_duplicates:
        duplicate_urls, duplicate_names = check_duplicates_in_readme(readme_content)
        
        if duplicate_urls or duplicate_names:
            if duplicate_urls:
                logger.error("Found duplicate URLs in README.md:")
                for url, names in duplicate_urls.items():
                    logger.error(f"  URL: {url} appears with names: {', '.join(names)}")
            
            if duplicate_names:
                logger.error("Found duplicate tool names in README.md:")
                for name, urls in duplicate_names.items():
                    logger.error(f"  Tool name: {name} appears with URLs: {', '.join(urls)}")
            
            sys.exit(1)
        else:
            logger.info("No duplicate entries found in README.md")
            sys.exit(0)
    
    # Extract tools and links
    data_tools = extract_tools_from_data_json(data)
    readme_links = extract_links_from_readme(readme_content)
    
    # Check if all tools are in README
    missing_tools = check_tools_in_readme(data_tools, readme_links)
    
    if missing_tools:
        logger.error("The following tools from data.json are missing in README.md:")
        for tool in missing_tools:
            logger.error(f"  - {tool['name']} ({tool['url']}) in category: {tool['category']}")
        sys.exit(1)
    else:
        logger.info("All tools from data.json are present in README.md")
        sys.exit(0)

if __name__ == "__main__":
    main()