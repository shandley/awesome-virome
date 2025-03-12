#!/usr/bin/env python3
"""
Script to verify repository information integrity in data.json and repo_updates.json.
"""

import os
import re
import sys
import json
import logging
import requests
from datetime import datetime
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# GitHub API base URL
GITHUB_API = "https://api.github.com/repos/"
GITLAB_API = "https://gitlab.com/api/v4/projects/"
BITBUCKET_API = "https://api.bitbucket.org/2.0/repositories/"

# GitHub API token (set as environment variable GITHUB_TOKEN)
GITHUB_TOKEN = os.environ.get("GITHUB_TOKEN", "")

# Maximum number of concurrent requests
MAX_WORKERS = 3

def load_data_json(data_json_path):
    """Load and parse the data.json file."""
    try:
        with open(data_json_path, 'r') as f:
            data = json.load(f)
        return data
    except (FileNotFoundError, json.JSONDecodeError) as e:
        logger.error(f"Error loading data.json: {e}")
        sys.exit(1)

def load_repo_updates(repo_updates_path):
    """Load the repo_updates.json file if it exists."""
    try:
        with open(repo_updates_path, 'r') as f:
            data = json.load(f)
        return data
    except FileNotFoundError:
        logger.warning(f"repo_updates.json not found at {repo_updates_path}. Creating a new file.")
        return []
    except json.JSONDecodeError as e:
        logger.error(f"Error parsing repo_updates.json: {e}")
        sys.exit(1)

def extract_tools_from_data_json(data):
    """Extract tool information from data.json."""
    tools = []
    for node in data.get('nodes', []):
        if node.get('type') == 'tool' and node.get('url'):
            tools.append({
                'name': node.get('name'),
                'url': node.get('url'),
                'stars': node.get('stars'),
                'category': node.get('category')
            })
    return tools

def get_github_repo_info(repo_url):
    """Get information for a GitHub repository."""
    try:
        # Extract owner/repo from URL
        match = re.search(r'github\.com/([^/]+/[^/]+)', repo_url)
        if not match:
            return None, "invalid_url", None
        
        repo_path = match.group(1).rstrip('/')
        
        # Remove any extra path components
        if "/" in repo_path.split("/", 2)[1]:
            repo_path = "/".join(repo_path.split("/", 2)[:2])
        
        # Create API URL
        api_url = f"{GITHUB_API}{repo_path}"
        
        headers = {"Accept": "application/vnd.github.v3+json"}
        if GITHUB_TOKEN:
            headers["Authorization"] = f"token {GITHUB_TOKEN}"
        
        response = requests.get(api_url, headers=headers, timeout=10)
        
        if response.status_code == 404:
            return None, "not_found", None
            
        response.raise_for_status()
        
        data = response.json()
        updated_at = data.get('updated_at') or data.get('pushed_at')
        stars = data.get('stargazers_count', 0)
        
        if updated_at:
            return datetime.strptime(updated_at, '%Y-%m-%dT%H:%M:%SZ'), None, stars
        return None, None, stars
    except Exception as e:
        logger.error(f"Error checking GitHub repo {repo_url}: {e}")
        return None, "error", None

def get_gitlab_repo_info(repo_url):
    """Get information for a GitLab repository."""
    try:
        # Extract project ID or path from URL
        match = re.search(r'gitlab\.com/([^/]+/[^/]+)', repo_url)
        if not match:
            return None, "invalid_url", None
        
        project_path = match.group(1).rstrip('/')
        
        # URL encode the project path
        encoded_path = requests.utils.quote(project_path, safe='')
        
        # Create API URL for project lookup
        api_url = f"{GITLAB_API}?search={encoded_path}"
        
        response = requests.get(api_url, timeout=10)
        
        if response.status_code == 404:
            return None, "not_found", None
            
        response.raise_for_status()
        
        projects = response.json()
        for project in projects:
            if project['path_with_namespace'].lower() == project_path.lower():
                project_id = project['id']
                project_url = f"{GITLAB_API}{project_id}"
                
                # Get detailed project info
                proj_response = requests.get(project_url, timeout=10)
                proj_response.raise_for_status()
                
                data = proj_response.json()
                updated_at = data.get('last_activity_at')
                stars = data.get('star_count', 0)
                
                if updated_at:
                    return datetime.strptime(updated_at, '%Y-%m-%dT%H:%M:%S.%fZ'), None, stars
                return None, None, stars
        return None, "not_found", None
    except Exception as e:
        logger.error(f"Error checking GitLab repo {repo_url}: {e}")
        return None, "error", None

def get_bitbucket_repo_info(repo_url):
    """Get information for a Bitbucket repository."""
    try:
        # Extract workspace/repo from URL
        match = re.search(r'bitbucket\.org/([^/]+/[^/]+)', repo_url)
        if not match:
            return None, "invalid_url", None
        
        repo_path = match.group(1).rstrip('/')
        
        # Create API URL
        api_url = f"{BITBUCKET_API}{repo_path}"
        
        response = requests.get(api_url, timeout=10)
        
        if response.status_code == 404:
            return None, "not_found", None
            
        response.raise_for_status()
        
        data = response.json()
        updated_on = data.get('updated_on')
        # Bitbucket doesn't have stars but has watchers
        watchers = data.get('watchers_count', 0)
        
        if updated_on:
            return datetime.strptime(updated_on, '%Y-%m-%dT%H:%M:%S.%f%z'), None, watchers
        return None, None, watchers
    except Exception as e:
        logger.error(f"Error checking Bitbucket repo {repo_url}: {e}")
        return None, "error", None

def verify_repo_info(tool):
    """Verify repository information for a single tool."""
    repo_url = tool['url']
    name = tool['name']
    
    if 'github.com' in repo_url:
        updated_at, status, stars = get_github_repo_info(repo_url)
    elif 'gitlab.com' in repo_url:
        updated_at, status, stars = get_gitlab_repo_info(repo_url)
    elif 'bitbucket.org' in repo_url:
        updated_at, status, stars = get_bitbucket_repo_info(repo_url)
    else:
        return {'name': name, 'url': repo_url, 'status': 'unsupported'}
    
    result = {
        'name': name,
        'url': repo_url,
        'last_updated': updated_at.isoformat() if updated_at else None,
        'status': status,
        'stars': stars,
        'data_json_stars': tool.get('stars')
    }
    
    # Check for discrepancies in star counts
    if stars is not None and tool.get('stars') is not None and abs(stars - tool.get('stars')) > 10:
        result['star_discrepancy'] = True
    
    return result

def batch_verify_repos(tools):
    """Verify repository information for multiple tools in parallel."""
    results = []
    
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        future_to_tool = {executor.submit(verify_repo_info, tool): tool for tool in tools}
        
        for future in as_completed(future_to_tool):
            tool = future_to_tool[future]
            try:
                result = future.result()
                results.append(result)
                logger.info(f"Verified {tool['name']} ({tool['url']})")
            except Exception as e:
                logger.error(f"Error verifying {tool['name']} ({tool['url']}): {e}")
                results.append({
                    'name': tool['name'],
                    'url': tool['url'],
                    'status': 'error',
                    'error': str(e)
                })
    
    return results

def analyze_verification_results(results):
    """Analyze verification results and report any issues."""
    issues = []
    
    # Check for unavailable repositories
    unavailable = [r for r in results if r.get('status') == 'not_found']
    if unavailable:
        issues.append("Unavailable repositories:")
        for repo in unavailable:
            issues.append(f"  - {repo['name']} ({repo['url']})")
    
    # Check for star count discrepancies
    star_discrepancies = [r for r in results if r.get('star_discrepancy')]
    if star_discrepancies:
        issues.append("Star count discrepancies:")
        for repo in star_discrepancies:
            issues.append(f"  - {repo['name']} ({repo['url']}): data.json: {repo['data_json_stars']}, actual: {repo['stars']}")
    
    # Check for errors
    errors = [r for r in results if r.get('status') == 'error']
    if errors:
        issues.append("Error checking repositories:")
        for repo in errors:
            issues.append(f"  - {repo['name']} ({repo['url']}): {repo.get('error', 'Unknown error')}")
    
    # Check for unsupported repository hosts
    unsupported = [r for r in results if r.get('status') == 'unsupported']
    if unsupported:
        issues.append("Unsupported repository hosts:")
        for repo in unsupported:
            issues.append(f"  - {repo['name']} ({repo['url']})")
    
    # Check for invalid URLs
    invalid_urls = [r for r in results if r.get('status') == 'invalid_url']
    if invalid_urls:
        issues.append("Invalid repository URLs:")
        for repo in invalid_urls:
            issues.append(f"  - {repo['name']} ({repo['url']})")
    
    return issues

def save_verification_results(results, output_path):
    """Save verification results to a JSON file."""
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2)
    logger.info(f"Verification results saved to {output_path}")

def main():
    """Main function to verify repository metadata."""
    # Check if GITHUB_TOKEN is set
    if not GITHUB_TOKEN:
        logger.warning("GITHUB_TOKEN environment variable not set. GitHub API rate limiting will be more restrictive.")
    
    repo_root = Path(__file__).parent.parent
    data_json_path = repo_root / 'data.json'
    repo_updates_path = repo_root / 'repo_updates.json'
    verification_results_path = repo_root / 'verification_results.json'
    
    # Load data
    data = load_data_json(data_json_path)
    tools = extract_tools_from_data_json(data)
    
    logger.info(f"Found {len(tools)} tools with URLs in data.json")
    
    # Verify repository information
    results = batch_verify_repos(tools)
    
    # Analyze results
    issues = analyze_verification_results(results)
    
    # Save results
    save_verification_results(results, verification_results_path)
    
    # Report issues
    if issues:
        for issue in issues:
            logger.error(issue)
        sys.exit(1)
    else:
        logger.info("No repository metadata issues found")
        sys.exit(0)

if __name__ == "__main__":
    main()