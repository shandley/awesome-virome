#!/usr/bin/env python3
"""
Script to fetch and store enhanced metadata for repositories in the Awesome-Virome list.
This script collects detailed metadata such as:
- License information
- Programming languages
- Repository topics/tags
- Latest release information
- Dependencies (where available)
- Issue count
- Repository creation date

The metadata is stored in structured JSON files in the metadata/ directory.

This script also supports smart cache invalidation, which ensures that when repository
information is updated, all related cached API responses are automatically invalidated.

Additionally, this module provides functions for extracting tool-specific metadata for
integration with the incremental metadata update process.
"""

import os
import re
import sys
import json
import time
import logging
import requests
from pathlib import Path
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed
from github import Github, RateLimitExceededException

# Import the cache manager
try:
    from apis.citations_api import cache_manager
    HAS_CACHE_MANAGER = True
except ImportError:
    HAS_CACHE_MANAGER = False
    print("Warning: Smart cache system not found. Cache invalidation will be disabled.")

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("metadata_enhancement.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# GitHub API token (set as environment variable GITHUB_TOKEN)
GITHUB_TOKEN = os.environ.get("GITHUB_TOKEN", "")
if not GITHUB_TOKEN:
    logger.warning("GITHUB_TOKEN environment variable not set. Rate limiting will be more restrictive.")

# Base paths
REPO_ROOT = Path(__file__).parent.parent
METADATA_DIR = REPO_ROOT / "metadata"
REPO_UPDATES_JSON = REPO_ROOT / "repo_updates.json"
README_PATH = REPO_ROOT / "README.md"

# Create metadata directory if it doesn't exist
METADATA_DIR.mkdir(exist_ok=True)

# Maximum number of concurrent requests
MAX_WORKERS = 3

# Rate limiting settings
RATE_LIMIT_DELAY = 2.0

class RateLimiter:
    """Simple rate limiter for API requests."""
    def __init__(self, requests_per_minute=30):
        self.delay = 60.0 / requests_per_minute
        self.last_request_time = 0
        
    def wait(self):
        """Wait if needed to respect rate limits."""
        current_time = time.time()
        time_since_last_request = current_time - self.last_request_time
        
        if time_since_last_request < self.delay:
            sleep_time = self.delay - time_since_last_request
            time.sleep(sleep_time)
            
        self.last_request_time = time.time()

# Create rate limiters for different APIs
github_limiter = RateLimiter(requests_per_minute=30)
gitlab_limiter = RateLimiter(requests_per_minute=30)
bitbucket_limiter = RateLimiter(requests_per_minute=30)

def load_repo_data():
    """Load repository data from repo_updates.json."""
    try:
        with open(REPO_UPDATES_JSON, 'r') as f:
            return json.load(f)
    except (FileNotFoundError, json.JSONDecodeError) as e:
        logger.error(f"Error loading repo_updates.json: {e}")
        logger.info("Falling back to extracting repos from README")
        return extract_repos_from_readme()

def extract_repos_from_readme():
    """Extract repository URLs from the README file."""
    try:
        with open(README_PATH, 'r', encoding='utf-8') as f:
            content = f.read()
        
        # Regular expression to match markdown links
        link_pattern = r'\[([^\]]+)\]\(([^)]+)\)'
        links = re.findall(link_pattern, content)
        
        repos = []
        for name, url in links:
            # Filter for GitHub, GitLab, and Bitbucket repositories
            if ('github.com/' in url or 'gitlab.com/' in url or 'bitbucket.org/' in url) and not url.endswith('.md'):
                repos.append({
                    "name": name,
                    "url": url.strip(),
                    "last_updated": None,
                    "status": None,
                    "stars": None
                })
        
        return repos
    except Exception as e:
        logger.error(f"Error extracting repos from README: {e}")
        return []

def sanitize_repo_name(repo_name):
    """Convert repo name to a valid filename."""
    sanitized = re.sub(r'[^\w\-\.]', '_', repo_name)
    return sanitized

def get_github_enhanced_metadata(repo_url, repo_name):
    """Get enhanced metadata for a GitHub repository."""
    try:
        github_limiter.wait()
        
        # Extract owner/repo from URL
        match = re.search(r'github\.com/([^/]+/[^/]+)', repo_url)
        if not match:
            logger.warning(f"Invalid GitHub URL format: {repo_url}")
            return None
        
        repo_path = match.group(1).rstrip('/')
        # Handle release URLs
        repo_path = repo_path.split('/releases/')[0]
        
        # Use GitHub REST API directly for better error handling without a token
        api_url = f"https://api.github.com/repos/{repo_path}"
        headers = {"Accept": "application/vnd.github.v3+json"}
        if GITHUB_TOKEN:
            headers["Authorization"] = f"Bearer {GITHUB_TOKEN}"
        
        logger.info(f"Fetching GitHub repo metadata: {repo_path}")
        response = requests.get(api_url, headers=headers)
        
        if response.status_code != 200:
            if response.status_code == 404:
                logger.warning(f"Repository not found: {repo_url}")
                return None
            elif response.status_code == 403:
                # Rate limit exceeded
                logger.warning(f"Rate limit exceeded for {repo_url}. Consider setting GITHUB_TOKEN.")
                return None
            else:
                logger.warning(f"Failed to fetch GitHub repo: {repo_url}, status code: {response.status_code}")
                return None
                
        repo_data = response.json()
        
        # Basic repository info
        metadata = {
            "name": repo_name,
            "url": repo_url,
            "provider": "github",
            "repo_path": repo_path,
            "description": repo_data.get('description'),
            "created_at": repo_data.get('created_at'),
            "updated_at": repo_data.get('updated_at'),
            "pushed_at": repo_data.get('pushed_at'),
            "stars": repo_data.get('stargazers_count'),
            "watchers": repo_data.get('watchers_count'),
            "forks": repo_data.get('forks_count'),
            "open_issues": repo_data.get('open_issues_count'),
            "default_branch": repo_data.get('default_branch'),
            "license": repo_data.get('license', {}).get('name') if repo_data.get('license') else None,
            "license_url": repo_data.get('license', {}).get('url') if repo_data.get('license') else None,
            "homepage": repo_data.get('homepage'),
            "is_archived": repo_data.get('archived', False),
            "is_disabled": repo_data.get('disabled', False),
            "fetch_time": datetime.now().isoformat(),
        }
        
        # Get languages
        try:
            github_limiter.wait()
            langs_url = f"{api_url}/languages"
            langs_response = requests.get(langs_url, headers=headers)
            if langs_response.status_code == 200:
                metadata["languages"] = langs_response.json()
            else:
                metadata["languages"] = {}
        except Exception as e:
            logger.warning(f"Error fetching languages for {repo_path}: {e}")
            metadata["languages"] = {}
        
        # Get topics/tags
        try:
            github_limiter.wait()
            topics_url = f"{api_url}/topics"
            topics_headers = headers.copy()
            topics_headers["Accept"] = "application/vnd.github.mercy-preview+json"
            topics_response = requests.get(topics_url, headers=topics_headers)
            if topics_response.status_code == 200:
                json_response = topics_response.json()
                if isinstance(json_response, dict) and 'names' in json_response:
                    metadata["topics"] = json_response.get('names', [])
                else:
                    metadata["topics"] = []
                    logger.warning(f"Unexpected format in topics response for {repo_path}")
            else:
                metadata["topics"] = []
                logger.warning(f"Failed to fetch topics for {repo_path}: Status {topics_response.status_code}")
        except Exception as e:
            logger.warning(f"Error fetching topics for {repo_path}: {e}")
            metadata["topics"] = []
        
        # Get latest release if available
        try:
            github_limiter.wait()
            release_url = f"{api_url}/releases/latest"
            release_response = requests.get(release_url, headers=headers)
            if release_response.status_code == 200:
                latest_release = release_response.json()
                metadata["latest_release"] = {
                    "name": latest_release.get('name'),
                    "tag": latest_release.get('tag_name'),
                    "published_at": latest_release.get('published_at'),
                    "url": latest_release.get('html_url')
                }
            else:
                metadata["latest_release"] = None
        except Exception as e:
            logger.warning(f"Error fetching latest release for {repo_path}: {e}")
            metadata["latest_release"] = None
        
        # Try to get dependencies if available through the dependency graph API
        # This requires additional permissions and may not always work
        if GITHUB_TOKEN:  # Only try this with a token
            try:
                github_limiter.wait()
                deps_url = f"{api_url}/dependency-graph/sbom"
                deps_response = requests.get(deps_url, headers=headers)
                if deps_response.status_code == 200:
                    metadata["dependencies"] = deps_response.json()
                else:
                    metadata["dependencies"] = None
            except Exception as e:
                logger.warning(f"Error fetching dependencies for {repo_path}: {e}")
                metadata["dependencies"] = None
        else:
            metadata["dependencies"] = None
        
        return metadata
    
    except RateLimitExceededException:
        logger.warning(f"GitHub API rate limit exceeded. Waiting and will retry later.")
        time.sleep(60)  # Wait a minute before allowing other operations
        return None
    except Exception as e:
        logger.error(f"Error fetching enhanced metadata for GitHub repo {repo_url}: {str(e)}")
        return None

def get_gitlab_enhanced_metadata(repo_url, repo_name):
    """Get enhanced metadata for a GitLab repository."""
    try:
        gitlab_limiter.wait()
        
        # Extract project path from URL
        match = re.search(r'gitlab\.com/([^/]+/[^/]+)', repo_url)
        if not match:
            logger.warning(f"Invalid GitLab URL format: {repo_url}")
            return None
        
        project_path = match.group(1).rstrip('/')
        
        # URL encode the project path
        encoded_path = requests.utils.quote(project_path, safe='')
        
        # Create API URL for project lookup
        api_url = f"https://gitlab.com/api/v4/projects?search={encoded_path}"
        
        logger.info(f"Fetching GitLab repo metadata: {project_path}")
        response = requests.get(api_url)
        
        if response.status_code != 200:
            logger.warning(f"Failed to fetch GitLab project: {repo_url}")
            return None
            
        projects = response.json()
        project_id = None
        
        for project in projects:
            if project['path_with_namespace'].lower() == project_path.lower():
                project_id = project['id']
                break
        
        if not project_id:
            logger.warning(f"GitLab project not found: {repo_url}")
            return None
        
        # Get detailed project info
        gitlab_limiter.wait()
        project_url = f"https://gitlab.com/api/v4/projects/{project_id}"
        proj_response = requests.get(project_url)
        
        if proj_response.status_code != 200:
            return None
            
        data = proj_response.json()
        
        # Get languages
        gitlab_limiter.wait()
        langs_url = f"https://gitlab.com/api/v4/projects/{project_id}/languages"
        langs_response = requests.get(langs_url)
        languages = langs_response.json() if langs_response.status_code == 200 else {}
        
        # Basic repository info
        metadata = {
            "name": repo_name,
            "url": repo_url,
            "provider": "gitlab",
            "repo_path": project_path,
            "description": data.get('description'),
            "created_at": data.get('created_at'),
            "updated_at": data.get('last_activity_at'),
            "stars": data.get('star_count', 0),
            "forks": data.get('forks_count', 0),
            "open_issues": data.get('open_issues_count', 0),
            "default_branch": data.get('default_branch'),
            "license": data.get('license', {}).get('name'),
            "license_url": data.get('license', {}).get('url'),
            "homepage": data.get('web_url'),
            "is_archived": data.get('archived', False),
            "languages": languages,
            "topics": data.get('topics', []),
            "fetch_time": datetime.now().isoformat(),
        }
        
        # Get latest release if available
        try:
            gitlab_limiter.wait()
            release_url = f"https://gitlab.com/api/v4/projects/{project_id}/releases"
            release_response = requests.get(release_url)
            
            if release_response.status_code == 200 and release_response.json():
                latest_release = release_response.json()[0]  # First release is latest
                metadata["latest_release"] = {
                    "name": latest_release.get('name'),
                    "tag": latest_release.get('tag_name'),
                    "published_at": latest_release.get('released_at'),
                    "url": latest_release.get('_links', {}).get('self')
                }
            else:
                metadata["latest_release"] = None
        except Exception:
            metadata["latest_release"] = None
        
        return metadata
    except Exception as e:
        logger.error(f"Error fetching enhanced metadata for GitLab repo {repo_url}: {e}")
        return None

def get_bitbucket_enhanced_metadata(repo_url, repo_name):
    """Get enhanced metadata for a Bitbucket repository."""
    try:
        bitbucket_limiter.wait()
        
        # Extract workspace/repo from URL
        match = re.search(r'bitbucket\.org/([^/]+/[^/]+)', repo_url)
        if not match:
            logger.warning(f"Invalid Bitbucket URL format: {repo_url}")
            return None
        
        repo_path = match.group(1).rstrip('/')
        
        # Create API URL
        api_url = f"https://api.bitbucket.org/2.0/repositories/{repo_path}"
        
        logger.info(f"Fetching Bitbucket repo metadata: {repo_path}")
        response = requests.get(api_url)
        
        if response.status_code != 200:
            logger.warning(f"Failed to fetch Bitbucket repository: {repo_url}")
            return None
            
        data = response.json()
        
        # Basic repository info
        metadata = {
            "name": repo_name,
            "url": repo_url,
            "provider": "bitbucket",
            "repo_path": repo_path,
            "description": data.get('description'),
            "created_at": data.get('created_on'),
            "updated_at": data.get('updated_on'),
            "default_branch": data.get('mainbranch', {}).get('name'),
            "language": data.get('language'),
            "is_private": data.get('is_private', False),
            "homepage": data.get('website'),
            "fetch_time": datetime.now().isoformat(),
        }
        
        # Get watchers and forks
        if 'links' in data:
            if 'watchers' in data['links']:
                bitbucket_limiter.wait()
                watchers_url = data['links']['watchers']['href']
                watchers_response = requests.get(watchers_url)
                if watchers_response.status_code == 200:
                    metadata["watchers"] = watchers_response.json().get('size', 0)
            
            if 'forks' in data['links']:
                bitbucket_limiter.wait()
                forks_url = data['links']['forks']['href']
                forks_response = requests.get(forks_url)
                if forks_response.status_code == 200:
                    metadata["forks"] = forks_response.json().get('size', 0)
        
        return metadata
    except Exception as e:
        logger.error(f"Error fetching enhanced metadata for Bitbucket repo {repo_url}: {e}")
        return None

def get_enhanced_metadata(repo):
    """Get enhanced metadata for a repository based on its URL."""
    repo_name = repo["name"]
    repo_url = repo["url"]
    
    logger.info(f"Getting enhanced metadata for {repo_name} ({repo_url})")
    
    # Skip repositories marked as unavailable
    if repo.get("status") == "not_found":
        logger.info(f"Skipping unavailable repository: {repo_name}")
        return None
    
    # Handle different repository hosts
    if 'github.com' in repo_url:
        return get_github_enhanced_metadata(repo_url, repo_name)
    elif 'gitlab.com' in repo_url:
        return get_gitlab_enhanced_metadata(repo_url, repo_name)
    elif 'bitbucket.org' in repo_url:
        return get_bitbucket_enhanced_metadata(repo_url, repo_name)
    else:
        logger.warning(f"Unsupported repository host for {repo_url}")
        return None

def save_metadata(metadata):
    """Save metadata to a JSON file in the metadata directory."""
    if not metadata:
        return
    
    repo_name = metadata["name"]
    repo_url = metadata.get("url", "")
    sanitized_name = sanitize_repo_name(repo_name)
    file_path = METADATA_DIR / f"{sanitized_name}.json"
    
    # Check if the repository was previously saved and has been updated
    if file_path.exists():
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                old_metadata = json.load(f)
                
            # If repository has been updated, invalidate related caches
            if HAS_CACHE_MANAGER and old_metadata.get("updated_at") != metadata.get("updated_at"):
                logger.info(f"Repository {repo_name} has been updated. Invalidating related caches...")
                invalidated = cache_manager.invalidate_repo_caches(repo_url)
                logger.info(f"Invalidated {invalidated} cache entries for {repo_name}")
        except Exception as e:
            logger.warning(f"Error comparing metadata for cache invalidation: {e}")
    
    with open(file_path, 'w', encoding='utf-8') as f:
        json.dump(metadata, f, indent=2)
    
    logger.info(f"Saved metadata for {repo_name} to {file_path}")
    return file_path

def generate_summary_json():
    """Generate a summary JSON file with essential metadata from all repositories."""
    metadata_files = list(METADATA_DIR.glob("*.json"))
    summary = {
        "repositories": [],
        "generated_at": datetime.now().isoformat(),
        "total_count": len(metadata_files)
    }
    
    for file_path in metadata_files:
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                metadata = json.load(f)
            
            # Extract essential information for the summary
            repo_summary = {
                "name": metadata.get("name"),
                "url": metadata.get("url"),
                "provider": metadata.get("provider"),
                "description": metadata.get("description"),
                "stars": metadata.get("stars"),
                "forks": metadata.get("forks"),
                "open_issues": metadata.get("open_issues"),
                "license": metadata.get("license"),
                "languages": metadata.get("languages"),
                "topics": metadata.get("topics"),
                "updated_at": metadata.get("updated_at"),
                "created_at": metadata.get("created_at"),
                "latest_release": metadata.get("latest_release", {}).get("tag") if metadata.get("latest_release") else None
            }
            
            summary["repositories"].append(repo_summary)
        except Exception as e:
            logger.error(f"Error processing metadata file {file_path}: {e}")
    
    # Sort repositories by stars (if available)
    summary["repositories"].sort(key=lambda x: x.get("stars", 0) if x.get("stars") is not None else 0, reverse=True)
    
    # Save the summary JSON
    summary_path = METADATA_DIR / "summary.json"
    with open(summary_path, 'w', encoding='utf-8') as f:
        json.dump(summary, f, indent=2)
    
    logger.info(f"Generated metadata summary at {summary_path}")
    return summary_path

def batch_process_repos(repos, batch_size=5, batch_delay=10):
    """Process repositories in batches to avoid rate limiting."""
    results = []
    total_repos = len(repos)
    
    for i in range(0, total_repos, batch_size):
        batch = repos[i:i + batch_size]
        logger.info(f"Processing batch {i//batch_size + 1}/{(total_repos + batch_size - 1)//batch_size} ({len(batch)} repos)")
        
        batch_results = []
        with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
            futures = {executor.submit(get_enhanced_metadata, repo): repo for repo in batch}
            
            for future in as_completed(futures):
                repo = futures[future]
                try:
                    metadata = future.result()
                    if metadata:
                        save_metadata(metadata)
                        batch_results.append(metadata)
                    logger.info(f"Processed {repo['name']} ({repo['url']})")
                except Exception as e:
                    logger.error(f"Error processing {repo['name']} ({repo['url']}): {e}")
        
        results.extend(batch_results)
        
        # If there are more batches to process, wait to avoid rate limiting
        if i + batch_size < total_repos:
            logger.info(f"Waiting {batch_delay} seconds before processing next batch...")
            time.sleep(batch_delay)
    
    return results

def main():
    """Main function to enhance metadata for repositories."""
    # Create metadata directory if it doesn't exist
    METADATA_DIR.mkdir(exist_ok=True)
    
    # Load repository data
    logger.info("Loading repository data...")
    repos = load_repo_data()
    
    if not repos:
        logger.error("No repositories found. Exiting.")
        return
    
    logger.info(f"Found {len(repos)} repositories for metadata enhancement")
    
    # Process repositories in batches - reduce batch size without a token
    batch_size = 3 if not GITHUB_TOKEN else 5
    batch_delay = 20 if not GITHUB_TOKEN else 15
    
    # Process repositories in batches
    logger.info("Starting enhanced metadata collection...")
    metadata_results = batch_process_repos(repos, batch_size=batch_size, batch_delay=batch_delay)
    
    # Generate summary JSON - even with partial or empty results
    logger.info("Generating metadata summary...")
    summary_path = generate_summary_json()
    
    successful_count = len([m for m in metadata_results if m is not None])
    logger.info(f"Enhanced metadata collection completed. Successfully processed {successful_count}/{len(repos)} repositories.")
    logger.info(f"Metadata summary available at {summary_path}")

# Functions for incremental metadata update integration
def get_repo_metadata(repo_url, github_api=None):
    """
    Get repository metadata for a given URL.
    This is a wrapper function used by the incremental update process.
    
    Args:
        repo_url: URL of the repository
        github_api: Optional GitHub API client
        
    Returns:
        Dictionary containing repository metadata
    """
    # Simple utility to extract repository name from URL
    def extract_repo_name(url):
        if 'github.com' in url:
            match = re.search(r'github\.com/([^/]+)/([^/]+)', url)
            if match:
                return match.group(2).split('.')[0]
        elif 'gitlab.com' in url:
            match = re.search(r'gitlab\.com/([^/]+)/([^/]+)', url)
            if match:
                return match.group(2).split('.')[0]
        elif 'bitbucket.org' in url:
            match = re.search(r'bitbucket\.org/([^/]+)/([^/]+)', url)
            if match:
                return match.group(2).split('.')[0]
        return None

    repo_name = extract_repo_name(repo_url)
    
    # Get metadata based on repository host
    if 'github.com' in repo_url:
        return get_github_enhanced_metadata(repo_url, repo_name or "")
    elif 'gitlab.com' in repo_url:
        return get_gitlab_enhanced_metadata(repo_url, repo_name or "")
    elif 'bitbucket.org' in repo_url:
        return get_bitbucket_enhanced_metadata(repo_url, repo_name or "")
    else:
        return {"name": repo_name, "url": repo_url, "provider": "unknown"}

def extract_tool_metadata(tool_id, repo_url, tool_name, github_api=None, semantic_scholar_api=None, crossref_api=None):
    """
    Extract tool-specific metadata for a given tool.
    This function integrates with academic_impact.py to collect citation data.
    
    Args:
        tool_id: Identifier for the tool
        repo_url: URL of the repository
        tool_name: Name of the tool
        github_api: GitHub API client (optional)
        semantic_scholar_api: Semantic Scholar API client (optional)
        crossref_api: CrossRef API client (optional)
        
    Returns:
        Dictionary containing tool-specific metadata
    """
    # Import the academic impact collector lazily to avoid circular imports
    try:
        from academic_impact import AcademicImpactCollector
        academic_collector = AcademicImpactCollector(
            github_token=os.environ.get("GITHUB_TOKEN"),
            semantic_scholar_key=os.environ.get("SEMANTIC_SCHOLAR_KEY"),
            contact_email=os.environ.get("CONTACT_EMAIL")
        )
        
        # Create a minimal tool object for the academic impact collector
        tool_obj = {
            "id": tool_id,
            "name": tool_name,
            "url": repo_url,
            "description": ""
        }
        
        # Collect academic impact data
        academic_impact = academic_collector.process_tool(tool_obj)
        
        # Extract relevant fields
        return {
            "academic_impact": {
                "doi": academic_impact.get("doi"),
                "citation_info": academic_impact.get("citation_info"),
                "citation_metrics": academic_impact.get("citation_metrics")
            }
        }
    except (ImportError, Exception) as e:
        logger.warning(f"Error collecting academic impact data: {e}")
        return {
            "academic_impact": {
                "doi": None,
                "citation_info": {},
                "citation_metrics": {}
            }
        }

if __name__ == "__main__":
    main()