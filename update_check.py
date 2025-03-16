#!/usr/bin/env python3
"""
Script to automatically check when GitHub/GitLab/Bitbucket repositories in the README were last updated.
Also identifies unavailable repositories (404 errors).
"""

import re
import os
import sys
import argparse
import time
import json
import logging
import requests
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Tuple, Optional, Dict, List, Any

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("update_check.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Define API base URLs
GITHUB_API = "https://api.github.com/repos/"
GITLAB_API = "https://gitlab.com/api/v4/projects/"
BITBUCKET_API = "https://api.bitbucket.org/2.0/repositories/"

# Maximum number of concurrent requests
MAX_WORKERS = 3

# Rate limiting delay (seconds) - increased to avoid GitHub rate limiting
RATE_LIMIT_DELAY = 2.0

# GitHub API token (set as environment variable GITHUB_TOKEN)
GITHUB_TOKEN = os.environ.get("GITHUB_TOKEN", "")

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
github_limiter = RateLimiter(requests_per_minute=30)  # GitHub allows 60 requests/hour without token
gitlab_limiter = RateLimiter(requests_per_minute=30)
bitbucket_limiter = RateLimiter(requests_per_minute=30)

def extract_repos_from_readme(readme_path):
    """Extract repository URLs from the README file."""
    with open(readme_path, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Regular expression to match markdown links
    link_pattern = r'\[([^\]]+)\]\(([^)]+)\)'
    links = re.findall(link_pattern, content)
    
    repos = []
    for name, url in links:
        # Filter for GitHub, GitLab, and Bitbucket repositories
        if ('github.com/' in url or 'gitlab.com/' in url or 'bitbucket.org/' in url) and not url.endswith('.md'):
            repos.append((name, url.strip()))
    
    return repos

def get_github_repo_info(repo_url):
    """Get the last updated time, star count, and version info for a GitHub repository."""
    try:
        # Use rate limiter
        github_limiter.wait()
        
        # Extract owner/repo from URL
        match = re.search(r'github\.com/([^/]+/[^/]+)', repo_url)
        if not match:
            return None, None, None, None
        
        repo_path = match.group(1).rstrip('/')
        # Handle branch or subdirectory in URL
        if '/tree/' in repo_path:
            repo_path = repo_path.split('/tree/')[0]
        
        # Create API URL
        api_url = f"{GITHUB_API}{repo_path}"
        
        headers = {"Accept": "application/vnd.github.v3+json"}
        if GITHUB_TOKEN:
            headers["Authorization"] = f"token {GITHUB_TOKEN}"
        
        logger.info(f"Fetching GitHub repo: {repo_path}")
        response = requests.get(api_url, headers=headers)
        
        # Check for rate limiting
        if response.status_code == 403 and 'X-RateLimit-Remaining' in response.headers:
            remaining = int(response.headers['X-RateLimit-Remaining'])
            reset_time = int(response.headers['X-RateLimit-Reset'])
            reset_datetime = datetime.fromtimestamp(reset_time)
            now = datetime.now()
            wait_time = (reset_datetime - now).total_seconds()
            
            if remaining == 0 and wait_time > 0:
                logger.warning(f"GitHub API rate limit exceeded. Reset at {reset_datetime}. Waiting {wait_time} seconds.")
                time.sleep(min(wait_time + 1, 3600))  # Wait until reset, but not more than an hour
                return get_github_repo_info(repo_url)  # Retry after waiting
        
        if response.status_code == 404:
            logger.warning(f"Repository not found: {repo_url}")
            return None, "not_found", None, None
        
        if response.status_code == 429:
            logger.warning(f"Rate limit hit for {repo_url}, waiting and retrying")
            time.sleep(10)  # Wait 10 seconds before retrying
            return get_github_repo_info(repo_url)  # Retry
            
        response.raise_for_status()
        
        data = response.json()
        updated_at = data.get('updated_at') or data.get('pushed_at')
        stars = data.get('stargazers_count', 0)
        
        # Get version information from releases
        version_info = None
        if GITHUB_TOKEN:  # Only fetch releases if we have a token to avoid rate limits
            try:
                github_limiter.wait()
                releases_url = f"{GITHUB_API}{repo_path}/releases"
                releases_response = requests.get(
                    releases_url, 
                    headers=headers, 
                    params={"per_page": 1}  # Only get the latest release
                )
                
                if releases_response.status_code == 200:
                    releases = releases_response.json()
                    if releases:
                        latest_release = releases[0]
                        tag_name = latest_release.get('tag_name', '')
                        # Clean up tag name to get version
                        version = tag_name
                        if version.startswith('v'):
                            version = version  # Keep the v prefix
                        elif version.startswith('release-'):
                            version = 'v' + version[8:]
                        
                        # Get release date
                        release_date = latest_release.get('published_at')
                        if release_date:
                            release_year = datetime.strptime(release_date, '%Y-%m-%dT%H:%M:%SZ').year
                            version_info = (version, release_year)
            except Exception as e:
                logger.warning(f"Error fetching releases for {repo_path}: {e}")
        
        if updated_at:
            return datetime.strptime(updated_at, '%Y-%m-%dT%H:%M:%SZ'), None, stars, version_info
        return None, None, stars, version_info
    except requests.exceptions.RequestException as e:
        if hasattr(e, 'response') and e.response is not None:
            if e.response.status_code == 404:
                logger.warning(f"Repository not found: {repo_url}")
                return None, "not_found", None, None
            elif e.response.status_code == 429:
                logger.warning(f"Rate limit hit for {repo_url}, waiting 60 seconds")
                time.sleep(60)  # Wait 60 seconds before retrying
                return get_github_repo_info(repo_url)  # Retry
        logger.error(f"Error fetching GitHub repo {repo_url}: {e}")
        return None, "error", None, None
    except Exception as e:
        logger.error(f"Unexpected error for {repo_url}: {e}")
        return None, "error", None, None

def get_gitlab_repo_info(repo_url):
    """Get the last updated time, star count, and version info for a GitLab repository."""
    try:
        # Use rate limiter
        gitlab_limiter.wait()
        
        # Extract project ID or path from URL
        match = re.search(r'gitlab\.com/([^/]+/[^/]+)', repo_url)
        if not match:
            return None, None, None, None
        
        project_path = match.group(1).rstrip('/')
        # Handle branch or subdirectory in URL
        if '/tree/' in project_path:
            project_path = project_path.split('/tree/')[0]
            
        # URL encode the project path
        encoded_path = requests.utils.quote(project_path, safe='')
        
        # Create API URL for project lookup
        api_url = f"{GITLAB_API}?search={encoded_path}"
        
        logger.info(f"Fetching GitLab repo: {project_path}")
        response = requests.get(api_url)
        
        if response.status_code == 404:
            logger.warning(f"Repository not found: {repo_url}")
            return None, "not_found", None, None
        
        if response.status_code == 429:
            logger.warning(f"Rate limit hit for {repo_url}, waiting 60 seconds")
            time.sleep(60)  # Wait 60 seconds before retrying
            return get_gitlab_repo_info(repo_url)  # Retry
            
        response.raise_for_status()
        
        projects = response.json()
        for project in projects:
            if project['path_with_namespace'].lower() == project_path.lower():
                project_id = project['id']
                project_url = f"{GITLAB_API}{project_id}"
                
                # Get detailed project info
                gitlab_limiter.wait()
                proj_response = requests.get(project_url)
                proj_response.raise_for_status()
                
                data = proj_response.json()
                updated_at = data.get('last_activity_at')
                stars = data.get('star_count', 0)
                
                # Get version information from tags
                version_info = None
                try:
                    gitlab_limiter.wait()
                    tags_url = f"{GITLAB_API}{project_id}/repository/tags"
                    tags_response = requests.get(
                        tags_url,
                        params={"per_page": 1, "order_by": "updated", "sort": "desc"}
                    )
                    
                    if tags_response.status_code == 200:
                        tags = tags_response.json()
                        if tags:
                            latest_tag = tags[0]
                            tag_name = latest_tag.get('name', '')
                            # Clean up tag name to get version
                            version = tag_name
                            if version.startswith('v'):
                                version = version  # Keep the v prefix
                            elif version.startswith('release-'):
                                version = 'v' + version[8:]
                            
                            # Get release date (use commit date as GitLab doesn't have release dates)
                            if 'commit' in latest_tag and 'committed_date' in latest_tag['commit']:
                                commit_date = latest_tag['commit']['committed_date']
                                if commit_date:
                                    commit_year = datetime.strptime(commit_date, '%Y-%m-%dT%H:%M:%S.%fZ').year
                                    version_info = (version, commit_year)
                except Exception as e:
                    logger.warning(f"Error fetching tags for {project_path}: {e}")
                
                if updated_at:
                    return datetime.strptime(updated_at, '%Y-%m-%dT%H:%M:%S.%fZ'), None, stars, version_info
                return None, None, stars, version_info
        return None, None, None, None
    except requests.exceptions.RequestException as e:
        if hasattr(e, 'response') and e.response is not None:
            if e.response.status_code == 404:
                logger.warning(f"Repository not found: {repo_url}")
                return None, "not_found", None, None
            elif e.response.status_code == 429:
                logger.warning(f"Rate limit hit for {repo_url}, waiting 60 seconds")
                time.sleep(60)  # Wait 60 seconds before retrying
                return get_gitlab_repo_info(repo_url)  # Retry
        logger.error(f"Error fetching GitLab repo {repo_url}: {e}")
        return None, "error", None, None
    except Exception as e:
        logger.error(f"Unexpected error for {repo_url}: {e}")
        return None, "error", None, None

def get_bitbucket_repo_info(repo_url):
    """Get the last updated time, watchers count, and version info for a Bitbucket repository."""
    try:
        # Use rate limiter
        bitbucket_limiter.wait()
        
        # Extract workspace/repo from URL
        match = re.search(r'bitbucket\.org/([^/]+/[^/]+)', repo_url)
        if not match:
            return None, None, None, None
        
        repo_path = match.group(1).rstrip('/')
        # Handle branch or subdirectory in URL
        if '/src/' in repo_path:
            repo_path = repo_path.split('/src/')[0]
        
        # Create API URL
        api_url = f"{BITBUCKET_API}{repo_path}"
        
        logger.info(f"Fetching Bitbucket repo: {repo_path}")
        response = requests.get(api_url)
        
        if response.status_code == 404:
            logger.warning(f"Repository not found: {repo_url}")
            return None, "not_found", None, None
        
        if response.status_code == 429:
            logger.warning(f"Rate limit hit for {repo_url}, waiting 60 seconds")
            time.sleep(60)  # Wait 60 seconds before retrying
            return get_bitbucket_repo_info(repo_url)  # Retry
            
        response.raise_for_status()
        
        data = response.json()
        updated_on = data.get('updated_on')
        # Bitbucket doesn't have stars but has watchers
        watchers = data.get('watchers_count', 0)
        
        # Get version information from tags
        version_info = None
        try:
            bitbucket_limiter.wait()
            tags_url = f"{BITBUCKET_API}{repo_path}/refs/tags?sort=-target.date&pagelen=1"
            tags_response = requests.get(tags_url)
            
            if tags_response.status_code == 200:
                tags_data = tags_response.json()
                if 'values' in tags_data and tags_data['values']:
                    latest_tag = tags_data['values'][0]
                    tag_name = latest_tag.get('name', '')
                    # Clean up tag name to get version
                    version = tag_name
                    if version.startswith('v'):
                        version = version  # Keep the v prefix
                    elif version.startswith('release-'):
                        version = 'v' + version[8:]
                    
                    # Get tag date
                    if 'target' in latest_tag and 'date' in latest_tag['target']:
                        tag_date = latest_tag['target']['date']
                        if tag_date:
                            # Bitbucket date format might vary
                            try:
                                tag_year = datetime.strptime(tag_date, '%Y-%m-%dT%H:%M:%S.%f%z').year
                                version_info = (version, tag_year)
                            except ValueError:
                                try:
                                    tag_year = datetime.strptime(tag_date, '%Y-%m-%dT%H:%M:%S%z').year
                                    version_info = (version, tag_year)
                                except ValueError:
                                    logger.warning(f"Couldn't parse date format for {repo_path} tag: {tag_date}")
        except Exception as e:
            logger.warning(f"Error fetching tags for {repo_path}: {e}")
        
        if updated_on:
            return datetime.strptime(updated_on, '%Y-%m-%dT%H:%M:%S.%f%z'), None, watchers, version_info
        return None, None, watchers, version_info
    except requests.exceptions.RequestException as e:
        if hasattr(e, 'response') and e.response is not None:
            if e.response.status_code == 404:
                logger.warning(f"Repository not found: {repo_url}")
                return None, "not_found", None, None
            elif e.response.status_code == 429:
                logger.warning(f"Rate limit hit for {repo_url}, waiting 60 seconds")
                time.sleep(60)  # Wait 60 seconds before retrying
                return get_bitbucket_repo_info(repo_url)  # Retry
        logger.error(f"Error fetching Bitbucket repo {repo_url}: {e}")
        return None, "error", None, None
    except Exception as e:
        logger.error(f"Unexpected error for {repo_url}: {e}")
        return None, "error", None, None

def get_repo_last_updated(repo_name, repo_url):
    """Get the last updated time, stars, and version for a repository based on its URL."""
    # Handle different repository hosts
    if 'github.com' in repo_url:
        last_updated, status, stars, version_info = get_github_repo_info(repo_url)
        return repo_name, repo_url, last_updated, status, stars, version_info
    elif 'gitlab.com' in repo_url:
        last_updated, status, stars, version_info = get_gitlab_repo_info(repo_url)
        return repo_name, repo_url, last_updated, status, stars, version_info
    elif 'bitbucket.org' in repo_url:
        last_updated, status, stars, version_info = get_bitbucket_repo_info(repo_url)
        return repo_name, repo_url, last_updated, status, stars, version_info
    return repo_name, repo_url, None, None, None, None

def batch_process_repos(repos, batch_size=5, batch_delay=10):
    """Process repositories in batches to avoid rate limiting."""
    results = []
    total_repos = len(repos)
    
    for i in range(0, total_repos, batch_size):
        batch = repos[i:i + batch_size]
        logger.info(f"Processing batch {i//batch_size + 1}/{(total_repos + batch_size - 1)//batch_size} ({len(batch)} repos)")
        
        batch_results = []
        with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
            futures = {executor.submit(get_repo_last_updated, name, url): (name, url) for name, url in batch}
            
            for future in as_completed(futures):
                name, url = futures[future]
                try:
                    result = future.result()
                    batch_results.append(result)
                    logger.info(f"Processed {name} ({url})")
                except Exception as e:
                    logger.error(f"Error processing {name} ({url}): {e}")
                    batch_results.append((name, url, None, "error", None, None))
        
        results.extend(batch_results)
        
        # If there are more batches to process, wait to avoid rate limiting
        if i + batch_size < total_repos:
            logger.info(f"Waiting {batch_delay} seconds before processing next batch...")
            time.sleep(batch_delay)
    
    return results

def extract_repo_categories(readme_path):
    """Extract repositories and their categories from the README file."""
    with open(readme_path, 'r', encoding='utf-8') as f:
        content = f.read()

    # Define the main categories we want to track
    main_categories = {
        "Virus and Phage Identification": [
            "## Virus and Phage Identification",
            "### Metagenome Analysis",
            "### Integrated Viruses",
            "### RNA Virus Identification"
        ],
        "Host Prediction": ["## Host Prediction"],
        "Genome Analysis": [
            "## Genome Analysis",
            "### Genome Annotation",
            "### Genome Assembly",
            "### Genome Completeness",
            "### Genome Comparison",
            "### Gene Finding"
        ],
        "Taxonomy": ["## Taxonomy"]
    }
    
    # Map each section header to its main category
    section_to_category = {}
    for main_cat, sections in main_categories.items():
        for section in sections:
            section_to_category[section] = main_cat
    
    # Initialize repo categories
    repo_categories = {}
    current_section = None
    current_category = None
    
    # Regular expression to match markdown links
    link_pattern = r'\[([^\]]+)\]\(([^)]+)\)'
    
    # Process the README line by line
    for line in content.split('\n'):
        # Check if this line is a section header
        if line.startswith('#'):
            for section, category in section_to_category.items():
                if line.strip() == section:
                    current_section = section
                    current_category = category
                    break
        
        # If we're in a tracked section, look for repository links
        if current_category:
            matches = re.findall(link_pattern, line)
            for name, url in matches:
                # Consider only GitHub, GitLab, and Bitbucket repositories
                if ('github.com/' in url or 'gitlab.com/' in url or 'bitbucket.org/' in url) and not url.endswith('.md'):
                    if name not in repo_categories:
                        repo_categories[name] = current_category
    
    return repo_categories

def update_readme_with_dates_status_and_stars(readme_path, repo_data):
    """Update the README with the last updated information, availability status, and create a popular packages section."""
    with open(readme_path, 'r', encoding='utf-8') as f:
        content = f.read()
    
    updated_content = content
    unavailable_repos = []
    github_repos_with_stars = []
    
    # Extract repository categories
    repo_categories = extract_repo_categories(readme_path)
    
    # Sort repositories by name for consistent processing
    sorted_repos = sorted(repo_data, key=lambda x: x[0])
    
    for data in sorted_repos:
        # Handle both old and new format of data tuple
        if len(data) == 6:
            repo_name, repo_url, last_updated, status, stars, version_info = data
        else:
            repo_name, repo_url, last_updated, status, stars = data
            version_info = None
            
        # Collect GitHub repos with stars for popular packages section
        if 'github.com' in repo_url and stars is not None and stars > 0:
            category = repo_categories.get(repo_name, "Other")
            github_repos_with_stars.append((repo_name, repo_url, stars, category))
            
        # Check if repository is unavailable
        if status == "not_found":
            unavailable_repos.append((repo_name, repo_url))
            
            # Check if the repo is already marked as unavailable
            repo_pattern = re.escape(f"[{repo_name}]({repo_url})") + r".*?(\[unavailable\])?"
            unavailable_match = re.search(repo_pattern, updated_content)
            
            if unavailable_match and unavailable_match.group(1):
                # Already marked as unavailable, no need to change
                pass
            else:
                # Mark as unavailable
                # Replace all existing tags with the new tag
                # First remove all existing tags
                temp_content = re.sub(
                    re.escape(f"[{repo_name}]({repo_url})") + r"(?:\s+\[.*?\])*",
                    f"[{repo_name}]({repo_url})",
                    updated_content
                )
                # Then add the new tag
                updated_content = temp_content.replace(
                    f"[{repo_name}]({repo_url})",
                    f"[{repo_name}]({repo_url}) [unavailable]"
                )
        else:
            # Get existing version tag from README if it exists
            version_pattern = re.escape(f"[{repo_name}]({repo_url})") + r".*?(\[v[\d\.]+,\s*\d{4}\])"
            existing_version_match = re.search(version_pattern, updated_content)
            
            # Keep track of tags to add
            tags = []
            
            # Add update date if available
            if last_updated:
                year = last_updated.year
                month = last_updated.month
                tags.append(f"[Updated: {month:02d}/{year}]")
            
            # Add version info if available from API and not already in README
            if version_info and not existing_version_match:
                version, year = version_info
                tags.append(f"[{version}, {year}]")
            elif existing_version_match:
                # Keep existing version tag
                tags.append(existing_version_match.group(1))
            
            if tags:
                # First remove all existing tags
                tag_removal_pattern = re.escape(f"[{repo_name}]({repo_url})") + r"(?:\s+\[.*?\])*"
                entry_without_tags = re.search(tag_removal_pattern, updated_content).group(0)
                base_entry = f"[{repo_name}]({repo_url})"
                
                # Build new entry with all tags
                new_entry = base_entry + " " + " ".join(tags)
                
                # Replace the entire entry
                updated_content = updated_content.replace(entry_without_tags, new_entry)
    
    # Create overall popular packages section
    if github_repos_with_stars:
        top_repos = sorted(github_repos_with_stars, key=lambda x: x[2], reverse=True)[:20]  # Top 20 repos by stars
        
        popular_packages_section = "## Popular Packages\n\nRanked by GitHub stars:\n\n"
        for i, (name, url, stars, _) in enumerate(top_repos, 1):
            popular_packages_section += f"{i}. [{name}]({url}) - ⭐ {stars} stars\n"
        popular_packages_section += "\n"
        
        # Check if the section already exists
        if "## Popular Packages" in updated_content:
            # Replace existing section
            pattern = r"## Popular Packages.*?(?=\n## |\Z)"
            updated_content = re.sub(pattern, popular_packages_section, updated_content, flags=re.DOTALL)
        else:
            # Add after the introduction
            updated_content = updated_content.replace(
                "## Getting Started", 
                f"{popular_packages_section}## Getting Started"
            )
    
    # Create "Top Packages by Category" section
    if github_repos_with_stars:
        # Group repositories by category
        repos_by_category = {}
        for name, url, stars, category in github_repos_with_stars:
            if category not in repos_by_category:
                repos_by_category[category] = []
            repos_by_category[category].append((name, url, stars))
        
        # Sort repositories within each category by stars
        for category in repos_by_category:
            repos_by_category[category] = sorted(repos_by_category[category], key=lambda x: x[2], reverse=True)
        
        # Create the section content
        top_packages_section = "## Top Packages by Category\n\nHere are the most starred packages in key categories:\n\n"
        
        # Include only specific categories of interest
        categories_of_interest = [
            "Virus and Phage Identification",
            "Host Prediction",
            "Genome Analysis",
            "Taxonomy"
        ]
        
        for category in categories_of_interest:
            if category in repos_by_category and repos_by_category[category]:
                top_packages_section += f"### {category}\n"
                # Take top 3 from each category
                for i, (name, url, stars) in enumerate(repos_by_category[category][:3], 1):
                    top_packages_section += f"{i}. [{name}]({url}) - ⭐ {stars} stars\n"
                top_packages_section += "\n"
        
        # Check if the section already exists
        if "## Top Packages by Category" in updated_content:
            # Replace existing section
            pattern = r"## Top Packages by Category.*?(?=\n## |\Z)"
            updated_content = re.sub(pattern, top_packages_section, updated_content, flags=re.DOTALL)
        else:
            # Add after Popular Packages section
            if "## Popular Packages" in updated_content:
                updated_content = updated_content.replace(
                    "## Getting Started", 
                    f"{top_packages_section}## Getting Started"
                )
            else:
                # Add before Getting Started if Popular Packages doesn't exist
                updated_content = updated_content.replace(
                    "## Getting Started", 
                    f"{top_packages_section}## Getting Started"
                )
    
    # Write the updated content back to the README
    with open(readme_path, 'w', encoding='utf-8') as f:
        f.write(updated_content)
    
    # Also create a file with the list of unavailable repositories
    with open(Path(__file__).parent / "unavailable_repos.md", 'w') as f:
        f.write("# Unavailable Repositories\n\n")
        f.write("The following repositories were not found (404 errors):\n\n")
        for name, url in unavailable_repos:
            f.write(f"- [{name}]({url})\n")
    
    # Create a file with the list of starred repositories
    with open(Path(__file__).parent / "starred_repos.md", 'w') as f:
        f.write("# GitHub Repositories by Stars\n\n")
        
        if github_repos_with_stars:
            sorted_by_stars = sorted(github_repos_with_stars, key=lambda x: x[2], reverse=True)
            f.write("| Repository | Category | Stars |\n")
            f.write("|------------|----------|-------|\n")
            for name, url, stars, category in sorted_by_stars:
                f.write(f"| [{name}]({url}) | {category} | {stars} |\n")
        else:
            f.write("No GitHub repositories with stars found.\n")
    
    return updated_content, unavailable_repos, github_repos_with_stars

def main(readme):
    """Main function to check and update repository information."""
    # Check if GITHUB_TOKEN is set
    if not GITHUB_TOKEN:
        logger.warning("GITHUB_TOKEN environment variable not set. Rate limiting will be more restrictive.")
    
    readme_path = os.path.join(Path(__file__).parent, readme)
    
    print(f"Checking and updating repository information in {readme_path}")
    if not os.path.exists(readme_path):
        logger.error(f"README file not found at {readme_path}")
        sys.exit(1)
    
    logger.info(f"Extracting repositories from {readme_path}...")
    repos = extract_repos_from_readme(readme_path)
    logger.info(f"Found {len(repos)} repository URLs in the README")
    
    # Extract repository categories
    repo_categories = extract_repo_categories(readme_path)
    logger.info(f"Categorized {len(repo_categories)} repositories into their respective sections")
    
    # Process repositories in batches
    results = batch_process_repos(repos, batch_size=10, batch_delay=15)
    
    # Count statistics
    updated_repos = [r for r in results if r[2] is not None]
    unavailable_repos = [r for r in results if r[3] == "not_found"]
    error_repos = [r for r in results if r[3] == "error"]
    github_repos_with_stars = [r for r in results if 'github.com' in r[1] and r[4] is not None and r[4] > 0]
    
    logger.info(f"Found update times for {len(updated_repos)}/{len(repos)} repositories")
    logger.info(f"Found {len(unavailable_repos)} unavailable repositories")
    logger.info(f"Encountered errors with {len(error_repos)} repositories")
    logger.info(f"Found {len(github_repos_with_stars)} GitHub repositories with stars")
    
    # Group repositories by category
    repos_by_category = {}
    for data in results:
        # Handle both old and new format of data tuple
        if len(data) == 6:
            name, url, _, _, stars, _ = data
        else:
            name, url, _, _, stars = data
            
        if name in repo_categories and 'github.com' in url and stars is not None and stars > 0:
            category = repo_categories[name]
            if category not in repos_by_category:
                repos_by_category[category] = []
            repos_by_category[category].append((name, url, stars))
    
    # Log information about the top repositories by category
    logger.info("Top repositories by category:")
    categories_of_interest = [
        "Virus and Phage Identification",
        "Host Prediction",
        "Genome Analysis",
        "Taxonomy"
    ]
    
    for category in categories_of_interest:
        if category in repos_by_category:
            top_repos = sorted(repos_by_category[category], key=lambda x: x[2], reverse=True)[:3]
            logger.info(f"  {category}:")
            for name, _, stars in top_repos:
                logger.info(f"    - {name}: {stars} stars")
    
    # Update the README with the dates, availability status, and stars
    updated_content, unavailable, repos_with_stars = update_readme_with_dates_status_and_stars(readme_path, results)
    logger.info(f"README updated with repository information")
    logger.info(f"Updated both Popular Packages and Top Packages by Category sections")
    logger.info(f"Created unavailable_repos.md with {len(unavailable)} unavailable repositories")
    logger.info(f"Created starred_repos.md with {len(repos_with_stars)} starred repositories")
    
    # Save results to JSON for later reference
    with open(Path(__file__).parent / "repo_updates.json", 'w') as f:
        json_results = []
        for data in results:
            # Handle both old and new format of data tuple
            if len(data) == 6:
                name, url, date, status, stars, version_info = data
            else:
                name, url, date, status, stars = data
                version_info = None
                
            category = repo_categories.get(name, "Other")
            result_dict = {
                "name": name,
                "url": url,
                "category": category,
                "last_updated": date.isoformat() if date else None,
                "status": status,
                "stars": stars
            }
            
            # Add version info if available
            if version_info:
                version, year = version_info
                result_dict["version"] = version
                result_dict["version_year"] = year
                
            json_results.append(result_dict)
        json.dump(json_results, f, indent=2)
    
    logger.info("Results saved to repo_updates.json")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Check and update repository information in README")
    parser.add_argument("-r", "--readme", type=str, help="Path to the README file", default="README.md")
    args = parser.parse_args()

    main(args.readme)
