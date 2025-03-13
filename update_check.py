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
    with open(readme_path, 'r') as f:
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
    """Get the last updated time and star count for a GitHub repository."""
    try:
        # Use rate limiter
        github_limiter.wait()
        
        # Extract owner/repo from URL
        match = re.search(r'github\.com/([^/]+/[^/]+)', repo_url)
        if not match:
            return None, None, None
        
        repo_path = match.group(1).rstrip('/')
        
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
            return None, "not_found", None
        
        if response.status_code == 429:
            logger.warning(f"Rate limit hit for {repo_url}, waiting and retrying")
            time.sleep(10)  # Wait 10 seconds before retrying
            return get_github_repo_info(repo_url)  # Retry
            
        response.raise_for_status()
        
        data = response.json()
        updated_at = data.get('updated_at') or data.get('pushed_at')
        stars = data.get('stargazers_count', 0)
        
        if updated_at:
            return datetime.strptime(updated_at, '%Y-%m-%dT%H:%M:%SZ'), None, stars
        return None, None, stars
    except requests.exceptions.RequestException as e:
        if hasattr(e, 'response') and e.response is not None:
            if e.response.status_code == 404:
                logger.warning(f"Repository not found: {repo_url}")
                return None, "not_found", None
            elif e.response.status_code == 429:
                logger.warning(f"Rate limit hit for {repo_url}, waiting 60 seconds")
                time.sleep(60)  # Wait 60 seconds before retrying
                return get_github_repo_info(repo_url)  # Retry
        logger.error(f"Error fetching GitHub repo {repo_url}: {e}")
        return None, "error", None
    except Exception as e:
        logger.error(f"Unexpected error for {repo_url}: {e}")
        return None, "error", None

def get_gitlab_repo_info(repo_url):
    """Get the last updated time and star count for a GitLab repository."""
    try:
        # Use rate limiter
        gitlab_limiter.wait()
        
        # Extract project ID or path from URL
        match = re.search(r'gitlab\.com/([^/]+/[^/]+)', repo_url)
        if not match:
            return None, None, None
        
        project_path = match.group(1).rstrip('/')
        
        # URL encode the project path
        encoded_path = requests.utils.quote(project_path, safe='')
        
        # Create API URL for project lookup
        api_url = f"{GITLAB_API}?search={encoded_path}"
        
        logger.info(f"Fetching GitLab repo: {project_path}")
        response = requests.get(api_url)
        
        if response.status_code == 404:
            logger.warning(f"Repository not found: {repo_url}")
            return None, "not_found", None
        
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
                
                if updated_at:
                    return datetime.strptime(updated_at, '%Y-%m-%dT%H:%M:%S.%fZ'), None, stars
                return None, None, stars
        return None, None, None
    except requests.exceptions.RequestException as e:
        if hasattr(e, 'response') and e.response is not None:
            if e.response.status_code == 404:
                logger.warning(f"Repository not found: {repo_url}")
                return None, "not_found", None
            elif e.response.status_code == 429:
                logger.warning(f"Rate limit hit for {repo_url}, waiting 60 seconds")
                time.sleep(60)  # Wait 60 seconds before retrying
                return get_gitlab_repo_info(repo_url)  # Retry
        logger.error(f"Error fetching GitLab repo {repo_url}: {e}")
        return None, "error", None
    except Exception as e:
        logger.error(f"Unexpected error for {repo_url}: {e}")
        return None, "error", None

def get_bitbucket_repo_info(repo_url):
    """Get the last updated time and watchers count for a Bitbucket repository."""
    try:
        # Use rate limiter
        bitbucket_limiter.wait()
        
        # Extract workspace/repo from URL
        match = re.search(r'bitbucket\.org/([^/]+/[^/]+)', repo_url)
        if not match:
            return None, None, None
        
        repo_path = match.group(1).rstrip('/')
        
        # Create API URL
        api_url = f"{BITBUCKET_API}{repo_path}"
        
        logger.info(f"Fetching Bitbucket repo: {repo_path}")
        response = requests.get(api_url)
        
        if response.status_code == 404:
            logger.warning(f"Repository not found: {repo_url}")
            return None, "not_found", None
        
        if response.status_code == 429:
            logger.warning(f"Rate limit hit for {repo_url}, waiting 60 seconds")
            time.sleep(60)  # Wait 60 seconds before retrying
            return get_bitbucket_repo_info(repo_url)  # Retry
            
        response.raise_for_status()
        
        data = response.json()
        updated_on = data.get('updated_on')
        # Bitbucket doesn't have stars but has watchers
        watchers = data.get('watchers_count', 0)
        
        if updated_on:
            return datetime.strptime(updated_on, '%Y-%m-%dT%H:%M:%S.%f%z'), None, watchers
        return None, None, watchers
    except requests.exceptions.RequestException as e:
        if hasattr(e, 'response') and e.response is not None:
            if e.response.status_code == 404:
                logger.warning(f"Repository not found: {repo_url}")
                return None, "not_found", None
            elif e.response.status_code == 429:
                logger.warning(f"Rate limit hit for {repo_url}, waiting 60 seconds")
                time.sleep(60)  # Wait 60 seconds before retrying
                return get_bitbucket_repo_info(repo_url)  # Retry
        logger.error(f"Error fetching Bitbucket repo {repo_url}: {e}")
        return None, "error", None
    except Exception as e:
        logger.error(f"Unexpected error for {repo_url}: {e}")
        return None, "error", None

def get_repo_last_updated(repo_name, repo_url):
    """Get the last updated time and stars for a repository based on its URL."""
    # Handle different repository hosts
    if 'github.com' in repo_url:
        last_updated, status, stars = get_github_repo_info(repo_url)
        return repo_name, repo_url, last_updated, status, stars
    elif 'gitlab.com' in repo_url:
        last_updated, status, stars = get_gitlab_repo_info(repo_url)
        return repo_name, repo_url, last_updated, status, stars
    elif 'bitbucket.org' in repo_url:
        last_updated, status, stars = get_bitbucket_repo_info(repo_url)
        return repo_name, repo_url, last_updated, status, stars
    return repo_name, repo_url, None, None, None

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
                    batch_results.append((name, url, None, "error", None))
        
        results.extend(batch_results)
        
        # If there are more batches to process, wait to avoid rate limiting
        if i + batch_size < total_repos:
            logger.info(f"Waiting {batch_delay} seconds before processing next batch...")
            time.sleep(batch_delay)
    
    return results

def update_readme_with_dates_status_and_stars(readme_path, repo_data):
    """Update the README with the last updated information, availability status, and create a popular packages section."""
    with open(readme_path, 'r') as f:
        content = f.read()
    
    updated_content = content
    unavailable_repos = []
    github_repos_with_stars = []
    
    # Sort repositories by name for consistent processing
    sorted_repos = sorted(repo_data, key=lambda x: x[0])
    
    for repo_name, repo_url, last_updated, status, stars in sorted_repos:
        # Collect GitHub repos with stars for popular packages section
        if 'github.com' in repo_url and stars is not None and stars > 0:
            github_repos_with_stars.append((repo_name, repo_url, stars))
            
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
                updated_content = updated_content.replace(
                    f"[{repo_name}]({repo_url})",
                    f"[{repo_name}]({repo_url}) [unavailable]"
                )
        
        # Add update date if available
        elif last_updated:
            year = last_updated.year
            month = last_updated.month
            
            # Format the date as [Updated: MM/YYYY]
            date_info = f"[Updated: {month:02d}/{year}]"
            
            # Check if the repo is already in the README with a date
            repo_pattern = re.escape(f"[{repo_name}]({repo_url})") + r".*?(\[Updated:.*?\])?"
            date_match = re.search(repo_pattern, updated_content)
            
            if date_match and date_match.group(1):
                # Replace existing date
                updated_content = updated_content.replace(
                    date_match.group(1),
                    date_info
                )
            else:
                # Add new date after the URL
                updated_content = updated_content.replace(
                    f"[{repo_name}]({repo_url})",
                    f"[{repo_name}]({repo_url}) {date_info}"
                )
    
    # Create popular packages section
    if github_repos_with_stars:
        top_repos = sorted(github_repos_with_stars, key=lambda x: x[2], reverse=True)[:20]  # Top 20 repos by stars
        
        popular_packages_section = "## Popular Packages\n\nRanked by GitHub stars:\n\n"
        for i, (name, url, stars) in enumerate(top_repos, 1):
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
    
    # Write the updated content back to the README
    with open(readme_path, 'w') as f:
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
            f.write("| Repository | Stars |\n")
            f.write("|------------|-------|\n")
            for name, url, stars in sorted_by_stars:
                f.write(f"| [{name}]({url}) | {stars} |\n")
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
    
    # Update the README with the dates, availability status, and stars
    updated_content, unavailable, repos_with_stars = update_readme_with_dates_status_and_stars(readme_path, results)
    logger.info(f"README updated with repository information")
    logger.info(f"Created unavailable_repos.md with {len(unavailable)} unavailable repositories")
    logger.info(f"Created starred_repos.md with {len(repos_with_stars)} starred repositories")
    
    # Save results to JSON for later reference
    with open(Path(__file__).parent / "repo_updates.json", 'w') as f:
        json_results = []
        for name, url, date, status, stars in results:
            json_results.append({
                "name": name,
                "url": url,
                "last_updated": date.isoformat() if date else None,
                "status": status,
                "stars": stars
            })
        json.dump(json_results, f, indent=2)
    
    logger.info("Results saved to repo_updates.json")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Check and update repository information in README")
    parser.add_argument("-r", "--readme", type=str, help="Path to the README file", default="README.md")
    args = parser.parse_args()

    main(args.readme)
