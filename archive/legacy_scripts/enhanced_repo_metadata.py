#!/usr/bin/env python3
"""
Enhanced Repository Metadata Collection Script

This script collects comprehensive repository-level metadata beyond basic stats:
- Citation information (CITATION.cff, DOIs, paper references)
- Documentation quality metrics
- Community health indicators
- Dependency analysis
- Installation methods

The script processes repositories from repo_updates.json, uses GitHub API 
with proper authentication and rate limiting, and stores results in a 
structured JSON format.

Usage:
    python enhanced_repo_metadata.py [--repo-file REPO_FILE] [--output-dir OUTPUT_DIR] [--incremental]
"""

import os
import re
import sys
import time
import json
import logging
import argparse
import datetime
import base64
import requests
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from urllib.parse import urlparse, quote

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("enhanced_metadata.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# GitHub API token (set as environment variable GITHUB_TOKEN)
GITHUB_TOKEN = os.environ.get("GITHUB_TOKEN", "")
if not GITHUB_TOKEN:
    logger.warning("GITHUB_TOKEN environment variable not set. Rate limiting will be more restrictive.")

# Define constants
MAX_WORKERS = 3  # Maximum number of concurrent requests
RATE_LIMIT_DELAY = 2.0  # Base delay between API requests (seconds)
MAX_RETRIES = 3  # Number of retries for API requests

# Documentation files to check for
DOC_FILES = [
    "README.md", "README", "readme.md", "readme",
    "CONTRIBUTING.md", "CONTRIBUTING", "contributing.md", "contributing",
    "CHANGELOG.md", "CHANGELOG", "changelog.md", "changelog",
    "DOCUMENTATION.md", "docs/README.md", "docs/index.md",
    "LICENSE", "LICENSE.md", "license", "license.md",
    "CODE_OF_CONDUCT.md", "SECURITY.md", "SUPPORT.md",
    "USAGE.md", "INSTALLATION.md"
]

# Dependency files to check for
DEP_FILES = [
    "requirements.txt", "pyproject.toml", "setup.py", "Pipfile", "environment.yml", 
    "package.json", "yarn.lock", "package-lock.json", "Gemfile", "Cargo.toml",
    "pom.xml", "build.gradle", "Podfile", "go.mod", "composer.json"
]

# Installation files to check for
INSTALL_FILES = [
    "setup.py", "pyproject.toml", "Makefile", "CMakeLists.txt", "configure",
    "Dockerfile", "docker-compose.yml", ".dockerignore", 
    "install.sh", "INSTALL", "INSTALL.md",
    "environment.yml", "Pipfile", "requirements.txt",
    "package.json"
]

# File patterns for each category (documentation, examples, tutorials)
DOC_PATTERNS = {
    "documentation": ["documentation", "docs", "doc", "manual", "guide", "handbook", "reference"],
    "examples": ["example", "examples", "demo", "demos", "sample", "samples"],
    "tutorials": ["tutorial", "tutorials", "howto", "how-to", "walkthrough", "step-by-step"],
    "tests": ["test", "tests", "specs", "spec"]
}

# Citation patterns for DOI and references
DOI_PATTERN = r'(doi|DOI):\s*10\.\d{4,}[\d\.]+/[^\s\]>),"]*'
ARXIV_PATTERN = r'arXiv:\s*\d{4}\.\d{4,}'
CITATION_PATTERNS = [
    # Citation or reference statements
    r'(?:cite|citing|citation|reference|refer to)(?:\s+this(?:\s+(?:code|software|tool|package|library|work|paper|article|project))?)?(?:\s+as)?\s*[:\-]?\s*(?P<citation>[^\.]{10,200}?\.\s*(?:\d{4}))',
    # Bibtext patterns
    r'@(?:article|book|inproceedings|proceedings|incollection|techreport|misc)\{(?P<citation>[^}]{10,300})\}',
    # Publication statements
    r'(?:published|appeared|presented)(?:\s+in)?\s+(?P<citation>[^\.]{10,200}?\.\s*(?:\d{4}))'
]

# Rate limiters
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

# Create rate limiters
github_limiter = RateLimiter(requests_per_minute=30)
gitlab_limiter = RateLimiter(requests_per_minute=30)
bitbucket_limiter = RateLimiter(requests_per_minute=30)

class GithubAPI:
    """GitHub API wrapper with authentication and rate limiting."""
    def __init__(self, token=None):
        self.token = token
        self.base_url = "https://api.github.com"
        self.headers = {
            "Accept": "application/vnd.github.v3+json"
        }
        if token:
            self.headers["Authorization"] = f"Bearer {token}"
        self.limiter = github_limiter

    def get(self, endpoint, params=None, retries=MAX_RETRIES):
        """Make a GET request to the GitHub API with rate limiting and retries."""
        self.limiter.wait()
        url = f"{self.base_url}/{endpoint.lstrip('/')}"
        
        for i in range(retries + 1):
            try:
                response = requests.get(url, headers=self.headers, params=params)
                
                # Check for rate limiting
                if response.status_code == 403 and 'X-RateLimit-Remaining' in response.headers:
                    remaining = int(response.headers['X-RateLimit-Remaining'])
                    reset_time = int(response.headers['X-RateLimit-Reset'])
                    reset_datetime = datetime.datetime.fromtimestamp(reset_time)
                    now = datetime.datetime.now()
                    wait_time = (reset_datetime - now).total_seconds()
                    
                    if remaining == 0 and wait_time > 0:
                        logger.warning(f"GitHub API rate limit exceeded. Reset at {reset_datetime}. Waiting {wait_time:.2f} seconds.")
                        if i < retries:
                            time.sleep(min(wait_time + 1, 3600))  # Wait until reset, but not more than an hour
                            continue
                
                if response.status_code == 404:
                    logger.warning(f"Resource not found: {url}")
                    return None
                
                elif response.status_code in [429, 500, 502, 503, 504]:
                    logger.warning(f"Temporary error {response.status_code} for {url}, retrying after delay")
                    if i < retries:
                        time.sleep(2 ** i)  # Exponential backoff
                        continue
                
                response.raise_for_status()
                return response.json()
            
            except requests.exceptions.RequestException as e:
                logger.error(f"Error accessing {url}: {e}")
                if i < retries:
                    time.sleep(2 ** i)  # Exponential backoff
                    continue
                return None
        
        return None

    def get_repository(self, owner, repo):
        """Get repository information."""
        return self.get(f"repos/{owner}/{repo}")
    
    def get_content(self, owner, repo, path, ref=None):
        """Get repository content."""
        params = {}
        if ref:
            params["ref"] = ref
        return self.get(f"repos/{owner}/{repo}/contents/{path}", params)
    
    def get_issues(self, owner, repo, state="all", per_page=100, page=1):
        """Get repository issues."""
        params = {"state": state, "per_page": per_page, "page": page}
        return self.get(f"repos/{owner}/{repo}/issues", params)
    
    def get_pulls(self, owner, repo, state="all", per_page=100, page=1):
        """Get repository pull requests."""
        params = {"state": state, "per_page": per_page, "page": page}
        return self.get(f"repos/{owner}/{repo}/pulls", params)
    
    def get_releases(self, owner, repo, per_page=100, page=1):
        """Get repository releases."""
        params = {"per_page": per_page, "page": page}
        return self.get(f"repos/{owner}/{repo}/releases", params)
    
    def get_contributors(self, owner, repo, per_page=100, page=1):
        """Get repository contributors."""
        params = {"per_page": per_page, "page": page}
        return self.get(f"repos/{owner}/{repo}/contributors", params)
    
    def get_file_content(self, owner, repo, file_path, ref=None):
        """Get file content from repository."""
        content_data = self.get_content(owner, repo, file_path, ref)
        if not content_data or not isinstance(content_data, dict) or "content" not in content_data:
            return None
        
        # Decode content from base64
        content = content_data.get("content", "")
        if content:
            try:
                return base64.b64decode(content.replace("\n", "")).decode('utf-8')
            except Exception as e:
                logger.error(f"Error decoding content for {owner}/{repo}/{file_path}: {e}")
        return None

def parse_github_url(url):
    """Extract owner and repo from GitHub URL."""
    patterns = [
        r"github\.com/([^/]+)/([^/]+)/?",  # https://github.com/owner/repo
        r"github\.com/([^/]+)/([^/]+)\.git"  # https://github.com/owner/repo.git
    ]
    
    for pattern in patterns:
        match = re.search(pattern, url)
        if match:
            owner, repo = match.groups()
            return owner, repo
    return None, None

def load_repo_data(repo_file_path):
    """Load repository data from JSON file."""
    try:
        with open(repo_file_path, 'r', encoding='utf-8') as f:
            return json.load(f)
    except (FileNotFoundError, json.JSONDecodeError) as e:
        logger.error(f"Error loading repo data: {e}")
        return []

def get_existing_metadata(output_dir):
    """Load existing metadata to support incremental updates."""
    metadata_dir = Path(output_dir)
    metadata_file = metadata_dir / "enhanced_repo_metadata.json"
    
    if metadata_file.exists():
        try:
            with open(metadata_file, 'r', encoding='utf-8') as f:
                return json.load(f)
        except (json.JSONDecodeError, Exception) as e:
            logger.error(f"Error loading existing metadata: {e}")
    
    return {"repositories": {}, "last_updated": "", "metadata_version": "1.0"}

def check_file_exists(github_api, owner, repo, file_path):
    """Check if a file exists in a repository."""
    content = github_api.get_content(owner, repo, file_path)
    return content is not None

def find_files_matching_pattern(github_api, owner, repo, patterns, search_dirs=None):
    """Find files matching any of the patterns in the repository."""
    if search_dirs is None:
        search_dirs = ["", "docs/", ".github/"]
    
    matches = []
    for directory in search_dirs:
        content = github_api.get_content(owner, repo, directory)
        if not content or not isinstance(content, list):
            continue
        
        for item in content:
            if item.get("type") == "file":
                file_name = item.get("name", "").lower()
                for pattern in patterns:
                    if pattern.lower() in file_name:
                        matches.append({
                            "name": item.get("name"),
                            "path": item.get("path"),
                            "url": item.get("html_url"),
                            "pattern": pattern
                        })
    
    return matches

def extract_citations(github_api, owner, repo):
    """Extract citation information from repository."""
    citations = {
        "citation_file": None,
        "doi": None,
        "references": [],
        "has_citation": False
    }
    
    # Check for CITATION.cff file
    citation_content = github_api.get_file_content(owner, repo, "CITATION.cff")
    if citation_content:
        citations["citation_file"] = "CITATION.cff"
        citations["has_citation"] = True
        
        # Extract DOI from CITATION.cff if present
        doi_match = re.search(DOI_PATTERN, citation_content)
        if doi_match:
            citations["doi"] = doi_match.group(0)
    
    # Check for DOI in README.md
    readme_content = github_api.get_file_content(owner, repo, "README.md")
    if readme_content:
        # Look for DOI pattern
        doi_match = re.search(DOI_PATTERN, readme_content)
        if doi_match and not citations["doi"]:
            citations["doi"] = doi_match.group(0)
            citations["has_citation"] = True
        
        # Look for arXiv pattern
        arxiv_match = re.search(ARXIV_PATTERN, readme_content)
        if arxiv_match:
            citations["references"].append({"type": "arxiv", "value": arxiv_match.group(0)})
            citations["has_citation"] = True
        
        # Look for citation patterns
        for pattern in CITATION_PATTERNS:
            matches = re.finditer(pattern, readme_content)
            for match in matches:
                if 'citation' in match.groupdict():
                    citation_text = match.group('citation').strip()
                    if citation_text and len(citation_text) > 10:  # Minimum length to be valid
                        citations["references"].append({"type": "text", "value": citation_text})
                        citations["has_citation"] = True
    
    # Check repository description for DOI
    repo_info = github_api.get_repository(owner, repo)
    if repo_info and "description" in repo_info and repo_info["description"]:
        doi_match = re.search(DOI_PATTERN, repo_info["description"])
        if doi_match and not citations["doi"]:
            citations["doi"] = doi_match.group(0)
            citations["has_citation"] = True
    
    return citations

def analyze_documentation(github_api, owner, repo):
    """Analyze documentation quality."""
    documentation = {
        "files": {},
        "completeness_score": 0,
        "has_examples": False,
        "has_tutorials": False,
        "has_user_guide": False,
        "total_doc_files": 0,
        "example_files": [],
        "tutorial_files": []
    }
    
    # Check for each documentation file
    for doc_file in DOC_FILES:
        file_exists = check_file_exists(github_api, owner, repo, doc_file)
        documentation["files"][doc_file] = file_exists
        if file_exists:
            documentation["total_doc_files"] += 1
    
    # Calculate documentation completeness score (0-100)
    essential_docs = ["README.md", "LICENSE", "CONTRIBUTING.md"]
    nice_to_have = ["CHANGELOG.md", "CODE_OF_CONDUCT.md", "SECURITY.md"]
    
    # Essential docs count for 60% of score
    essential_score = sum(1 for doc in essential_docs if documentation["files"].get(doc) or documentation["files"].get(doc.lower()))
    essential_percent = (essential_score / len(essential_docs)) * 60
    
    # Nice to have docs count for 40% of score
    nice_score = sum(1 for doc in nice_to_have if documentation["files"].get(doc) or documentation["files"].get(doc.lower()))
    nice_percent = (nice_score / len(nice_to_have)) * 40
    
    documentation["completeness_score"] = round(essential_percent + nice_percent)
    
    # Check for examples, tutorials, and user guides
    for category, patterns in DOC_PATTERNS.items():
        files = find_files_matching_pattern(github_api, owner, repo, patterns)
        if files:
            if category == "examples":
                documentation["has_examples"] = True
                documentation["example_files"] = files
            elif category == "tutorials":
                documentation["has_tutorials"] = True
                documentation["tutorial_files"] = files
            elif category == "documentation" and len(files) > 0:
                documentation["has_user_guide"] = True
    
    return documentation

def analyze_community_health(github_api, owner, repo):
    """Analyze community health metrics."""
    community = {
        "contributor_count": 0,
        "issue_resolution_time_days": None,
        "open_issues_count": 0,
        "closed_issues_count": 0,
        "open_prs_count": 0,
        "merged_prs_count": 0,
        "release_count": 0,
        "release_frequency_days": None,
        "last_activity_date": None,
        "activity_score": 0  # 0-100
    }
    
    # Get contributors
    contributors = github_api.get_contributors(owner, repo)
    if contributors:
        community["contributor_count"] = len(contributors)
    
    # Get issues
    closed_issues = github_api.get_issues(owner, repo, state="closed")
    open_issues = github_api.get_issues(owner, repo, state="open")
    
    if closed_issues:
        community["closed_issues_count"] = len(closed_issues)
        
        # Calculate average time to close issues
        resolution_times = []
        for issue in closed_issues:
            if "created_at" in issue and "closed_at" in issue:
                created = datetime.datetime.fromisoformat(issue["created_at"].replace("Z", "+00:00"))
                closed = datetime.datetime.fromisoformat(issue["closed_at"].replace("Z", "+00:00"))
                resolution_time = (closed - created).total_seconds() / (24 * 3600)  # in days
                resolution_times.append(resolution_time)
        
        if resolution_times:
            community["issue_resolution_time_days"] = round(sum(resolution_times) / len(resolution_times), 1)
    
    if open_issues:
        community["open_issues_count"] = len(open_issues)
    
    # Get pull requests
    open_prs = github_api.get_pulls(owner, repo, state="open")
    closed_prs = github_api.get_pulls(owner, repo, state="closed")
    merged_prs = [pr for pr in closed_prs if pr.get("merged_at")]
    
    if open_prs:
        community["open_prs_count"] = len(open_prs)
    if merged_prs:
        community["merged_prs_count"] = len(merged_prs)
    
    # Get releases
    releases = github_api.get_releases(owner, repo)
    if releases:
        community["release_count"] = len(releases)
        
        # Calculate release frequency
        if len(releases) >= 2:
            release_dates = []
            for release in releases:
                if "published_at" in release:
                    date = datetime.datetime.fromisoformat(release["published_at"].replace("Z", "+00:00"))
                    release_dates.append(date)
            
            release_dates.sort(reverse=True)
            time_diffs = [(release_dates[i] - release_dates[i+1]).total_seconds() / (24 * 3600)
                          for i in range(len(release_dates)-1)]
            
            if time_diffs:
                community["release_frequency_days"] = round(sum(time_diffs) / len(time_diffs), 1)
    
    # Determine last activity date
    activity_dates = []
    
    # Check latest commit date
    repo_info = github_api.get_repository(owner, repo)
    if repo_info and "pushed_at" in repo_info:
        pushed_date = datetime.datetime.fromisoformat(repo_info["pushed_at"].replace("Z", "+00:00"))
        activity_dates.append(pushed_date)
    
    # Check latest issue/PR activity
    latest_issues = github_api.get_issues(owner, repo, state="all", per_page=5)
    if latest_issues:
        for issue in latest_issues:
            if "updated_at" in issue:
                updated_date = datetime.datetime.fromisoformat(issue["updated_at"].replace("Z", "+00:00"))
                activity_dates.append(updated_date)
    
    if activity_dates:
        community["last_activity_date"] = max(activity_dates).isoformat()
    
    # Calculate activity score (0-100)
    activity_score = 0
    
    # Contributors: 0-25 points
    if community["contributor_count"] >= 10:
        activity_score += 25
    elif community["contributor_count"] >= 5:
        activity_score += 15
    elif community["contributor_count"] >= 2:
        activity_score += 10
    elif community["contributor_count"] >= 1:
        activity_score += 5
    
    # Issues and PRs: 0-25 points
    issue_pr_count = community["closed_issues_count"] + community["merged_prs_count"]
    if issue_pr_count >= 100:
        activity_score += 25
    elif issue_pr_count >= 50:
        activity_score += 20
    elif issue_pr_count >= 20:
        activity_score += 15
    elif issue_pr_count >= 10:
        activity_score += 10
    elif issue_pr_count >= 1:
        activity_score += 5
    
    # Releases: 0-25 points
    if community["release_count"] >= 10:
        activity_score += 25
    elif community["release_count"] >= 5:
        activity_score += 20
    elif community["release_count"] >= 3:
        activity_score += 15
    elif community["release_count"] >= 1:
        activity_score += 10
    
    # Recent activity: 0-25 points
    if community["last_activity_date"]:
        last_activity = datetime.datetime.fromisoformat(community["last_activity_date"].replace("Z", "+00:00"))
        days_since_activity = (datetime.datetime.now(datetime.timezone.utc) - last_activity).days
        
        if days_since_activity <= 7:
            activity_score += 25
        elif days_since_activity <= 30:
            activity_score += 20
        elif days_since_activity <= 90:
            activity_score += 15
        elif days_since_activity <= 180:
            activity_score += 10
        elif days_since_activity <= 365:
            activity_score += 5
    
    community["activity_score"] = activity_score
    
    return community

def analyze_dependencies(github_api, owner, repo):
    """Analyze dependency information."""
    dependencies = {
        "dependency_files": {},
        "has_docker_support": False,
        "platforms": [],
        "language_specific": {},
        "dependencies_found": False,
        "total_dependencies": 0
    }
    
    # Check for each dependency file
    for dep_file in DEP_FILES:
        content = github_api.get_file_content(owner, repo, dep_file)
        if content:
            dependencies["dependency_files"][dep_file] = True
            dependencies["dependencies_found"] = True
            
            # Parse specific dependency files
            if dep_file == "requirements.txt":
                python_deps = parse_requirements_txt(content)
                dependencies["language_specific"]["python"] = python_deps
                dependencies["total_dependencies"] += len(python_deps)
            
            elif dep_file == "package.json":
                node_deps = parse_package_json(content)
                dependencies["language_specific"]["javascript"] = node_deps
                dependencies["total_dependencies"] += len(node_deps.get("dependencies", [])) + len(node_deps.get("devDependencies", []))
            
            elif dep_file == "Pipfile":
                pipfile_deps = parse_pipfile(content)
                dependencies["language_specific"]["python_pipfile"] = pipfile_deps
                dependencies["total_dependencies"] += len(pipfile_deps)
            
            elif dep_file == "environment.yml":
                conda_deps = parse_environment_yml(content)
                dependencies["language_specific"]["conda"] = conda_deps
                dependencies["total_dependencies"] += len(conda_deps)
            
            elif dep_file == "Cargo.toml":
                rust_deps = parse_cargo_toml(content)
                dependencies["language_specific"]["rust"] = rust_deps
                dependencies["total_dependencies"] += len(rust_deps)
            
            # Add more parsers as needed for other dependency file formats
        else:
            dependencies["dependency_files"][dep_file] = False
    
    # Check for Docker support
    docker_files = ["Dockerfile", "docker-compose.yml", ".dockerignore"]
    for docker_file in docker_files:
        if github_api.get_file_content(owner, repo, docker_file):
            dependencies["has_docker_support"] = True
            break
    
    # Determine platform compatibility
    # This is a heuristic based on repository content
    repo_info = github_api.get_repository(owner, repo)
    
    if repo_info and "language" in repo_info:
        language = repo_info.get("language", "").lower()
        
        # Platform compatibility based on primary language
        if language in ["python", "ruby", "javascript", "typescript", "php", "go"]:
            dependencies["platforms"].extend(["linux", "macos", "windows"])
        elif language in ["c#", "vb.net"]:
            dependencies["platforms"].extend(["windows", "linux"])
        elif language in ["objective-c", "swift"]:
            dependencies["platforms"].append("macos")
        elif language in ["java", "kotlin"]:
            dependencies["platforms"].extend(["linux", "macos", "windows"])
        elif language in ["c", "c++"]:
            dependencies["platforms"].extend(["linux", "macos", "windows"])
        
    # Check for platform-specific keywords in files
    readme_content = github_api.get_file_content(owner, repo, "README.md")
    if readme_content:
        platforms = {
            "linux": ["linux", "ubuntu", "debian", "fedora", "centos", "redhat"],
            "macos": ["macos", "mac os", "osx", "mac", "apple"],
            "windows": ["windows", "win32", "win64"],
            "android": ["android"],
            "ios": ["ios", "iphone", "ipad"]
        }
        
        for platform, keywords in platforms.items():
            for keyword in keywords:
                if re.search(r'\b' + re.escape(keyword) + r'\b', readme_content, re.IGNORECASE):
                    if platform not in dependencies["platforms"]:
                        dependencies["platforms"].append(platform)
    
    # Deduplicate platforms
    dependencies["platforms"] = list(set(dependencies["platforms"]))
    
    return dependencies

def parse_requirements_txt(content):
    """Parse requirements.txt file to extract dependencies."""
    dependencies = []
    lines = content.split('\n')
    
    for line in lines:
        line = line.strip()
        # Skip comments and empty lines
        if not line or line.startswith('#'):
            continue
        
        # Handle basic requirements
        # Example: package==1.0.0, package>=1.0.0, package~=1.0.0
        match = re.match(r'^([A-Za-z0-9_\-\.]+)([=~<>!]+.+)?.*$', line)
        if match:
            package_name = match.group(1)
            version_spec = match.group(2) if match.group(2) else "any"
            dependencies.append({"name": package_name, "version": version_spec})
    
    return dependencies

def parse_package_json(content):
    """Parse package.json file to extract dependencies."""
    try:
        data = json.loads(content)
        return {
            "dependencies": data.get("dependencies", {}),
            "devDependencies": data.get("devDependencies", {}),
            "peerDependencies": data.get("peerDependencies", {})
        }
    except json.JSONDecodeError:
        return {"dependencies": {}, "devDependencies": {}, "peerDependencies": {}}

def parse_pipfile(content):
    """Parse Pipfile to extract dependencies."""
    dependencies = []
    in_packages_section = False
    
    lines = content.split('\n')
    for line in lines:
        line = line.strip()
        
        if line == "[packages]":
            in_packages_section = True
            continue
        elif line.startswith('[') and line.endswith(']'):
            in_packages_section = False
            continue
        
        if in_packages_section and "=" in line:
            parts = line.split('=', 1)
            if len(parts) == 2:
                package_name = parts[0].strip().strip('"')
                version = parts[1].strip().strip('"')
                dependencies.append({"name": package_name, "version": version})
    
    return dependencies

def parse_environment_yml(content):
    """Parse environment.yml to extract dependencies."""
    dependencies = []
    in_dependencies_section = False
    
    lines = content.split('\n')
    for line in lines:
        line = line.strip()
        
        if line == "dependencies:":
            in_dependencies_section = True
            continue
        elif line.startswith('name:') or line.startswith('channels:'):
            in_dependencies_section = False
            continue
        
        if in_dependencies_section and line.startswith('-'):
            dep = line[1:].strip()
            # Handle pip dependencies separately
            if dep == "pip:":
                continue
            if "=" in dep:
                parts = dep.split('=')
                if len(parts) >= 2:
                    package_name = parts[0].strip()
                    version = "=".join(parts[1:]).strip()
                    dependencies.append({"name": package_name, "version": version})
            else:
                dependencies.append({"name": dep, "version": "any"})
    
    return dependencies

def parse_cargo_toml(content):
    """Parse Cargo.toml to extract dependencies."""
    dependencies = []
    in_dependencies_section = False
    
    lines = content.split('\n')
    for line in lines:
        line = line.strip()
        
        if line == "[dependencies]":
            in_dependencies_section = True
            continue
        elif line.startswith('[') and line.endswith(']'):
            in_dependencies_section = False
            continue
        
        if in_dependencies_section and "=" in line:
            parts = line.split('=', 1)
            if len(parts) == 2:
                package_name = parts[0].strip()
                version = parts[1].strip().strip('"')
                dependencies.append({"name": package_name, "version": version})
    
    return dependencies

def analyze_installation_methods(github_api, owner, repo):
    """Analyze installation methods."""
    installation = {
        "methods": [],
        "files": {},
        "complexity_score": 0,  # 0-100, lower is simpler
        "has_install_docs": False
    }
    
    # Check for each installation file
    for install_file in INSTALL_FILES:
        content = github_api.get_file_content(owner, repo, install_file)
        installation["files"][install_file] = content is not None
        
        if content:
            # Identify installation methods based on files
            if install_file == "setup.py":
                installation["methods"].append("pip")
            elif install_file == "pyproject.toml":
                installation["methods"].append("pip")
                # Check if using poetry
                if "tool.poetry" in content:
                    installation["methods"].append("poetry")
            elif install_file == "environment.yml":
                installation["methods"].append("conda")
            elif install_file == "Pipfile":
                installation["methods"].append("pipenv")
            elif install_file == "requirements.txt":
                installation["methods"].append("pip")
            elif install_file == "package.json":
                installation["methods"].append("npm")
                if "yarn" in content:
                    installation["methods"].append("yarn")
            elif install_file == "Makefile":
                installation["methods"].append("make")
            elif install_file == "CMakeLists.txt":
                installation["methods"].append("cmake")
            elif install_file == "configure":
                installation["methods"].append("configure")
            elif install_file in ["Dockerfile", "docker-compose.yml"]:
                installation["methods"].append("docker")
    
    # Check README.md for installation instructions
    readme_content = github_api.get_file_content(owner, repo, "README.md")
    if readme_content:
        install_section_pattern = re.search(r'#+\s*Install(?:ation)?\s*(?:Instructions)?.*?\n(.*?)(?=#+|\Z)', 
                                          readme_content, re.DOTALL | re.IGNORECASE)
        
        if install_section_pattern:
            installation["has_install_docs"] = True
            install_section = install_section_pattern.group(1)
            
            # Check for various installation methods mentioned in the README
            install_methods = {
                "pip": [r'pip install', r'python -m pip install'],
                "conda": [r'conda install', r'conda env create'],
                "npm": [r'npm install', r'yarn add'],
                "yarn": [r'yarn add', r'yarn install'],
                "docker": [r'docker pull', r'docker run', r'docker-compose', r'docker build'],
                "make": [r'make install', r'make build'],
                "cmake": [r'cmake', r'cmake \.', r'cmake build'],
                "manual": [r'git clone', r'compile', r'build from source'],
                "homebrew": [r'brew install'],
                "apt": [r'apt install', r'apt-get install'],
                "yum": [r'yum install'],
                "pacman": [r'pacman -S']
            }
            
            for method, patterns in install_methods.items():
                for pattern in patterns:
                    if re.search(pattern, install_section, re.IGNORECASE):
                        if method not in installation["methods"]:
                            installation["methods"].append(method)
        
        # Otherwise check whole README for installation mentions
        else:
            install_keywords = [
                r'pip install', r'conda install', r'npm install', r'yarn add',
                r'docker pull', r'docker run', r'make install'
            ]
            
            for keyword in install_keywords:
                if re.search(keyword, readme_content, re.IGNORECASE):
                    installation["has_install_docs"] = True
                    break
    
    # Check for INSTALL.md or similar file
    install_files = ["INSTALL", "INSTALL.md", "INSTALLATION.md"]
    for install_file in install_files:
        if github_api.get_file_content(owner, repo, install_file):
            installation["has_install_docs"] = True
            break
    
    # Deduplicate methods
    installation["methods"] = list(set(installation["methods"]))
    
    # Calculate installation complexity score
    # Lower score means simpler installation
    complexity_score = 0
    
    # Method-based factors (0-50 points)
    if "pip" in installation["methods"] or "conda" in installation["methods"]:
        complexity_score += 10  # Very easy
    elif "npm" in installation["methods"] or "yarn" in installation["methods"]:
        complexity_score += 15  # Pretty easy
    elif "docker" in installation["methods"]:
        complexity_score += 20  # Moderately easy
    elif "make" in installation["methods"] or "cmake" in installation["methods"]:
        complexity_score += 30  # Moderate
    elif "manual" in installation["methods"]:
        complexity_score += 40  # Complex
    else:
        complexity_score += 50  # No clear method
    
    # Documentation factors (0-50 points)
    if installation["has_install_docs"]:
        if readme_content and len(readme_content) > 0:
            readme_install_ratio = len(re.findall(r'install', readme_content, re.IGNORECASE)) / len(readme_content.split())
            if readme_install_ratio > 0.02:  # More than 2% of words about installation
                complexity_score += 10  # Well documented
            else:
                complexity_score += 25  # Somewhat documented
        else:
            complexity_score += 30  # Basic documentation
    else:
        complexity_score += 50  # No documentation
    
    installation["complexity_score"] = min(100, complexity_score)
    
    return installation

def analyze_github_repo(repo_data, github_api, existing_metadata=None):
    """Analyze a GitHub repository to collect enhanced metadata."""
    repo_name = repo_data.get("name", "")
    repo_url = repo_data.get("url", "")
    
    # Extract owner and repo from GitHub URL
    owner, repo = parse_github_url(repo_url)
    if not owner or not repo:
        logger.warning(f"Invalid GitHub URL: {repo_url}")
        return None
    
    logger.info(f"Analyzing GitHub repository: {owner}/{repo}")
    
    # Check if we need to update this repository's metadata
    should_update = True
    if existing_metadata and "repositories" in existing_metadata:
        if repo_url in existing_metadata["repositories"]:
            existing = existing_metadata["repositories"][repo_url]
            # Check if data is less than 30 days old
            if "fetch_time" in existing:
                fetch_time = datetime.datetime.fromisoformat(existing["fetch_time"])
                days_since_fetch = (datetime.datetime.now() - fetch_time).days
                if days_since_fetch < 30:
                    logger.info(f"Using existing metadata for {repo_url} (fetched {days_since_fetch} days ago)")
                    return existing
    
    try:
        # Get base repository information
        repo_info = github_api.get_repository(owner, repo)
        if not repo_info:
            logger.warning(f"Repository not found: {repo_url}")
            return None
        
        # Collect enhanced metadata
        # 1. Citation information
        logger.info(f"Extracting citation information for {owner}/{repo}")
        citation_info = extract_citations(github_api, owner, repo)
        
        # 2. Documentation quality
        logger.info(f"Analyzing documentation for {owner}/{repo}")
        documentation_info = analyze_documentation(github_api, owner, repo)
        
        # 3. Community health
        logger.info(f"Analyzing community health for {owner}/{repo}")
        community_info = analyze_community_health(github_api, owner, repo)
        
        # 4. Dependency analysis
        logger.info(f"Analyzing dependencies for {owner}/{repo}")
        dependency_info = analyze_dependencies(github_api, owner, repo)
        
        # 5. Installation methods
        logger.info(f"Analyzing installation methods for {owner}/{repo}")
        installation_info = analyze_installation_methods(github_api, owner, repo)
        
        # Combine all metadata
        metadata = {
            "name": repo_name,
            "url": repo_url,
            "owner": owner,
            "repo": repo,
            "fetch_time": datetime.datetime.now().isoformat(),
            "citation": citation_info,
            "documentation": documentation_info,
            "community": community_info,
            "dependencies": dependency_info,
            "installation": installation_info,
            "metadata_version": "1.0"
        }
        
        return metadata
    
    except Exception as e:
        logger.error(f"Error analyzing repository {repo_url}: {e}")
        return None

def process_repositories(repo_data, output_dir, incremental=False):
    """Process repositories to collect enhanced metadata."""
    # Create output directory if it doesn't exist
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    # Load existing metadata if incremental update
    existing_metadata = None
    if incremental:
        existing_metadata = get_existing_metadata(output_dir)
    
    # Initialize GitHub API
    github_api = GithubAPI(token=GITHUB_TOKEN)
    
    # Process repositories
    results = {
        "repositories": {},
        "last_updated": datetime.datetime.now().isoformat(),
        "metadata_version": "1.0",
        "total_repositories": len(repo_data),
        "successful_repositories": 0,
        "failed_repositories": 0
    }
    
    # Group repositories by platform
    repositories = {
        "github": [],
        "gitlab": [],
        "bitbucket": [],
        "other": []
    }
    
    for repo in repo_data:
        url = repo.get("url", "").lower()
        if "github.com" in url:
            repositories["github"].append(repo)
        elif "gitlab.com" in url:
            repositories["gitlab"].append(repo)
        elif "bitbucket.org" in url:
            repositories["bitbucket"].append(repo)
        else:
            repositories["other"].append(repo)
    
    # Process GitHub repositories (currently only GitHub is supported)
    logger.info(f"Processing {len(repositories['github'])} GitHub repositories")
    
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {}
        for repo in repositories["github"]:
            future = executor.submit(analyze_github_repo, repo, github_api, existing_metadata)
            futures[future] = repo
        
        for future in as_completed(futures):
            repo = futures[future]
            repo_url = repo.get("url", "")
            
            try:
                metadata = future.result()
                if metadata:
                    results["repositories"][repo_url] = metadata
                    results["successful_repositories"] += 1
                    logger.info(f"Successfully processed {repo_url}")
                else:
                    results["failed_repositories"] += 1
                    logger.warning(f"Failed to process {repo_url}")
            except Exception as e:
                results["failed_repositories"] += 1
                logger.error(f"Error processing {repo_url}: {e}")
    
    # Save results
    output_file = output_path / "enhanced_repo_metadata.json"
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(results, f, indent=2)
    
    logger.info(f"Enhanced metadata saved to {output_file}")
    logger.info(f"Successfully processed {results['successful_repositories']}/{results['total_repositories']} repositories")
    
    # Create metadata summary
    create_metadata_summary(results, output_path)
    
    return results

def create_metadata_summary(results, output_path):
    """Create a summary of the enhanced metadata."""
    summary = {
        "last_updated": results["last_updated"],
        "total_repositories": results["total_repositories"],
        "successful_repositories": results["successful_repositories"],
        "repository_summaries": [],
        "citation_stats": {
            "has_citation": 0,
            "has_doi": 0
        },
        "documentation_stats": {
            "avg_completeness_score": 0,
            "has_examples": 0,
            "has_tutorials": 0,
            "has_user_guide": 0
        },
        "community_stats": {
            "avg_contributor_count": 0,
            "avg_activity_score": 0
        },
        "dependency_stats": {
            "has_docker_support": 0,
            "platforms": {
                "linux": 0,
                "macos": 0,
                "windows": 0,
                "android": 0,
                "ios": 0
            }
        },
        "installation_stats": {
            "methods": {},
            "avg_complexity_score": 0
        }
    }
    
    if results["successful_repositories"] == 0:
        # Save empty summary
        output_file = output_path / "enhanced_metadata_summary.json"
        with open(output_file, 'w', encoding='utf-8') as f:
            json.dump(summary, f, indent=2)
        return
    
    # Collect statistics
    completeness_scores = []
    activity_scores = []
    contributor_counts = []
    complexity_scores = []
    
    for repo_url, metadata in results["repositories"].items():
        # Create repository summary
        repo_summary = {
            "name": metadata.get("name", ""),
            "url": repo_url,
            "citation": {
                "has_citation": metadata.get("citation", {}).get("has_citation", False),
                "doi": metadata.get("citation", {}).get("doi")
            },
            "documentation": {
                "completeness_score": metadata.get("documentation", {}).get("completeness_score", 0),
                "has_examples": metadata.get("documentation", {}).get("has_examples", False),
                "has_tutorials": metadata.get("documentation", {}).get("has_tutorials", False),
                "has_user_guide": metadata.get("documentation", {}).get("has_user_guide", False)
            },
            "community": {
                "contributor_count": metadata.get("community", {}).get("contributor_count", 0),
                "activity_score": metadata.get("community", {}).get("activity_score", 0),
                "last_activity_date": metadata.get("community", {}).get("last_activity_date")
            },
            "dependencies": {
                "has_docker_support": metadata.get("dependencies", {}).get("has_docker_support", False),
                "platforms": metadata.get("dependencies", {}).get("platforms", [])
            },
            "installation": {
                "methods": metadata.get("installation", {}).get("methods", []),
                "complexity_score": metadata.get("installation", {}).get("complexity_score", 0)
            }
        }
        
        summary["repository_summaries"].append(repo_summary)
        
        # Update citation stats
        if metadata.get("citation", {}).get("has_citation", False):
            summary["citation_stats"]["has_citation"] += 1
            if metadata.get("citation", {}).get("doi"):
                summary["citation_stats"]["has_doi"] += 1
        
        # Update documentation stats
        doc_score = metadata.get("documentation", {}).get("completeness_score", 0)
        completeness_scores.append(doc_score)
        if metadata.get("documentation", {}).get("has_examples", False):
            summary["documentation_stats"]["has_examples"] += 1
        if metadata.get("documentation", {}).get("has_tutorials", False):
            summary["documentation_stats"]["has_tutorials"] += 1
        if metadata.get("documentation", {}).get("has_user_guide", False):
            summary["documentation_stats"]["has_user_guide"] += 1
        
        # Update community stats
        contributor_count = metadata.get("community", {}).get("contributor_count", 0)
        contributor_counts.append(contributor_count)
        activity_score = metadata.get("community", {}).get("activity_score", 0)
        activity_scores.append(activity_score)
        
        # Update dependency stats
        if metadata.get("dependencies", {}).get("has_docker_support", False):
            summary["dependency_stats"]["has_docker_support"] += 1
        platforms = metadata.get("dependencies", {}).get("platforms", [])
        for platform in platforms:
            if platform in summary["dependency_stats"]["platforms"]:
                summary["dependency_stats"]["platforms"][platform] += 1
        
        # Update installation stats
        methods = metadata.get("installation", {}).get("methods", [])
        for method in methods:
            if method not in summary["installation_stats"]["methods"]:
                summary["installation_stats"]["methods"][method] = 0
            summary["installation_stats"]["methods"][method] += 1
        complexity_score = metadata.get("installation", {}).get("complexity_score", 0)
        complexity_scores.append(complexity_score)
    
    # Calculate averages
    if completeness_scores:
        summary["documentation_stats"]["avg_completeness_score"] = round(sum(completeness_scores) / len(completeness_scores), 1)
    if activity_scores:
        summary["community_stats"]["avg_activity_score"] = round(sum(activity_scores) / len(activity_scores), 1)
    if contributor_counts:
        summary["community_stats"]["avg_contributor_count"] = round(sum(contributor_counts) / len(contributor_counts), 1)
    if complexity_scores:
        summary["installation_stats"]["avg_complexity_score"] = round(sum(complexity_scores) / len(complexity_scores), 1)
    
    # Sort repository summaries by activity score (highest first)
    summary["repository_summaries"].sort(key=lambda x: x["community"]["activity_score"], reverse=True)
    
    # Save summary
    output_file = output_path / "enhanced_metadata_summary.json"
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump(summary, f, indent=2)
    
    logger.info(f"Enhanced metadata summary saved to {output_file}")

def main():
    """Main function to collect enhanced repository metadata."""
    parser = argparse.ArgumentParser(description="Collect enhanced repository metadata.")
    parser.add_argument("--repo-file", default="repo_updates.json", help="Path to repository updates JSON file")
    parser.add_argument("--output-dir", default="metadata", help="Directory to store metadata results")
    parser.add_argument("--incremental", action="store_true", help="Perform incremental update instead of full scan")
    args = parser.parse_args()
    
    # Load repository data
    logger.info(f"Loading repository data from {args.repo_file}")
    repo_data = load_repo_data(args.repo_file)
    if not repo_data:
        logger.error("No repository data found or file is empty")
        sys.exit(1)
    
    logger.info(f"Found {len(repo_data)} repositories to process")
    
    # Process repositories
    results = process_repositories(repo_data, args.output_dir, args.incremental)
    
    logger.info(f"Enhanced metadata collection completed for {results['successful_repositories']}/{results['total_repositories']} repositories")

if __name__ == "__main__":
    main()