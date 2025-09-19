#!/usr/bin/env python3
"""
GitHub Metrics Enhancer - Phase 1 Citation System V2

A focused, reliable system for collecting GitHub repository metrics.
Uses only free GitHub API - no premium API dependencies.

Features:
- Collects comprehensive GitHub repository data
- Handles rate limiting gracefully
- Updates existing data.json with enhanced metrics
- Safe error handling and logging
"""

import json
import requests
import time
import re
from typing import Dict, List, Optional
from datetime import datetime
import os
import sys

class GitHubMetricsCollector:
    """Collects comprehensive GitHub metrics for repository tools."""

    def __init__(self, github_token: Optional[str] = None):
        """Initialize with optional GitHub token for higher rate limits."""
        self.github_token = github_token
        self.session = requests.Session()
        if github_token:
            self.session.headers.update({
                'Authorization': f'token {github_token}',
                'Accept': 'application/vnd.github.v3+json'
            })

        self.base_url = "https://api.github.com"
        self.rate_limit_remaining = 60  # Default for unauthenticated

    def extract_github_repo(self, url: str) -> Optional[str]:
        """Extract owner/repo from GitHub URL."""
        patterns = [
            r'github\.com/([^/]+/[^/]+)',
            r'github\.com/([^/]+/[^/]+)\.git',
            r'github\.com/([^/]+/[^/]+)/.*'
        ]

        for pattern in patterns:
            match = re.search(pattern, url)
            if match:
                repo_path = match.group(1)
                # Remove any trailing path components
                if '/' in repo_path[repo_path.find('/') + 1:]:
                    repo_path = '/'.join(repo_path.split('/')[:2])
                return repo_path
        return None

    def check_rate_limit(self):
        """Check GitHub API rate limit status."""
        try:
            response = self.session.get(f"{self.base_url}/rate_limit")
            if response.status_code == 200:
                data = response.json()
                self.rate_limit_remaining = data['rate']['remaining']
                reset_time = data['rate']['reset']
                return {
                    'remaining': self.rate_limit_remaining,
                    'reset_time': reset_time,
                    'limit': data['rate']['limit']
                }
        except Exception as e:
            print(f"Warning: Could not check rate limit: {e}")
        return None

    def get_repo_metrics(self, repo_path: str) -> Dict:
        """Get comprehensive repository metrics from GitHub API."""
        if self.rate_limit_remaining <= 1:
            print("Rate limit approaching, waiting...")
            time.sleep(60)
            self.check_rate_limit()

        try:
            # Get repository data
            response = self.session.get(f"{self.base_url}/repos/{repo_path}")
            if response.status_code != 200:
                return {'error': f"HTTP {response.status_code}"}

            repo_data = response.json()

            # Extract comprehensive metrics
            metrics = {
                'stars': repo_data.get('stargazers_count', 0),
                'forks': repo_data.get('forks_count', 0),
                'watchers': repo_data.get('watchers_count', 0),
                'open_issues': repo_data.get('open_issues_count', 0),
                'size_kb': repo_data.get('size', 0),
                'language': repo_data.get('language'),
                'license': repo_data.get('license', {}).get('name') if repo_data.get('license') else None,
                'created_at': repo_data.get('created_at'),
                'updated_at': repo_data.get('updated_at'),
                'pushed_at': repo_data.get('pushed_at'),
                'archived': repo_data.get('archived', False),
                'disabled': repo_data.get('disabled', False),
                'topics': repo_data.get('topics', []),
                'description': repo_data.get('description', ''),
                'default_branch': repo_data.get('default_branch', 'main'),
                'has_issues': repo_data.get('has_issues', False),
                'has_wiki': repo_data.get('has_wiki', False),
                'has_pages': repo_data.get('has_pages', False)
            }

            # Get language breakdown
            try:
                lang_response = self.session.get(f"{self.base_url}/repos/{repo_path}/languages")
                if lang_response.status_code == 200:
                    metrics['languages'] = lang_response.json()
                    metrics['all_languages'] = list(metrics['languages'].keys())
            except Exception as e:
                print(f"Warning: Could not get languages for {repo_path}: {e}")

            # Get recent activity metrics
            try:
                # Get recent commits (last 30 days activity indicator)
                since_date = datetime.now().replace(day=1).isoformat() + 'Z'
                commits_response = self.session.get(
                    f"{self.base_url}/repos/{repo_path}/commits",
                    params={'since': since_date, 'per_page': 1}
                )
                if commits_response.status_code == 200:
                    metrics['recent_activity'] = len(commits_response.json()) > 0
                else:
                    metrics['recent_activity'] = None
            except Exception as e:
                print(f"Warning: Could not get activity for {repo_path}: {e}")
                metrics['recent_activity'] = None

            # Update rate limit counter
            self.rate_limit_remaining -= 2  # We made 2+ requests

            return metrics

        except Exception as e:
            return {'error': str(e)}

    def enhance_tool_metrics(self, tool: Dict) -> Dict:
        """Enhance a single tool's GitHub metrics."""
        url = tool.get('url', '')
        if 'github.com' not in url:
            return tool  # Skip non-GitHub tools

        repo_path = self.extract_github_repo(url)
        if not repo_path:
            print(f"Warning: Could not extract repo path from {url}")
            return tool

        print(f"Enhancing metrics for {repo_path}...")
        metrics = self.get_repo_metrics(repo_path)

        if 'error' in metrics:
            print(f"Error getting metrics for {repo_path}: {metrics['error']}")
            return tool

        # Update tool with enhanced metrics
        enhanced_tool = tool.copy()
        enhanced_tool.update({
            'github_metrics': metrics,
            'repo_path': repo_path,
            'provider': 'github',
            'last_metrics_update': datetime.now().isoformat(),
            # Update existing fields with fresh data
            'stars': metrics.get('stars', tool.get('stars', 0)),
            'forks': metrics.get('forks', tool.get('forks', 0)),
            'language': metrics.get('language', tool.get('language')),
            'is_archived': metrics.get('archived', False),
            'languages': metrics.get('languages', {}),
            'all_languages': metrics.get('all_languages', [])
        })

        return enhanced_tool

def main():
    """Main function to enhance GitHub metrics in data.json."""

    # Check for GitHub token
    github_token = os.environ.get('GITHUB_TOKEN')
    if not github_token:
        print("Warning: No GITHUB_TOKEN found. Using unauthenticated requests (lower rate limit).")
        print("Set GITHUB_TOKEN environment variable for higher rate limits.")

    # Initialize collector
    collector = GitHubMetricsCollector(github_token)

    # Check rate limit
    rate_info = collector.check_rate_limit()
    if rate_info:
        print(f"GitHub API Rate Limit: {rate_info['remaining']}/{rate_info['limit']} remaining")

    # Load existing data
    try:
        with open('data.json', 'r') as f:
            data = json.load(f)
    except Exception as e:
        print(f"Error loading data.json: {e}")
        return 1

    # Find GitHub tools to enhance
    github_tools = []
    for node in data.get('nodes', []):
        if node.get('type') == 'tool' and 'github.com' in node.get('url', ''):
            github_tools.append(node)

    print(f"Found {len(github_tools)} GitHub tools to enhance")

    # Enhance metrics for each tool
    enhanced_count = 0
    for i, tool in enumerate(github_tools):
        print(f"Processing {i+1}/{len(github_tools)}: {tool['name']}")

        enhanced_tool = collector.enhance_tool_metrics(tool)
        if 'github_metrics' in enhanced_tool:
            enhanced_count += 1
            # Update the tool in the data structure
            for j, node in enumerate(data['nodes']):
                if node.get('id') == tool['id']:
                    data['nodes'][j] = enhanced_tool
                    break

        # Rate limiting pause
        if i > 0 and i % 10 == 0:
            print("Pausing for rate limiting...")
            time.sleep(2)

    # Save enhanced data
    try:
        with open('data.json', 'w') as f:
            json.dump(data, f, indent=2, ensure_ascii=False)
        print(f"Successfully enhanced {enhanced_count} tools")
        return 0
    except Exception as e:
        print(f"Error saving data.json: {e}")
        return 1

if __name__ == "__main__":
    sys.exit(main())