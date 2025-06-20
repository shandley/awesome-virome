#!/usr/bin/env python3
"""
Process Citation Validation Results

This script processes the results of citation validation and creates
a pull request with any fixes or updated metrics.

It handles:
1. Checking validation results
2. Updating metrics files
3. Committing changes
4. Creating a PR with appropriate description

This script is designed to be called from the GitHub Actions workflow.
"""

import os
import sys
import json
import datetime
import subprocess
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple
from github import Github

# Constants
REPO_ROOT = Path(__file__).parent.parent
REPORTS_DIR = REPO_ROOT / "reports"
VALIDATION_REPORT_PATH = REPORTS_DIR / "citation_validation.json"
METRICS_DIR = REPORTS_DIR / "citations"
METRICS_FILE_PATH = METRICS_DIR / "citation_validation_metrics.json"

def run_cmd(cmd: str) -> Tuple[int, str, str]:
    """Run a shell command and return exit code, stdout, and stderr"""
    process = subprocess.Popen(
        cmd, 
        shell=True, 
        stdout=subprocess.PIPE, 
        stderr=subprocess.PIPE,
        universal_newlines=True
    )
    stdout, stderr = process.communicate()
    return process.returncode, stdout, stderr

def setup_git_config() -> None:
    """Set up git configuration for GitHub Actions"""
    print("Setting up git configuration...")
    run_cmd('git config --local user.email "action@github.com"')
    run_cmd('git config --local user.name "GitHub Action"')

def load_validation_results() -> Dict[str, Any]:
    """Load validation results from JSON file"""
    print(f"Loading validation results from {VALIDATION_REPORT_PATH}...")
    
    if not VALIDATION_REPORT_PATH.exists():
        print("Validation report not found!")
        return {}
    
    try:
        with open(VALIDATION_REPORT_PATH, 'r') as f:
            return json.load(f)
    except json.JSONDecodeError:
        print("Invalid JSON in validation report!")
        return {}

def update_metrics_file(validation_results: Dict[str, Any]) -> bool:
    """Update metrics file with validation results"""
    print(f"Updating metrics file at {METRICS_FILE_PATH}...")
    
    # Create directory if it doesn't exist
    METRICS_DIR.mkdir(parents=True, exist_ok=True)
    
    if not validation_results:
        # Create an empty metrics file
        empty_metrics = {
            "timestamp": datetime.datetime.now().isoformat(),
            "total_tools": 0,
            "tools_with_doi": 0,
            "doi_percentage": 0,
            "valid_doi_percentage": 0,
            "consistent_doi_percentage": 0,
            "tools_with_pubmed_data": 0,
            "pubmed_percentage": 0,
            "tools_with_academic_impact": 0,
            "academic_impact_percentage": 0,
            "citation_format_counts": {"bibtex": 0, "apa": 0, "mla": 0},
            "issues_count": 1,
            "error": "No validation data available"
        }
        
        with open(METRICS_FILE_PATH, 'w') as f:
            json.dump(empty_metrics, f, indent=2)
        
        return False
    
    # Extract key metrics
    metrics = {
        "timestamp": validation_results.get("timestamp", datetime.datetime.now().isoformat()),
        "total_tools": validation_results.get("total_tools", 0),
        "tools_with_doi": validation_results.get("tools_with_doi", 0),
        "doi_percentage": validation_results.get("doi_percentage", 0),
        "valid_doi_percentage": validation_results.get("valid_doi_percentage", 0),
        "consistent_doi_percentage": validation_results.get("consistent_doi_percentage", 0),
        "tools_with_pubmed_data": validation_results.get("tools_with_pubmed_data", 0),
        "pubmed_percentage": validation_results.get("pubmed_percentage", 0),
        "tools_with_academic_impact": validation_results.get("tools_with_academic_impact", 0),
        "academic_impact_percentage": validation_results.get("academic_impact_percentage", 0),
        "citation_format_counts": validation_results.get("citation_format_counts", {"bibtex": 0, "apa": 0, "mla": 0}),
        "issues_count": len(validation_results.get("all_issues", [])),
        "total_dois_fixed": validation_results.get("total_dois_fixed", 0)
    }
    
    with open(METRICS_FILE_PATH, 'w') as f:
        json.dump(metrics, f, indent=2)
    
    return True

def commit_changes(fixed_dois: int) -> Tuple[bool, str]:
    """Commit changes to a new branch and return success and branch name"""
    # Create a unique branch name with timestamp
    timestamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    branch_name = f"citation-validation-{timestamp}"
    
    print(f"Creating branch {branch_name}...")
    status, stdout, stderr = run_cmd(f"git checkout -b {branch_name}")
    if status != 0:
        print(f"Error creating branch: {stderr}")
        return False, ""
    
    # Add changes
    print("Adding metadata files...")
    run_cmd("git add metadata/**/*.json || true")
    
    print("Adding metrics file...")
    run_cmd("git add reports/citations/citation_validation_metrics.json || true")
    
    # Check if there are changes to commit
    status, stdout, stderr = run_cmd("git diff --staged --quiet")
    if status == 0:
        print("No changes to commit.")
        return False, ""
    
    # Create commit message
    if fixed_dois > 0:
        commit_message = f"Update citation validation metrics and fix {fixed_dois} DOI format issues"
    else:
        commit_message = "Update citation validation metrics"
    
    print(f"Committing changes: {commit_message}")
    status, stdout, stderr = run_cmd(f'git commit -m "{commit_message}"')
    if status != 0:
        print(f"Error committing changes: {stderr}")
        return False, ""
    
    # Push branch
    print(f"Pushing branch {branch_name}...")
    status, stdout, stderr = run_cmd(f"git push origin {branch_name}")
    if status != 0:
        print(f"Error pushing branch: {stderr}")
        return False, ""
    
    return True, branch_name

def create_pull_request(branch_name: str, fixed_dois: int) -> bool:
    """Create a pull request for the changes"""
    # Get repository info from environment
    github_repository = os.environ.get("GITHUB_REPOSITORY", "")
    if not github_repository:
        print("GITHUB_REPOSITORY environment variable not set!")
        return False
    
    # Create PR title and body
    if fixed_dois > 0:
        title = f"Update citation validation metrics and fix {fixed_dois} DOI format issues"
        body = f"""## Automated PR with citation validation results

This PR includes:
- Updated citation validation metrics in `reports/citations/citation_validation_metrics.json`
- Automatically fixed {fixed_dois} DOI formatting issues in metadata files

Generated automatically by the Citation Data Validation workflow."""
    else:
        title = "Update citation validation metrics"
        body = """## Automated PR with citation validation results

This PR includes:
- Updated citation validation metrics in `reports/citations/citation_validation_metrics.json`

Generated automatically by the Citation Data Validation workflow."""
    
    # Create PR using GitHub API
    github_token = os.environ.get("GITHUB_TOKEN")
    if not github_token:
        print("GITHUB_TOKEN environment variable not set!")
        return False
    
    try:
        g = Github(github_token)
        repo = g.get_repo(github_repository)
        pr = repo.create_pull(
            title=title,
            body=body,
            head=branch_name,
            base="main"
        )
        print(f"Created PR #{pr.number}: {pr.html_url}")
        return True
    except Exception as e:
        print(f"Error creating PR: {e}")
        
        # Fallback to gh CLI
        print("Trying with gh CLI...")
        status, stdout, stderr = run_cmd(f'gh pr create --base main --head {branch_name} --title "{title}" --body "{body}"')
        if status != 0:
            print(f"Error creating PR with gh CLI: {stderr}")
            return False
        
        print(f"Created PR with gh CLI: {stdout}")
        return True

def main() -> int:
    """Main function to process validation results"""
    # Set up git config
    setup_git_config()
    
    # Load validation results
    validation_results = load_validation_results()
    fixed_dois = validation_results.get("total_dois_fixed", 0)
    print(f"Found {fixed_dois} fixed DOIs in validation results.")
    
    # Update metrics file
    if not update_metrics_file(validation_results):
        print("Failed to update metrics file.")
        return 1
    
    # Commit changes
    changes_committed, branch_name = commit_changes(fixed_dois)
    if not changes_committed:
        print("No changes committed.")
        return 0
    
    # Create PR
    pr_created = create_pull_request(branch_name, fixed_dois)
    if not pr_created:
        print("Failed to create PR.")
        return 1
    
    print("Successfully processed validation results and created PR.")
    return 0

if __name__ == "__main__":
    sys.exit(main())