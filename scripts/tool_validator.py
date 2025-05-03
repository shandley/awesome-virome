#!/usr/bin/env python3
"""
Tool Validator for Awesome-Virome Repository

This script provides enhanced validation for user-contributed tools,
ensuring they meet data quality standards before being added to the collection.
"""

import os
import re
import sys
import json
import logging
import requests
import argparse
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple, Set
import jsonschema
from urllib.parse import urlparse

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Constants
GITHUB_API = "https://api.github.com/repos/"
GITLAB_API = "https://gitlab.com/api/v4/projects/"
GITHUB_TOKEN = os.environ.get("GITHUB_TOKEN", "")
DOI_PATTERN = r"^10\.\d{4,9}/[-._;()/:A-Z0-9]+$"
CROSSREF_API = "https://api.crossref.org/works/"
QUALITY_THRESHOLD = 70  # Minimum quality score (out of 100) for tool entries
ROOT_DIR = Path(__file__).parent.parent

class ToolValidator:
    """Validates tool entries for the Awesome-Virome collection."""
    
    def __init__(self, schema_path: Optional[str] = None):
        """Initialize the validator with the tool schema."""
        if schema_path is None:
            schema_path = ROOT_DIR / "scripts" / "schema.json"
        
        try:
            with open(schema_path, 'r') as f:
                self.schema = json.load(f)
        except (FileNotFoundError, json.JSONDecodeError) as e:
            logger.error(f"Error loading schema: {e}")
            self.schema = None
            
        # Initialize result containers
        self.errors = []
        self.warnings = []
        self.suggestions = []
        self.quality_score = 0
    
    def validate_tool(self, tool_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Validate a tool entry against all quality checks.
        
        Args:
            tool_data: Dictionary containing tool data
            
        Returns:
            Dictionary with validation results
        """
        # Reset validation state
        self.errors = []
        self.warnings = []
        self.suggestions = []
        self.quality_score = 0
        
        # Basic schema validation
        schema_valid = self.validate_schema(tool_data)
        
        # URL checks
        url_valid, url_suggestions = self.validate_url(tool_data.get('url', ''))
        self.suggestions.extend(url_suggestions)
        
        # Description quality
        description_valid, desc_suggestions = self.validate_description(tool_data.get('description', ''))
        self.suggestions.extend(desc_suggestions)
        
        # Citation validation
        citation_valid, citation_suggestions = self.validate_citation(
            tool_data.get('doi', ''), 
            tool_data.get('citation', '')
        )
        self.suggestions.extend(citation_suggestions)
        
        # Category relevance
        category_valid, category_suggestions = self.validate_category(
            tool_data.get('category', ''),
            tool_data.get('description', '')
        )
        self.suggestions.extend(category_suggestions)
        
        # Calculate overall quality score
        self.calculate_quality_score(tool_data)
        
        # Compile validation results
        validation_result = {
            'valid': len(self.errors) == 0,
            'quality_score': self.quality_score,
            'meets_threshold': self.quality_score >= QUALITY_THRESHOLD,
            'errors': self.errors,
            'warnings': self.warnings,
            'suggestions': self.suggestions
        }
        
        return validation_result
    
    def validate_schema(self, tool_data: Dict[str, Any]) -> bool:
        """
        Validate tool data against JSON schema.
        
        Args:
            tool_data: Dictionary containing tool data
            
        Returns:
            True if valid, False otherwise
        """
        if self.schema is None:
            self.warnings.append("Schema validation skipped - schema not available")
            return True
        
        try:
            jsonschema.validate(instance=tool_data, schema=self.schema)
            return True
        except jsonschema.exceptions.ValidationError as e:
            self.errors.append(f"Schema validation failed: {e.message}")
            return False
    
    def validate_url(self, url: str) -> Tuple[bool, List[str]]:
        """
        Validate tool URL for accessibility and correctness.
        
        Args:
            url: Tool URL to validate
            
        Returns:
            Tuple of (valid, suggestions)
        """
        suggestions = []
        
        if not url:
            self.errors.append("URL is required")
            return False, suggestions
        
        # Check URL format
        try:
            parsed_url = urlparse(url)
            if not all([parsed_url.scheme, parsed_url.netloc]):
                self.errors.append(f"Invalid URL format: {url}")
                return False, suggestions
        except Exception:
            self.errors.append(f"Invalid URL format: {url}")
            return False, suggestions
        
        # Check if URL is accessible
        try:
            response = requests.head(
                url, 
                timeout=10,
                allow_redirects=True,
                headers={'User-Agent': 'Mozilla/5.0 Awesome-Virome-Validator/1.0'}
            )
            
            if response.status_code >= 400:
                self.errors.append(f"URL not accessible (HTTP {response.status_code}): {url}")
                return False, suggestions
            
            # Check for redirects
            if len(response.history) > 0:
                suggestions.append(f"URL redirects to {response.url} - consider updating to the final URL")
                
        except requests.RequestException as e:
            self.errors.append(f"Error checking URL {url}: {str(e)}")
            return False, suggestions
        
        # Repository-specific validations
        if 'github.com' in url:
            github_valid, github_suggestions = self.validate_github_repo(url)
            suggestions.extend(github_suggestions)
            if not github_valid:
                return False, suggestions
        
        return True, suggestions
    
    def validate_github_repo(self, url: str) -> Tuple[bool, List[str]]:
        """
        Validate GitHub repository characteristics.
        
        Args:
            url: GitHub repository URL
            
        Returns:
            Tuple of (valid, suggestions)
        """
        suggestions = []
        
        # Extract owner/repo from URL
        match = re.search(r'github\.com/([^/]+/[^/]+)', url)
        if not match:
            return True, suggestions  # Not a standard GitHub URL, skip additional checks
        
        repo_path = match.group(1).rstrip('/')
        
        # Remove any extra path components
        if "/" in repo_path.split("/", 2)[1]:
            repo_path = "/".join(repo_path.split("/", 2)[:2])
        
        # Create API URL
        api_url = f"{GITHUB_API}{repo_path}"
        
        headers = {"Accept": "application/vnd.github.v3+json"}
        if GITHUB_TOKEN:
            headers["Authorization"] = f"token {GITHUB_TOKEN}"
        
        try:
            response = requests.get(api_url, headers=headers, timeout=10)
            
            if response.status_code == 404:
                self.errors.append(f"GitHub repository not found: {url}")
                return False, suggestions
                
            response.raise_for_status()
            
            data = response.json()
            
            # Check if repository is archived
            if data.get('archived', False):
                self.warnings.append(f"Repository is archived: {url}")
                suggestions.append("Consider if this archived repository should still be included")
            
            # Check for README
            readme_url = f"{GITHUB_API}{repo_path}/readme"
            readme_response = requests.get(readme_url, headers=headers, timeout=10)
            
            if readme_response.status_code == 404:
                self.warnings.append(f"Repository has no README: {url}")
                suggestions.append("A repository without a README may be hard for users to understand")
            
            # Check for recent activity
            updated_at = data.get('updated_at') or data.get('pushed_at')
            if updated_at:
                from datetime import datetime, timezone
                last_updated = datetime.strptime(updated_at, '%Y-%m-%dT%H:%M:%SZ')
                now = datetime.now(timezone.utc)
                days_since_update = (now - last_updated.replace(tzinfo=timezone.utc)).days
                
                if days_since_update > 365:
                    self.warnings.append(f"Repository has not been updated in {days_since_update} days")
                    suggestions.append("This tool may not be actively maintained")
                elif days_since_update > 180:
                    self.warnings.append(f"Repository has not been updated in {days_since_update} days")
            
            # Check for issues/PRs
            if data.get('open_issues_count', 0) > 0:
                open_issues = data.get('open_issues_count', 0)
                # Calculate approximate PR count (GitHub counts PRs as issues in this API)
                issues_url = f"{GITHUB_API}{repo_path}/issues?state=open&per_page=1"
                pr_url = f"{GITHUB_API}{repo_path}/pulls?state=open&per_page=1"
                
                try:
                    issues_response = requests.get(issues_url, headers=headers, timeout=10)
                    pr_response = requests.get(pr_url, headers=headers, timeout=10)
                    
                    if issues_response.ok and pr_response.ok:
                        issue_count = int(issues_response.headers.get('Link', '').split('page=')[-1].split('>', 1)[0]) if 'Link' in issues_response.headers else 0
                        pr_count = int(pr_response.headers.get('Link', '').split('page=')[-1].split('>', 1)[0]) if 'Link' in pr_response.headers else 0
                        
                        if issue_count > 20 and pr_count < issue_count * 0.2:
                            self.warnings.append(f"Repository has a high number of open issues ({issue_count})")
                            suggestions.append("High number of open issues might indicate maintenance concerns")
                except Exception:
                    pass  # Skip detailed issue analysis if it fails
                    
        except Exception as e:
            self.warnings.append(f"Error checking GitHub repository details: {str(e)}")
        
        return True, suggestions
    
    def validate_description(self, description: str) -> Tuple[bool, List[str]]:
        """
        Validate tool description quality.
        
        Args:
            description: Tool description
            
        Returns:
            Tuple of (valid, suggestions)
        """
        suggestions = []
        
        if not description:
            self.errors.append("Description is required")
            return False, suggestions
        
        # Check description length
        if len(description) < 10:
            self.errors.append("Description is too short (minimum 10 characters)")
            return False, suggestions
        
        if len(description) > 500:
            self.warnings.append("Description is very long (over 500 characters)")
            suggestions.append("Consider shortening the description for better readability")
        
        # Check for complete sentences
        if not description[0].isupper():
            suggestions.append("Description should start with a capital letter")
        
        # Check for missing periods at the end of sentences
        if not description.endswith('.'):
            suggestions.append("Description should end with a period")
        
        # Check for common keywords expected in tool descriptions
        virome_keywords = ['virus', 'viral', 'phage', 'virome', 'metagenome', 'sequenc']
        has_virome_keywords = any(keyword in description.lower() for keyword in virome_keywords)
        
        if not has_virome_keywords:
            self.warnings.append("Description doesn't mention common virome analysis terms")
            suggestions.append("Consider adding terms like 'virus', 'phage', or 'virome' to clarify relevance")
        
        # Check for placeholder text
        placeholder_patterns = [
            r'Lorem ipsum', 
            r'TODO', 
            r'to be completed', 
            r'add description'
        ]
        if any(re.search(pattern, description, re.IGNORECASE) for pattern in placeholder_patterns):
            self.errors.append("Description contains placeholder text")
            return False, suggestions
        
        return True, suggestions
    
    def validate_citation(self, doi: str, citation: str) -> Tuple[bool, List[str]]:
        """
        Validate DOI and citation information.
        
        Args:
            doi: Digital Object Identifier
            citation: Citation text
            
        Returns:
            Tuple of (valid, suggestions)
        """
        suggestions = []
        
        # If neither DOI nor citation is provided, that's okay but suggest adding them
        if not doi and not citation:
            suggestions.append("Adding a DOI or citation can help users cite your tool")
            return True, suggestions
        
        # Validate DOI format if provided
        if doi:
            if not re.match(DOI_PATTERN, doi, re.IGNORECASE):
                self.errors.append(f"Invalid DOI format: {doi}")
                suggestions.append("DOI should follow the format 10.NNNN/XXX")
                return False, suggestions
            
            # Check if DOI is resolvable (if rate limiting allows)
            try:
                headers = {
                    'Accept': 'application/json'
                }
                response = requests.get(f"{CROSSREF_API}{doi}", headers=headers, timeout=10)
                
                if response.status_code != 200:
                    self.warnings.append(f"DOI cannot be resolved: {doi}")
                    suggestions.append("Check that the DOI is correct and currently active")
                else:
                    # Check if we have citation but it doesn't match the DOI
                    if citation:
                        # Simple check - see if the DOI appears in the citation
                        if doi not in citation:
                            suggestions.append("Citation text should include the DOI for consistency")
            except requests.RequestException:
                # Don't fail validation if we can't check the DOI (could be rate limiting)
                self.warnings.append(f"Could not verify DOI: {doi}")
        
        # Validate citation format if provided
        if citation:
            # Check for common citation formats
            citation_patterns = {
                'apa': r'.+\(\d{4}\)\..+\..+\.',  # Author(s) (Year). Title. Source.
                'mla': r'.+\..+\..+\d{4}\.',  # Author(s). Title. Source, Year.
                'bibtex': r'@\w+\s*{.+,.+}'  # @article{key, field = {value}, ...}
            }
            
            matches_format = any(
                re.match(pattern, citation, re.DOTALL)
                for pattern in citation_patterns.values()
            )
            
            if not matches_format and not doi:
                self.warnings.append("Citation doesn't match common formats (APA, MLA, BibTeX)")
                suggestions.append("Consider using a standard citation format")
        
        return True, suggestions
    
    def validate_category(self, category: str, description: str) -> Tuple[bool, List[str]]:
        """
        Validate category relevance to description.
        
        Args:
            category: Tool category
            description: Tool description
            
        Returns:
            Tuple of (valid, suggestions)
        """
        suggestions = []
        
        if not category:
            self.errors.append("Category is required")
            return False, suggestions
        
        # Simple keyword matching for category relevance
        category_keywords = {
            "Virus and Phage Identification": [
                'identif', 'detect', 'recogni', 'phage', 'virus', 'viral', 'discovery'
            ],
            "Metagenome Analysis": [
                'metagenom', 'microbiome', 'communit', 'environment', 'ecolog'
            ],
            "Host Prediction": [
                'host', 'predict', 'interaction', 'infect', 'target'
            ],
            "Genome Analysis": [
                'genom', 'analys', 'sequenc', 'feature'
            ],
            "Genome Annotation": [
                'annotat', 'feature', 'gene', 'predict', 'function'
            ],
            "Genome Assembly": [
                'assembl', 'contig', 'scaffold', 'reconstruct'
            ],
            "Taxonomy": [
                'taxonom', 'classif', 'phyl', 'evolution', 'clade'
            ],
            "Visualization and Infrastructure": [
                'visuali', 'plot', 'graph', 'dashboard', 'interface', 'infrastructure'
            ],
            "Machine Learning Models": [
                'machine learning', 'ml', 'ai', 'model', 'neural', 'predict', 'deep learning'
            ]
        }
        
        # Check if category keywords appear in description
        if category in category_keywords and description:
            category_match = any(
                keyword in description.lower()
                for keyword in category_keywords[category]
            )
            
            if not category_match:
                self.warnings.append(f"Description doesn't contain keywords related to category '{category}'")
                suggestions.append(f"Consider adding terms related to {category} in the description")
        
        # Check for potentially better category matches
        better_matches = []
        if description:
            for cat, keywords in category_keywords.items():
                if cat != category:  # Skip current category
                    match_score = sum(
                        1 for keyword in keywords
                        if keyword in description.lower()
                    )
                    if match_score >= 2:  # At least 2 keywords match
                        better_matches.append(cat)
        
        if better_matches:
            suggestions.append(f"Based on description, also consider these categories: {', '.join(better_matches)}")
        
        return True, suggestions
    
    def calculate_quality_score(self, tool_data: Dict[str, Any]) -> None:
        """
        Calculate an overall quality score for the tool entry.
        
        Args:
            tool_data: Tool data dictionary
        """
        # Start with base score of 100
        score = 100
        
        # Deduct for each error (major issues)
        score -= len(self.errors) * 20
        
        # Deduct for warnings (minor issues)
        score -= len(self.warnings) * 5
        
        # Deduct for missing recommended fields
        recommended_fields = ['license', 'language', 'installation', 'citation']
        for field in recommended_fields:
            if not tool_data.get(field):
                score -= 5
                self.suggestions.append(f"Adding '{field}' information would improve the quality score")
        
        # Bonus for comprehensive descriptions
        if tool_data.get('description') and len(tool_data.get('description', '')) > 100:
            score += 5
        
        # Bonus for having a DOI
        if tool_data.get('doi'):
            score += 10
        
        # Bonus for having both citation and DOI
        if tool_data.get('doi') and tool_data.get('citation'):
            score += 5
        
        # Bonus for having installation methods
        if tool_data.get('installation'):
            score += 5
        
        # Ensure score is between 0 and 100
        self.quality_score = max(0, min(100, score))


# Command-line interface
def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Validate tool entries for the Awesome-Virome collection')
    parser.add_argument('--input', '-i', help='Path to JSON file containing tool data')
    parser.add_argument('--schema', '-s', help='Path to JSON schema file')
    parser.add_argument('--output', '-o', help='Path to write validation results')
    parser.add_argument('--url', '-u', help='Tool URL to validate')
    parser.add_argument('--verbose', '-v', action='store_true', help='Enable verbose output')
    
    return parser.parse_args()


def main():
    """Main function."""
    args = parse_arguments()
    
    # Configure logging level
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    
    # Initialize validator
    validator = ToolValidator(args.schema)
    
    # Process input based on arguments
    if args.input:
        try:
            with open(args.input, 'r') as f:
                tool_data = json.load(f)
        except (FileNotFoundError, json.JSONDecodeError) as e:
            logger.error(f"Error loading input file: {e}")
            sys.exit(1)
            
        # Validate the tool
        validation_result = validator.validate_tool(tool_data)
        
        # Print results
        if validation_result['valid']:
            logger.info(f"Validation successful with quality score: {validation_result['quality_score']}/100")
        else:
            logger.error(f"Validation failed with quality score: {validation_result['quality_score']}/100")
            for error in validation_result['errors']:
                logger.error(f"Error: {error}")
        
        for warning in validation_result['warnings']:
            logger.warning(f"Warning: {warning}")
            
        for suggestion in validation_result['suggestions']:
            logger.info(f"Suggestion: {suggestion}")
        
        # Write results to output file if specified
        if args.output:
            try:
                with open(args.output, 'w') as f:
                    json.dump(validation_result, f, indent=2)
                logger.info(f"Validation results written to {args.output}")
            except IOError as e:
                logger.error(f"Error writing output file: {e}")
    
    # Quick URL validation
    elif args.url:
        # Create minimal tool data with just the URL
        tool_data = {
            "name": "URL Validation Check",
            "url": args.url,
            "description": "Temporary description for URL validation check",
            "category": "Other"
        }
        
        # Validate just the URL
        _, url_suggestions = validator.validate_url(args.url)
        
        # Print results
        if validator.errors:
            logger.error("URL validation failed")
            for error in validator.errors:
                logger.error(f"Error: {error}")
        else:
            logger.info("URL is valid")
            
        for suggestion in url_suggestions:
            logger.info(f"Suggestion: {suggestion}")
    
    else:
        logger.error("No input provided. Use --input, --url, or --help for usage information.")
        sys.exit(1)


if __name__ == "__main__":
    main()