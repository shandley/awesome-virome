#!/usr/bin/env python3
"""
Validate GitHub workflow files for syntax errors

This script checks all workflow files in the .github/workflows directory
for basic YAML syntax and common workflow syntax errors.
"""

import os
import sys
import yaml
import glob
from pathlib import Path
import re

def validate_workflow_file(file_path):
    """
    Validate a GitHub workflow file for syntax errors
    
    Args:
        file_path: Path to the workflow file
        
    Returns:
        List of error messages, empty if no errors found
    """
    errors = []
    
    try:
        # Check if file exists
        if not os.path.exists(file_path):
            return [f"File not found: {file_path}"]
        
        # Read the file
        with open(file_path, 'r') as f:
            content = f.read()
        
        # Check YAML syntax
        try:
            workflow = yaml.safe_load(content)
        except yaml.YAMLError as e:
            errors.append(f"YAML syntax error in {file_path}: {e}")
            return errors
        
        # Check for required sections
        if not workflow:
            errors.append(f"Empty workflow file: {file_path}")
            return errors
        
        # Check for 'name' field
        if 'name' not in workflow:
            errors.append(f"Missing 'name' field in {file_path}")
        
        # Check for 'on' field (triggers)
        if 'on' not in workflow:
            errors.append(f"Missing 'on' field in {file_path}")
        
        # Check for 'jobs' field
        if 'jobs' not in workflow:
            errors.append(f"Missing 'jobs' field in {file_path}")
        else:
            # Check that jobs is a dictionary
            if not isinstance(workflow['jobs'], dict):
                errors.append(f"'jobs' field must be a dictionary in {file_path}")
            else:
                # Check each job
                for job_id, job in workflow['jobs'].items():
                    # Check that job has 'runs-on'
                    if 'runs-on' not in job:
                        errors.append(f"Missing 'runs-on' field in job '{job_id}' in {file_path}")
                    
                    # Check that job has 'steps'
                    if 'steps' not in job:
                        errors.append(f"Missing 'steps' field in job '{job_id}' in {file_path}")
                    elif not isinstance(job['steps'], list):
                        errors.append(f"'steps' field must be a list in job '{job_id}' in {file_path}")
        
        # Check for common issues in expression syntax
        # This is a basic check and may not catch all issues
        expression_pattern = r'\$\{\{.*?\}\}'
        expressions = re.findall(expression_pattern, content)
        for expr in expressions:
            # Check for missing closing brackets
            if expr.count('{{') != expr.count('}}'):
                errors.append(f"Mismatched curly braces in expression: {expr}")
            
            # Check for unbalanced parentheses
            if expr.count('(') != expr.count(')'):
                errors.append(f"Unbalanced parentheses in expression: {expr}")
            
            # Check for common errors in if conditions
            if 'if: ' in expr:
                if '= ' in expr and not '== ' in expr:
                    errors.append(f"Possible assignment instead of comparison in condition: {expr}")
        
        return errors
    
    except Exception as e:
        return [f"Error validating {file_path}: {e}"]

def validate_all_workflows(workflows_dir):
    """
    Validate all workflow files in the given directory
    
    Args:
        workflows_dir: Directory containing workflow files
        
    Returns:
        Dictionary of file paths to error lists
    """
    validation_results = {}
    
    # Find all YAML files in the workflows directory
    workflow_files = glob.glob(os.path.join(workflows_dir, "*.yml"))
    
    print(f"Found {len(workflow_files)} workflow files")
    
    # Validate each file
    for file_path in workflow_files:
        print(f"Validating {file_path}...")
        errors = validate_workflow_file(file_path)
        validation_results[file_path] = errors
    
    return validation_results

def main():
    """Main entry point for the script"""
    # Find the .github/workflows directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.dirname(os.path.dirname(script_dir))
    workflows_dir = os.path.join(repo_root, ".github", "workflows")
    
    if not os.path.exists(workflows_dir):
        print(f"Workflows directory not found: {workflows_dir}")
        sys.exit(1)
    
    # Validate all workflows
    validation_results = validate_all_workflows(workflows_dir)
    
    # Report results
    print("\nWorkflow Validation Results:")
    print("============================")
    
    has_errors = False
    for file_path, errors in validation_results.items():
        file_name = os.path.basename(file_path)
        if errors:
            has_errors = True
            print(f"\n{file_name}: {len(errors)} error(s)")
            for error in errors:
                print(f"  - {error}")
        else:
            print(f"{file_name}: No errors")
    
    # Exit with appropriate status code
    if has_errors:
        print("\nValidation failed with errors.")
        sys.exit(1)
    else:
        print("\nAll workflows validated successfully!")
        sys.exit(0)

if __name__ == "__main__":
    main()