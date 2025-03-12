#!/usr/bin/env python3
"""
Script to validate example workflows defined in the README.md file.
This checks that workflow tool combinations make sense and are consistent.
"""

import os
import re
import sys
import json
import logging
from pathlib import Path
from collections import defaultdict

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

def load_readme(readme_path):
    """Load the README.md file."""
    try:
        with open(readme_path, 'r') as f:
            readme_content = f.read()
        return readme_content
    except FileNotFoundError as e:
        logger.error(f"Error loading README.md: {e}")
        sys.exit(1)

def load_data_json(data_json_path):
    """Load and parse the data.json file."""
    try:
        with open(data_json_path, 'r') as f:
            data = json.load(f)
        return data
    except (FileNotFoundError, json.JSONDecodeError) as e:
        logger.error(f"Error loading data.json: {e}")
        sys.exit(1)

def extract_workflows_from_readme(readme_content):
    """Extract workflow definitions from README.md."""
    # Look for workflow sections
    workflow_sections = re.findall(r'### ([^#\n]+)Workflow:?\n([\s\S]+?)(?=\n##|\Z)', readme_content)
    
    workflows = []
    for section_name, section_content in workflow_sections:
        # Extract workflow steps
        steps = []
        for line in section_content.splitlines():
            # Look for numbered steps
            step_match = re.match(r'\s*\d+\.\s*(.*)', line)
            if step_match:
                step_text = step_match.group(1).strip()
                steps.append(step_text)
        
        # If we found steps, record this workflow
        if steps:
            workflows.append({
                'name': section_name.strip(),
                'steps': steps
            })
    
    return workflows

def extract_tools_from_data_json(data):
    """Extract tools and categories from data.json."""
    tools_by_name = {}
    tools_by_category = defaultdict(list)
    
    for node in data.get('nodes', []):
        if node.get('type') == 'tool':
            name = node.get('name')
            category = node.get('category')
            subcategory = node.get('subcategory')
            url = node.get('url')
            
            tools_by_name[name] = {
                'category': category,
                'subcategory': subcategory,
                'url': url
            }
            
            if category:
                tools_by_category[category].append(name)
    
    return tools_by_name, tools_by_category

def extract_tool_references_from_step(step_text):
    """Extract tool references from a workflow step."""
    # Extract text within parentheses as these often contain tool lists
    paren_matches = re.findall(r'\(e\.g\.,\s*([^)]+)\)', step_text)
    
    tools = []
    for match in paren_matches:
        # Split by commas or 'or' to get individual tool names
        parts = re.split(r',\s*|\s+or\s+', match)
        tools.extend([part.strip() for part in parts if part.strip()])
    
    # Also look for bold tool names outside parentheses
    bold_matches = re.findall(r'\*\*([^*]+)\*\*', step_text)
    tools.extend([match.strip() for match in bold_matches])
    
    return tools

def validate_workflow(workflow, tools_by_name, tools_by_category):
    """Validate a workflow definition."""
    issues = []
    referenced_tools = []
    
    for i, step in enumerate(workflow['steps']):
        # Extract tools referenced in this step
        step_tools = extract_tool_references_from_step(step)
        
        # Check if referenced tools exist in data.json
        for tool in step_tools:
            if tool not in tools_by_name:
                # Check if it might be a category instead
                found_in_category = False
                for category, category_tools in tools_by_category.items():
                    if any(tool.lower() in t.lower() for t in category_tools):
                        found_in_category = True
                        break
                
                if not found_in_category:
                    issues.append(f"Step {i+1}: Tool '{tool}' referenced but not found in data.json")
            else:
                referenced_tools.append(tool)
    
    # Check for duplicate steps
    step_texts = [step.lower() for step in workflow['steps']]
    for i, step in enumerate(step_texts):
        if step_texts.count(step) > 1:
            issues.append(f"Duplicate step found: '{workflow['steps'][i]}'")
    
    # Check for logical ordering issues
    # This is a simplified example - in a real implementation, you might have
    # more sophisticated logic to check for specific ordering requirements
    quality_control_keywords = ['quality', 'qc', 'filter']
    assembly_keywords = ['assembly', 'assembl', 'spades', 'megahit']
    annotation_keywords = ['annotat', 'function', 'DRAM', 'Pharokka']
    
    quality_control_steps = []
    assembly_steps = []
    annotation_steps = []
    
    for i, step in enumerate(workflow['steps']):
        step_lower = step.lower()
        if any(kw in step_lower for kw in quality_control_keywords):
            quality_control_steps.append(i)
        if any(kw in step_lower for kw in assembly_keywords):
            assembly_steps.append(i)
        if any(kw in step_lower for kw in annotation_keywords):
            annotation_steps.append(i)
    
    # Check if quality control comes before assembly
    if quality_control_steps and assembly_steps and min(assembly_steps) < max(quality_control_steps):
        issues.append("Workflow ordering issue: Assembly step appears before quality control step")
    
    # Check if assembly comes before annotation
    if assembly_steps and annotation_steps and min(annotation_steps) < max(assembly_steps):
        issues.append("Workflow ordering issue: Annotation step appears before assembly step")
    
    return issues, referenced_tools

def validate_workflows(workflows, tools_by_name, tools_by_category):
    """Validate all workflow definitions."""
    all_issues = {}
    all_referenced_tools = set()
    
    for workflow in workflows:
        issues, referenced_tools = validate_workflow(workflow, tools_by_name, tools_by_category)
        if issues:
            all_issues[workflow['name']] = issues
        all_referenced_tools.update(referenced_tools)
    
    return all_issues, all_referenced_tools

def main():
    """Main function to validate workflow definitions."""
    repo_root = Path(__file__).parent.parent
    readme_path = repo_root / 'README.md'
    data_json_path = repo_root / 'data.json'
    
    # Load files
    readme_content = load_readme(readme_path)
    data = load_data_json(data_json_path)
    
    # Extract workflows from README
    workflows = extract_workflows_from_readme(readme_content)
    logger.info(f"Found {len(workflows)} workflow definitions in README.md")
    
    # Extract tools from data.json
    tools_by_name, tools_by_category = extract_tools_from_data_json(data)
    logger.info(f"Found {len(tools_by_name)} tools in data.json")
    
    # Validate workflows
    workflow_issues, referenced_tools = validate_workflows(workflows, tools_by_name, tools_by_category)
    
    # Report results
    if workflow_issues:
        logger.error(f"Found issues in {len(workflow_issues)} workflows:")
        for workflow_name, issues in workflow_issues.items():
            logger.error(f"\nWorkflow '{workflow_name}':")
            for issue in issues:
                logger.error(f"  - {issue}")
        sys.exit(1)
    else:
        logger.info(f"All {len(workflows)} workflows are valid")
        logger.info(f"Referenced {len(referenced_tools)} unique tools in workflows")
        sys.exit(0)

if __name__ == "__main__":
    main()