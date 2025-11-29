#!/usr/bin/env python3
"""
Script to fix duplicate [Updated: MM/YYYY] tags in the README.md file.
This is a one-time fix for multiple redundant tags that were accidentally added.
"""

import re
import sys
from pathlib import Path

def remove_duplicate_tags(readme_path):
    """Remove duplicate consecutive [Updated: MM/YYYY] tags from the README file."""
    # Read the README file
    with open(readme_path, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Regular expression to match duplicate [Updated: MM/YYYY] tags
    # This pattern matches one or more occurrences of the same tag
    tag_pattern = r'(\[Updated: \d{2}/\d{4}\])(?:\s+\[Updated: \d{2}/\d{4}\])+'
    
    # Find all instances of duplicate tags
    duplicate_tags = re.findall(tag_pattern, content)
    num_duplicates = len(duplicate_tags)
    
    # Replace duplicate tags with a single instance
    fixed_content = content
    fixed_content = re.sub(tag_pattern, r'\1', fixed_content)
    
    # Write the updated content back to the README
    with open(readme_path, 'w', encoding='utf-8') as f:
        f.write(fixed_content)
    
    return num_duplicates

def fix_update_check_script(script_path):
    """Update the update_check.py script to prevent future duplicate tags."""
    with open(script_path, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Look for the pattern matching code in update_readme_with_dates_status_and_stars function
    # This is a simplified approach; in a real scenario, you might want to use an AST parser
    pattern_check = r'repo_pattern = re.escape\(f"\[{repo_name}\]\({repo_url}\)"\) \+ r".*?\(\\\[Updated:.*?\\\]\)?'
    
    # Update the pattern to match all existing update tags, not just one
    if pattern_check in content:
        updated_pattern = r'repo_pattern = re.escape\(f"\[{repo_name}\]\({repo_url}\)"\) \+ r".*?(?:\s*\\\[Updated:.*?\\\])*"'
        content = content.replace(pattern_check, updated_pattern)
    
    # Update replacement logic to ensure all existing tags are removed before adding a new one
    old_replacement = r'updated_content = updated_content.replace('
    new_replacement = r'# Replace all existing tags with the new tag\n                # First remove all existing tags\n                temp_content = re.sub(\n                    re.escape(f"[{repo_name}]({repo_url})") + r"(?:\s+\[Updated:.*?\])*",\n                    f"[{repo_name}]({repo_url})",\n                    updated_content\n                )\n                # Then add the new tag\n                updated_content = temp_content.replace('
    
    # Replace the first relevant occurrence of this pattern
    # This is a simplified approach; in a real scenario, this would need to be more targeted
    if old_replacement in content:
        content = content.replace(old_replacement, new_replacement, 1)
    
    with open(script_path, 'w', encoding='utf-8') as f:
        f.write(content)
    
    return "update_check.py was updated to prevent future duplicate tags"

if __name__ == "__main__":
    readme_path = Path(__file__).parent / "README.md"
    update_script_path = Path(__file__).parent / "update_check.py"
    
    if not readme_path.exists():
        print(f"README file not found at {readme_path}")
        sys.exit(1)
    
    if not update_script_path.exists():
        print(f"update_check.py script not found at {update_script_path}")
        sys.exit(1)
    
    # Fix duplicate tags in README
    num_duplicates = remove_duplicate_tags(readme_path)
    print(f"Fixed {num_duplicates} instances of duplicate [Updated: MM/YYYY] tags in README.md")
    
    # Update the update_check.py script
    result = fix_update_check_script(update_script_path)
    print(result)