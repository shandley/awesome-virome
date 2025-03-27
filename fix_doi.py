#!/usr/bin/env python3
"""
Script to fix DOIs with trailing parentheses in the data.json file
"""

import json
import sys
from pathlib import Path

def fix_dois(data_file="data.json"):
    """Fix DOIs with trailing parentheses in data.json"""
    # Load the data
    with open(data_file, 'r') as f:
        data = json.load(f)
    
    # Track changes
    changes_made = 0
    
    # Fix DOIs
    for node in data.get("nodes", []):
        if "doi" in node and node["doi"] and node["doi"].endswith(")"):
            old_doi = node["doi"]
            # Remove trailing parenthesis
            node["doi"] = node["doi"].rstrip(")")
            print(f"Fixed DOI for {node['name']}: {old_doi} -> {node['doi']}")
            changes_made += 1
    
    # Save the updated data
    if changes_made > 0:
        with open(data_file, 'w') as f:
            json.dump(data, f, indent=2)
        print(f"Fixed {changes_made} DOIs in {data_file}")
    else:
        print("No DOIs needed fixing")

if __name__ == "__main__":
    data_file = "data.json"
    if len(sys.argv) > 1:
        data_file = sys.argv[1]
    
    fix_dois(data_file)