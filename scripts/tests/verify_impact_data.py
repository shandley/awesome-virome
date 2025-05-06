#!/usr/bin/env python3
"""
Verify the impact_data.json file for synthetic data

This script analyzes the impact_data.json file to detect any synthetic or weighted data.
It checks for patterns that would indicate synthetic data generation, such as unrealistic
citation distributions or estimates.
"""

import os
import sys
import json
import argparse
from pathlib import Path
from collections import Counter
import statistics
import math

def load_impact_data(file_path):
    """
    Load the impact_data.json file
    
    Args:
        file_path: Path to the impact_data.json file
        
    Returns:
        The loaded JSON data as a dictionary
    """
    try:
        with open(file_path, 'r') as f:
            return json.load(f)
    except (json.JSONDecodeError, IOError) as e:
        print(f"Error loading impact_data.json: {e}")
        sys.exit(1)

def detect_synthetic_data_patterns(impact_data):
    """
    Detect patterns that would indicate synthetic data generation
    
    Args:
        impact_data: The impact_data.json content as a dictionary
        
    Returns:
        List of suspicious patterns detected
    """
    suspicious_patterns = []
    
    # Check for tools data
    if "tools" not in impact_data or not isinstance(impact_data["tools"], list):
        suspicious_patterns.append("No tools list found in impact_data.json")
        return suspicious_patterns
    
    # Skip empty tools list
    if not impact_data["tools"]:
        suspicious_patterns.append("Empty tools list in impact_data.json")
        return suspicious_patterns
    
    # Analyze citation distributions
    yearly_pattern_checks = []
    influential_ratio_checks = []
    
    for tool in impact_data["tools"]:
        name = tool.get("name", "Unknown tool")
        
        # Skip tools without citation data
        if "total_citations" not in tool or not tool["total_citations"]:
            continue
        
        total_citations = tool["total_citations"]
        
        # Check for suspicious influential citations
        if "influential_citations" in tool and tool["influential_citations"]:
            influential = tool["influential_citations"]
            ratio = influential / total_citations
            
            # Check if ratio is suspiciously consistent
            if ratio == 0.1 or ratio == 0.2:
                suspicious_patterns.append(
                    f"Tool '{name}' has influential citations exactly {ratio*100}% of total citations, "
                    f"suggesting a synthetic estimate"
                )
            
            influential_ratio_checks.append(ratio)
        
        # Check for yearly citation patterns
        if "citations_by_year" in tool and tool["citations_by_year"]:
            yearly = tool["citations_by_year"]
            
            # Skip if no yearly data
            if not yearly:
                continue
            
            # Check if yearly citations sum to total
            yearly_sum = sum(yearly.values())
            if abs(yearly_sum - total_citations) > 1:  # Allow small rounding differences
                suspicious_patterns.append(
                    f"Tool '{name}' yearly citations sum ({yearly_sum}) differs significantly "
                    f"from total citations ({total_citations})"
                )
            
            # Check for suspiciously smooth distributions
            years = sorted(yearly.keys())
            if len(years) >= 3:
                values = [yearly[y] for y in years]
                
                # Calculate how linear/smooth the progression is
                diffs = [values[i+1] - values[i] for i in range(len(values)-1)]
                
                # Check if differences are very consistent
                if len(set(diffs)) == 1:
                    suspicious_patterns.append(
                        f"Tool '{name}' has perfectly linear citation growth "
                        f"({', '.join(f'{y}:{yearly[y]}' for y in years)}), "
                        f"suggesting synthetic data"
                    )
                
                # Check for perfect exponential growth (common in synthetic data)
                ratios = [values[i+1]/values[i] for i in range(len(values)-1) if values[i] > 0]
                if len(ratios) >= 2 and max(ratios) - min(ratios) < 0.1:
                    suspicious_patterns.append(
                        f"Tool '{name}' has suspiciously consistent growth rate "
                        f"({', '.join(f'{y}:{yearly[y]}' for y in years)}), "
                        f"suggesting synthetic data"
                    )
                
                yearly_pattern_checks.append(len(set(diffs)) / len(diffs))
    
    # Check for consistent patterns across the entire dataset
    if influential_ratio_checks:
        # If all influential citation ratios are the same, that's suspicious
        if len(set([round(r, 1) for r in influential_ratio_checks])) == 1:
            suspicious_patterns.append(
                f"All tools have the same influential citation ratio ({influential_ratio_checks[0]:.1f}), "
                f"suggesting synthetic data"
            )
    
    if yearly_pattern_checks:
        # If all yearly citation patterns are similarly smooth, that's suspicious
        avg_pattern = sum(yearly_pattern_checks) / len(yearly_pattern_checks)
        if avg_pattern < 0.3:  # Very smooth patterns
            suspicious_patterns.append(
                "Most tools have suspiciously smooth yearly citation patterns, "
                "suggesting synthetic data"
            )
    
    return suspicious_patterns

def verify_citation_sources(impact_data, metadata_dir):
    """
    Verify that citation data matches source data
    
    Args:
        impact_data: The impact_data.json content as a dictionary
        metadata_dir: Directory containing metadata files
        
    Returns:
        List of discrepancies found
    """
    discrepancies = []
    
    # Check for tools data
    if "tools" not in impact_data or not isinstance(impact_data["tools"], list):
        discrepancies.append("No tools list found in impact_data.json")
        return discrepancies
    
    # Skip empty tools list
    if not impact_data["tools"]:
        discrepancies.append("Empty tools list in impact_data.json")
        return discrepancies
    
    # Create a lookup of tools by name
    tools_by_name = {tool.get("name"): tool for tool in impact_data["tools"] if "name" in tool}
    
    # Check pubmed_citations directory if it exists
    pubmed_dir = os.path.join(metadata_dir, "pubmed_citations")
    if os.path.exists(pubmed_dir):
        for filename in os.listdir(pubmed_dir):
            if not filename.endswith(".json") or filename in ["summary.json", "pubmed_citations.json"]:
                continue
            
            file_path = os.path.join(pubmed_dir, filename)
            try:
                with open(file_path, 'r') as f:
                    citation_data = json.load(f)
                
                name = citation_data.get("name")
                if not name or name not in tools_by_name:
                    continue
                
                tool = tools_by_name[name]
                
                # Get citation count from the source
                citation_info = citation_data.get("citation_info", {})
                publication = citation_info.get("publication", {})
                source_citation_count = publication.get("citation_count", 0)
                
                # Compare with impact_data.json
                impact_citation_count = tool.get("total_citations", 0)
                
                # Check for significant differences
                if source_citation_count > 0 and impact_citation_count > 0:
                    if source_citation_count != impact_citation_count:
                        discrepancies.append(
                            f"Tool '{name}' has different citation counts in source ({source_citation_count}) "
                            f"and impact_data.json ({impact_citation_count})"
                        )
            
            except (json.JSONDecodeError, IOError) as e:
                discrepancies.append(f"Error reading {file_path}: {e}")
    
    return discrepancies

def get_category_consistency(impact_data):
    """
    Check for category consistency and possible heuristic categorization
    
    Args:
        impact_data: The impact_data.json content as a dictionary
        
    Returns:
        List of suspicious category patterns
    """
    suspicious_patterns = []
    
    # Check for tools data
    if "tools" not in impact_data or not isinstance(impact_data["tools"], list):
        suspicious_patterns.append("No tools list found in impact_data.json")
        return suspicious_patterns
    
    # Skip empty tools list
    if not impact_data["tools"]:
        suspicious_patterns.append("Empty tools list in impact_data.json")
        return suspicious_patterns
    
    # Count categories by name patterns
    categories_by_pattern = {}
    
    for tool in impact_data["tools"]:
        name = tool.get("name", "").lower()
        category = tool.get("category", "Other Tools")
        
        # Look for common patterns in names
        patterns = []
        if "phage" in name:
            patterns.append("phage")
        if "host" in name:
            patterns.append("host")
        if any(term in name for term in ["tax", "class", "sort"]):
            patterns.append("taxonomy")
        if any(term in name for term in ["find", "ident", "detect"]):
            patterns.append("identification")
        
        # Skip if no patterns found
        if not patterns:
            continue
        
        # Check if category matches patterns
        for pattern in patterns:
            if pattern not in categories_by_pattern:
                categories_by_pattern[pattern] = {"categories": Counter(), "total": 0}
            
            categories_by_pattern[pattern]["categories"][category] += 1
            categories_by_pattern[pattern]["total"] += 1
    
    # Check for suspicious consistency in categories
    for pattern, data in categories_by_pattern.items():
        most_common = data["categories"].most_common(1)
        if most_common:
            category, count = most_common[0]
            consistency = count / data["total"]
            
            # If >90% of tools with a pattern have the same category, that's suspicious
            if consistency > 0.9 and data["total"] >= 3:
                suspicious_patterns.append(
                    f"{count}/{data['total']} tools with '{pattern}' in the name "
                    f"are categorized as '{category}', suggesting heuristic categorization"
                )
    
    return suspicious_patterns

def main():
    """Main entry point for the script"""
    parser = argparse.ArgumentParser(description="Verify impact_data.json for synthetic data")
    parser.add_argument("--file", "-f", help="Path to impact_data.json (default: repo root)")
    args = parser.parse_args()
    
    # Find impact_data.json
    if args.file:
        impact_data_path = args.file
    else:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        repo_root = os.path.dirname(os.path.dirname(script_dir))
        impact_data_path = os.path.join(repo_root, "impact_data.json")
    
    metadata_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "metadata")
    
    # Check if file exists
    if not os.path.exists(impact_data_path):
        print(f"Error: impact_data.json not found at {impact_data_path}")
        sys.exit(1)
    
    print(f"Verifying {impact_data_path} for synthetic data...")
    
    # Load the file
    impact_data = load_impact_data(impact_data_path)
    
    # Check for synthetic data patterns
    suspicious_patterns = detect_synthetic_data_patterns(impact_data)
    
    # Check for discrepancies with source data
    discrepancies = verify_citation_sources(impact_data, metadata_dir)
    
    # Check for category consistency
    category_patterns = get_category_consistency(impact_data)
    
    # Report results
    print("\nVerification Results:")
    print("====================")
    
    if not suspicious_patterns and not discrepancies and not category_patterns:
        print("\nNo synthetic data detected in impact_data.json! ✅")
        print("All citation data appears to be legitimate and matches source data.")
        sys.exit(0)
    
    if suspicious_patterns:
        print("\nSuspicious citation patterns detected:")
        for pattern in suspicious_patterns:
            print(f"  - {pattern}")
    
    if discrepancies:
        print("\nDiscrepancies with source data found:")
        for discrepancy in discrepancies:
            print(f"  - {discrepancy}")
    
    if category_patterns:
        print("\nSuspicious category patterns detected:")
        for pattern in category_patterns:
            print(f"  - {pattern}")
    
    sys.exit(1)

if __name__ == "__main__":
    main()