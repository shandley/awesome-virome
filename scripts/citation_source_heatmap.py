#!/usr/bin/env python3
"""
Citation Source Heatmap Generator

This script generates a heatmap visualization showing which citation sources
were used for which tools. It helps to visualize the coverage and availability
of citation data across different sources.

Usage:
    python citation_source_heatmap.py [--output FILENAME]
"""

import os
import json
import argparse
import logging
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Any, Optional, Set, Tuple
import matplotlib.pyplot as plt
import seaborn as sns

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

# Constants
BASE_DIR = Path(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
IMPACT_DATA_PATH = BASE_DIR / "impact_data.json"
OUTPUT_DIR = BASE_DIR
KNOWN_SOURCES = ["scopus", "wos", "icite", "crossref"]
SOURCE_COLORS = {
    "scopus": "#5470c6",     # Blue
    "wos": "#91cc75",        # Green
    "icite": "#fac858",      # Yellow
    "crossref": "#ee6666",   # Red
    "unknown": "#cccccc"     # Gray
}


def load_impact_data() -> Dict[str, Any]:
    """Load the impact_data.json file."""
    try:
        with open(IMPACT_DATA_PATH, 'r') as f:
            data = json.load(f)
        return data
    except Exception as e:
        logger.error(f"Error loading impact data: {e}")
        return {}


def get_source_data(impact_data: Dict[str, Any]) -> Tuple[pd.DataFrame, Dict[str, List[str]]]:
    """
    Extract citation source data for each tool.
    
    Args:
        impact_data: The impact data dictionary
        
    Returns:
        DataFrame with tools and their citation sources,
        Categories dictionary mapping category names to tool lists
    """
    tools = impact_data.get('tools', [])
    
    # Prepare data for the heatmap
    data = []
    categories = {}
    
    for tool in tools:
        name = tool.get('name', 'Unknown')
        total_citations = tool.get('total_citations', 0)
        
        # Skip tools with zero citations
        if total_citations <= 0:
            continue
            
        # Get citation sources
        citation_source = tool.get('citation_source') or tool.get('primary_source')
        yearly_source = tool.get('yearly_citation_source') or tool.get('yearly_data_source')
        
        # Standardize source names
        if citation_source and citation_source.lower() in [s.lower() for s in KNOWN_SOURCES]:
            citation_source = next(s for s in KNOWN_SOURCES if s.lower() == citation_source.lower())
        else:
            citation_source = "unknown"
            
        if yearly_source and yearly_source.lower() in [s.lower() for s in KNOWN_SOURCES]:
            yearly_source = next(s for s in KNOWN_SOURCES if s.lower() == yearly_source.lower())
        else:
            yearly_source = "unknown"
        
        # Add to category dictionary
        category = tool.get('category', 'Other Tools')
        if category not in categories:
            categories[category] = []
        categories[category].append(name)
        
        # Add to data
        data.append({
            'Tool': name,
            'Total Citations': total_citations,
            'Citation Source': citation_source,
            'Yearly Citation Source': yearly_source,
            'Category': category
        })
    
    # Convert to DataFrame
    df = pd.DataFrame(data)
    return df, categories


def generate_heatmap(df: pd.DataFrame, output_path: str):
    """
    Generate a heatmap of tools and their citation sources.
    
    Args:
        df: DataFrame with citation source data
        output_path: Path to save the output image
    """
    if df.empty:
        logger.error("No data available for heatmap generation")
        return
    
    # Create a pivot table for the heatmap
    # Count the occurrences of each citation source for each category
    category_source_pivot = pd.crosstab(
        df['Category'], 
        df['Citation Source'],
        values=df['Total Citations'],
        aggfunc='sum'
    )
    
    # Ensure all known sources are in the columns
    for source in KNOWN_SOURCES:
        if source not in category_source_pivot.columns:
            category_source_pivot[source] = 0
    
    # Fill NaN with 0
    category_source_pivot = category_source_pivot.fillna(0)
    
    # Sort by total citations
    category_source_pivot['total'] = category_source_pivot.sum(axis=1)
    category_source_pivot = category_source_pivot.sort_values('total', ascending=False)
    category_source_pivot = category_source_pivot.drop('total', axis=1)
    
    # Keep only known sources
    category_source_pivot = category_source_pivot[KNOWN_SOURCES]
    
    # Set up the plot
    plt.figure(figsize=(12, max(8, len(category_source_pivot) * 0.4)))
    
    # Generate the heatmap
    ax = sns.heatmap(
        category_source_pivot,
        cmap="YlGnBu",
        annot=True,
        fmt=".0f",
        linewidths=0.5,
        cbar_kws={"label": "Total Citations"}
    )
    
    # Customize the plot
    plt.title("Citation Sources by Category (Total Citations)")
    plt.tight_layout()
    
    # Save the figure
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    logger.info(f"Heatmap saved to {output_path}")
    
    # Close the figure to free memory
    plt.close()


def generate_source_distribution_chart(df: pd.DataFrame, output_path: str):
    """
    Generate a bar chart showing the distribution of citation sources.
    
    Args:
        df: DataFrame with citation source data
        output_path: Path to save the output image
    """
    if df.empty:
        logger.error("No data available for chart generation")
        return
    
    # Count tools using each citation source
    citation_source_counts = df['Citation Source'].value_counts()
    yearly_source_counts = df['Yearly Citation Source'].value_counts()
    
    # Combine into a new dataframe for plotting
    source_data = []
    
    for source in KNOWN_SOURCES:
        source_data.append({
            'Source': source,
            'Citation Count': citation_source_counts.get(source, 0),
            'Yearly Data': yearly_source_counts.get(source, 0)
        })
    
    source_df = pd.DataFrame(source_data)
    
    # Set up the plot
    plt.figure(figsize=(10, 6))
    
    # Create the grouped bar chart
    x = np.arange(len(source_df))
    width = 0.35
    
    plt.bar(x - width/2, source_df['Citation Count'], width, label='Citation Count Source')
    plt.bar(x + width/2, source_df['Yearly Data'], width, label='Yearly Data Source')
    
    # Customize the plot
    plt.xlabel('Citation Source')
    plt.ylabel('Number of Tools')
    plt.title('Distribution of Citation Sources')
    plt.xticks(x, source_df['Source'])
    plt.legend()
    
    # Add value labels on top of each bar
    for i, v in enumerate(source_df['Citation Count']):
        plt.text(i - width/2, v + 0.5, str(v), ha='center')
    
    for i, v in enumerate(source_df['Yearly Data']):
        plt.text(i + width/2, v + 0.5, str(v), ha='center')
    
    plt.tight_layout()
    
    # Save the figure
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    logger.info(f"Source distribution chart saved to {output_path}")
    
    # Close the figure to free memory
    plt.close()


def main():
    """Main function to run the heatmap generator."""
    parser = argparse.ArgumentParser(description="Generate citation source heatmap visualization")
    parser.add_argument(
        "--output-dir",
        default=str(OUTPUT_DIR),
        help="Directory to save output files (defaults to repository root)"
    )
    args = parser.parse_args()
    
    # Load impact data
    impact_data = load_impact_data()
    if not impact_data:
        logger.error("Failed to load impact data")
        return 1
    
    # Extract source data
    df, categories = get_source_data(impact_data)
    if df.empty:
        logger.error("No citation data found in impact_data.json")
        return 1
    
    # Ensure output directory exists
    output_dir = Path(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)
    
    # Generate visualizations
    heatmap_path = output_dir / "citation_source_heatmap.png"
    generate_heatmap(df, str(heatmap_path))
    
    distribution_path = output_dir / "citation_source_distribution.png"
    generate_source_distribution_chart(df, str(distribution_path))
    
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())