"""
Create a timeline image from the json file
"""

import os
import sys
import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

__author__ = 'Rob Edwards'

def create_timeline_image(data, outputfile, verbose=False):
    df = pd.read_json(data, orient='records')
    # Convert the date columns to datetime
    df['start_date'] = pd.to_datetime(df['created_at'], format='mixed', utc=True)
    df['end_date'] = pd.to_datetime(df['last_updated'], format='mixed', utc=True)
    df = df.dropna(subset=['start_date', 'end_date'])
    df = df.sort_values(by='start_date')

    # generate a colour map
    palette = sns.color_palette("husl", len(df['classification'].unique()))
    color_map = {classification: palette[i] for i, classification in enumerate(df['classification'].unique())}
    
    # Create the plot
    plt.figure(figsize=(24, df.shape[0]/4)) # this adds 1/4 inch per row
    for i, row in df.iterrows():
        plt.plot([row['start_date'], row['end_date']], [row['name'], row['name']], marker='o')

    plt.grid(axis='y', linewidth=0)
    # Customize the plot
    plt.xlabel('')
    plt.ylabel('')
    plt.title('Virus bioinformatics tools')
    plt.grid(True)
    # Add a legend
    handles = [plt.Line2D([0], [0], color=color_map[cat], lw=2) for cat in sorted(color_map)]
    labels = [cat for cat in sorted(color_map)]
    plt.legend(handles, labels, title='Category')

    plt.tight_layout()
    # Show the plot
    plt.savefig(outputfile)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=' ')
    parser.add_argument('-j', '--json', help='json file with the repo dates', default='repo_updates.json')
    parser.add_argument('-o', '--output', help='output file', default='timeline.png')
    parser.add_argument('-v', '--verbose', help='verbose output', action='store_true')
    args = parser.parse_args()

    create_timeline_image(args.json, args.output, args.verbose)







