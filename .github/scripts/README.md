# Repository Metrics Tracking System

## Overview

This directory contains scripts that power the repository metrics tracking system, which collects and visualizes important metrics about the repository health over time.

## Metrics Tracked

- **Performance metrics**: Load times for data files and parsing times for README
- **Validation success rate**: Percentage of validation tests that pass
- **Link health**: Tracking of broken vs working links over time

## Scripts

### `store_metrics_history.py`

- Collects metrics from GitHub Actions workflow runs
- Stores them in a structured JSON format
- Commits them to the `metrics-history` branch

### `generate_metrics_charts.py`

- Reads metrics data from the history files
- Generates SVG charts visualizing trends
- Creates an HTML dashboard for viewing metrics

## How to View Metrics

Metrics are automatically published to GitHub Pages from the `metrics-history` branch. You can view them at:

`https://[org-name].github.io/awesome-virome/metrics/`

## How It Works

1. After each workflow run, the metrics-history job executes
2. It checks out the metrics-history branch
3. Collects data from the workflow run
4. Updates JSON files with the latest metrics
5. Generates new charts and visualizations
6. Commits and pushes changes back to the metrics-history branch
7. GitHub Pages automatically deploys the updated dashboard

## Setting Up Metrics Tracking

To set up metrics tracking for the first time:

1. Run the repository testing workflow with the "Update metrics" option enabled
   - This will automatically create the `metrics-history` branch if it doesn't exist
   - The initial metrics will be collected and stored

2. Enable GitHub Pages for the repository:
   - Go to repository Settings > Pages
   - Set the Source to "Deploy from a branch"
   - Select the `metrics-history` branch
   - Set the folder to `/metrics_history/charts`
   - Click Save

3. After GitHub Pages is enabled, metrics visualizations will be automatically published with each workflow run

The URL for the metrics dashboard will be: `https://[org-name].github.io/awesome-virome/`

> Note: The metrics history system is designed to be self-initializing. Even if you don't have any data at first, it will start tracking metrics from the first workflow run where it's enabled.