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

Metrics are automatically integrated with the main GitHub Pages site. You can view them at:

`https://[org-name].github.io/awesome-virome/metrics_dashboard/`

A link to the metrics dashboard is automatically added to the main site's homepage.

## How It Works

1. After each workflow run, the metrics-history job executes in the repository-testing workflow
2. It checks out the metrics-history branch
3. Collects data from the workflow run
4. Updates JSON files with the latest metrics
5. Generates new charts and visualizations
6. Commits and pushes changes back to the metrics-history branch

7. When the GitHub Pages workflow runs (either automatically or manually), it:
   - Checks out the main branch for the primary content
   - Checks out the metrics-history branch into a temporary directory
   - Copies the metrics dashboard files to a `/metrics_dashboard` folder in the main content
   - Adds a link to the metrics dashboard in the main site's homepage
   - Deploys everything together as a single GitHub Pages site

## Setting Up Metrics Tracking

To set up metrics tracking for the first time:

1. Run the repository testing workflow with the "Update metrics" option enabled
   - This will automatically create the `metrics-history` branch if it doesn't exist
   - The initial metrics will be collected and stored

2. The modified GitHub Pages workflow will:
   - Automatically detect the metrics-history branch
   - Copy the metrics dashboard to the main site
   - Deploy everything together

3. No separate GitHub Pages configuration is needed for the metrics dashboard

The integrated approach allows for a single GitHub Pages deployment while keeping the metrics data stored in a separate branch for better organization.

The URL for the metrics dashboard will be: `https://[org-name].github.io/awesome-virome/metrics_dashboard/`

> Note: The metrics history system is designed to be self-initializing. Even if you don't have any data at first, it will start tracking metrics from the first workflow run where it's enabled.