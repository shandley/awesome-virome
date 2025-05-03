# Repository Health Metrics

This directory contains historical repository health metrics and visualizations that track the overall health, activity, and freshness of the awesome-virome repository.

## Directory Structure

- `data/` - Contains raw data files used to generate metrics
  - `git_stats.json` - Basic Git repository statistics
  - `freshness.json` - Data freshness metrics for different repository components
  - `activity_metrics.json` - Code activity metrics (commits, additions, deletions)
  - `contribution_metrics.json` - Contributor analysis metrics
  - Various CSV files with raw data for visualization

- `charts/` - Contains visualization charts generated from the metrics data
  - `monthly_commits.svg/.png` - Monthly commit activity visualization
  - `code_churn.svg/.png` - Code additions and deletions over time 
  - `top_contributors.svg/.png` - Top contributors by commit count
  - `contribution_distribution.svg/.png` - Distribution of contributions
  - `data_freshness.svg/.png` - Freshness of key data components

- `metrics_summary.json` - Summary of current repository health metrics

## Metrics Collection

These metrics are automatically collected by the Repository Health Metrics workflow which runs:
- Weekly (Sunday at midnight UTC)
- After major data updates
- When manually triggered

The metrics help maintainers track:
1. Code activity and contribution patterns
2. Data freshness and update frequency
3. Repository growth trends
4. Contribution distribution

## Visualizations

Visualizations are automatically generated to make the metrics more accessible. These include:
- Time-series charts showing activity over time
- Bar charts showing top contributors
- Pie charts showing contribution distribution
- Gauge charts showing data freshness

## Health Report

A comprehensive health report that includes these metrics and visualizations is available in the `reports/repository_health.md` file.