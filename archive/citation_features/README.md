# Archived Citation Features

**Archived Date:** 2025-11-29

## Why These Files Were Archived

These HTML files provided citation analytics and publication impact visualization features that were removed from the project. The citation tracking system relied on external APIs (CrossRef, PubMed, Semantic Scholar) that proved unreliable and required ongoing maintenance.

## Archived Files

- **citations.html** - Main citation analytics dashboard showing growth trends and impact metrics
- **publication_impact.html** - Publication impact visualization and network analysis
- **direct_publication_impact.html** - Tool citation network visualization using CrossRef data

## What Replaced This

The project now focuses on:
- Interactive tool dashboard (dashboard.html) with categorization and timeline visualizations
- Tool comparison matrix (comparison.html)
- Tool selection guide (selection-guide.html)
- GitHub-based metrics (stars, forks, activity) which are more reliable than citation APIs

## Related Changes

- Citation scripts archived to `archive/citation_scripts/`
- Citation workflows disabled in `.github/workflows/disabled/`
- Citation references removed from:
  - dashboard.html
  - README.md
  - .github/workflows/simplified-update-workflow.yml

## Preservation Note

These files are preserved for historical reference and could be restored if:
- Citation API reliability improves
- Alternative citation data sources become available
- Community requests citation tracking features
