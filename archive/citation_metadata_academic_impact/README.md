# Archived Academic Impact Citation Metadata

**Archived Date:** 2025-11-29

## About This Directory

This directory contains academic impact and citation metadata collected from multiple external APIs. The data was used to power citation analytics features that have been removed from the project.

**Size:** ~3.2MB (298 JSON files)

## Data Sources

The metadata in this directory was collected from:
- **CrossRef API** - Publication citation counts and relationships
- **Semantic Scholar API** - Academic paper citations and influence metrics
- **PubMed API** - Medical/biological literature citations

## File Structure

Each JSON file contains citation metadata for a specific tool or category:
- Individual tool files (e.g., `PHANOTATE.json`, `PhageBoost.json`)
- Category aggregation files (e.g., `Host Prediction.json`, `Genome Assembly.json`)
- Master aggregation file (`academic_impact.json`)

## Why This Was Archived

This citation metadata collection was discontinued because:
1. **API Reliability Issues** - Frequent rate limits, timeouts, and incomplete data
2. **Maintenance Burden** - Required constant monitoring and error handling
3. **Data Quality** - Inconsistent citation counts across different sources
4. **Limited Value** - GitHub metrics (stars, forks) proved more reliable indicators of tool adoption

## Related Archives

- `archive/citation_features/` - HTML pages that visualized this data
- `archive/citation_scripts/` - Python scripts that collected this data
- `archive/citation_metadata_pubmed/` - PubMed-specific citation metadata
- `.github/workflows/disabled/` - Citation data collection workflows

## Preservation Note

This data is preserved for historical reference and potential future use if:
- Citation API reliability significantly improves
- Alternative data sources become available
- Research projects need historical citation data
