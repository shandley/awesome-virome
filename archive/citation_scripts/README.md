# Archived Citation Scripts

**Archived Date:** 2025-11-29

## Why These Scripts Were Archived

These Python scripts powered the citation analytics system that was removed from the project. The citation tracking infrastructure relied on multiple external APIs that proved unreliable:

- **CrossRef API** - Publication citation data (rate limits, incomplete data)
- **PubMed API** - Medical literature citations (limited coverage of bioinformatics tools)
- **Semantic Scholar API** - Academic paper citations (inconsistent results)

The maintenance burden and API unreliability outweighed the value provided by citation metrics.

## Archived Scripts

### Citation Data Collection
- **academic_impact.py** - Main script for collecting academic impact metrics from multiple sources
- **pubmed_citations.py** - PubMed citation data collector
- **comprehensive_citation_data.py** - Comprehensive citation data aggregator from all sources
- **real_citation_data.py** - Real-time citation data fetcher
- **direct_citation_collector.py** - Direct citation relationship collector using CrossRef

### Citation Reporting
- **citation_report.py** - Generated citation analytics reports
- **dashboard_citation_data.py** - Prepared citation data for dashboard visualization
- **generate_citation_data.py** - Citation data generator for static pages

### Citation Analysis
- **citation_source_heatmap.py** - Visualization of citation sources and relationships
- **validate_citations.py** - Citation data validation and quality checks

### Impact Tracking
- **update_impact_data.py** - Updated impact_data.json with citation metrics
- **simplified_update_impact_data.py** - Simplified version for workflow automation

## Dependencies Removed

These scripts required additional dependencies that are no longer needed:
- Specific API client libraries for CrossRef, PubMed, Semantic Scholar
- Citation data parsing and normalization libraries
- Network analysis libraries for citation relationship graphs

## What Replaced This

The project now relies on:
- **scripts/github_metrics_workflow.py** - Reliable GitHub API-based metrics
- **scripts/enhance_metadata.py** - Enhanced metadata from repository sources
- **scripts/bioinformatics_metadata.py** - Bioinformatics-specific metadata

## Related Changes

- Citation HTML pages archived to `archive/citation_features/`
- Citation workflows disabled in `.github/workflows/disabled/`
- `impact_data.json` generation removed from workflows
- Citation-related flags removed from `update_data_json.py`

## Preservation Note

These scripts are preserved for historical reference. They could be restored if:
- More reliable citation APIs become available
- Alternative approaches to citation tracking are developed
- The community specifically requests citation analytics features
- DOI-based citation tracking becomes more robust

## Historical Context

The citation system went through two design iterations:
1. **Version 1** - Basic citation counting from CrossRef
2. **Version 2** (Partially Implemented) - Multi-source citation aggregation with PubMed, Semantic Scholar, and CrossRef

Version 2 Phase 1 (GitHub metrics enhancement) was successfully deployed and remains active. Phases 2-4 (citation tracking) were abandoned due to API reliability issues.
