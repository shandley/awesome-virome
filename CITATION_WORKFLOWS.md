# Citation Data Collection and Management

This document provides an overview of how citation data is collected, validated, and managed in the awesome-virome repository.

## Overview

The citation system collects real citation data from published papers about virome analysis tools. This data is used to:

1. Provide accurate citation metrics for each tool
2. Show tool impact in the citation dashboard
3. Support analysis of the virome tool ecosystem

## Citation Systems

The repository has two citation collection systems:

### 1. Real Data System

- **Real data only**: All citation data is collected from real sources (PubMed, CrossRef, etc.)
- **Single source of truth**: Each task (data collection, DOI fixing, etc.) has one authoritative source
- **Maintainable workflows**: Clear separation of responsibilities between workflows

### 2. Hybrid System (New)

- **Prioritizes real data**: Uses real citation data from PubMed, CrossRef, etc. when available
- **Transparent source attribution**: Clearly indicates which source provided each citation count 
- **Visualization support**: For tools without real citation data, generates synthetic data for visualization purposes only
- **Clear marking**: All synthetic data is clearly marked with `citation_source: "synthetic"`

## Citation Sources

Citations are collected from the following sources, in order of priority:

1. **PubMed** - Citation data from PubMed Central and PubMed, stored in `metadata/pubmed_citations/*.json`
2. **iCite** - NIH's iCite API provides citation data for papers with PMIDs
3. **Scopus** - When configured with API credentials, provides citation data from Elsevier's Scopus database
4. **Web of Science** - When configured with API credentials, provides citation data from Clarivate's Web of Science
5. **CrossRef** - Provides basic citation data for papers with DOIs

## Attribution Fields

Each tool in `impact_data.json` includes the following attribution fields:

- `citation_source`: The source of the total citation count (e.g., "pubmed", "icite", "scopus", "synthetic")
- `yearly_citation_source`: The source of yearly citation breakdown data
- `total_citations`: Total citation count from the specified source
- `citations_by_year`: Yearly breakdown of citations, if available

## Authoritative Components

### Scripts

1. **comprehensive_citation_data.py** - The authoritative script for generating `impact_data.json` with real data only
   - Collects data from multiple sources: academic_impact, bioinformatics, and citation reports
   - Uses only real citation data with no synthetic data generation
   - Consolidates all available citation information

2. **update_citation_data.py** - The hybrid approach script for citation data
   - Combines real citation data with synthetic data for visualization
   - Clearly attributes all data to its source
   - Prioritizes real data over synthetic data

3. **add_citation_attribution.py** - Adds source attribution to citation data
   - Ensures all citation data has proper source attribution
   - Works with both real and synthetic data

4. **citation_source_heatmap.py** - Generates visualizations of citation sources
   - Creates heatmaps showing which citation sources were used for which tools
   - Creates distribution charts of citation sources

5. **auto_fix_dois.py** - The authoritative script for DOI format fixing
   - Used by the `fix-dois.yml` workflow
   - Automatically fixes common DOI format issues

### Workflows

1. **Citation Data Collection and Impact Data Generation** (pubmed-citations.yml)
   - Collects citation data from PubMed and other sources
   - Generates citation reports
   - Creates the authoritative `impact_data.json` file using comprehensive_citation_data.py

2. **Hybrid Citation System** (hybrid-citation-system.yml)
   - Runs update_citation_data.py to combine real and synthetic data
   - Adds source attribution with add_citation_attribution.py
   - Generates visualizations with citation_source_heatmap.py
   - Commits changes with clear indication of data sources

3. **DOI Format Fixing** (fix-dois.yml)
   - Runs the authoritative DOI fixing script auto_fix_dois.py
   - Fixes common DOI formatting issues across all metadata files

4. **Citation Data Validation** (citation-validation-direct.yml)
   - Validates DOI consistency and data quality
   - Creates validation reports
   - Does NOT modify DOIs (this is handled by the DOI Format Fixing workflow)

## Data Flow

1. **Citation Collection**:
   - Citation data is collected from PubMed, CrossRef, and other sources
   - Data is stored in JSON files in `metadata/pubmed_citations/`

2. **DOI Validation and Fixing**:
   - DOIs are validated for consistency and format
   - Any format issues are fixed by the DOI Format Fixing workflow
   - Validation reports are generated in `reports/citations/`

3. **Impact Data Generation**:
   - All citation data is consolidated into `impact_data.json`
   - With the hybrid approach, synthetic data is added for visualization purposes
   - All data points are attributed to their source
   - This file is used by the dashboard for visualization

## Visualization

Citation source data is visualized in two primary ways:

1. **Citation Source Heatmap** (`citation_source_heatmap.png`) - Shows citation sources by category
2. **Citation Source Distribution** (`citation_source_distribution.png`) - Shows the distribution of citation sources across all tools

## Choosing the Right System

### When to use the Real Data System
- For scientific analysis of citation impact
- When transparency about data limitations is most important
- When you want to guarantee all data comes from academic sources

### When to use the Hybrid System
- For comprehensive visualizations that include all tools
- When a complete picture is more important than 100% real data
- For dashboards and user interfaces that need data for all tools

## Testing Before Pushing Changes

Before pushing changes to GitHub, run the test suite to verify that workflows function correctly:

```bash
# Run all tests
python scripts/tests/run_all_tests.py

# Run individual tests
python scripts/tests/test_comprehensive_citation_data.py
python scripts/tests/test_auto_fix_dois.py
python scripts/tests/validate_workflows.py
python scripts/tests/run_workflows_dryrun.py
python scripts/tests/verify_impact_data.py
```

The test suite verifies:
1. Scripts function correctly
2. DOI fixing works correctly
3. Workflow files don't have syntax errors
4. Workflows can run without errors
5. impact_data.json contains correctly attributed data