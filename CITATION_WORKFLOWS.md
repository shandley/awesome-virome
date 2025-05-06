# Citation Data Collection and Management

This document provides an overview of how citation data is collected, validated, and managed in the awesome-virome repository.

## Overview

The citation system collects real citation data from published papers about virome analysis tools. This data is used to:

1. Provide accurate citation metrics for each tool
2. Show tool impact in the citation dashboard
3. Support analysis of the virome tool ecosystem

## Key Principles

- **Real data only**: All citation data is collected from real sources (PubMed, CrossRef, etc.) with no synthetic data generation
- **Single source of truth**: Each task (data collection, DOI fixing, etc.) has one authoritative source
- **Maintainable workflows**: Clear separation of responsibilities between workflows

## Authoritative Components

### Scripts

1. **comprehensive_citation_data.py** - The authoritative script for generating `impact_data.json`
   - Collects data from multiple sources: academic_impact, bioinformatics, and citation reports
   - Uses only real citation data with no synthetic data generation
   - Consolidates all available citation information

2. **auto_fix_dois.py** - The authoritative script for DOI format fixing
   - Used by the `fix-dois.yml` workflow
   - Automatically fixes common DOI format issues

### Workflows

1. **Citation Data Collection and Impact Data Generation** (pubmed-citations.yml)
   - Collects citation data from PubMed and other sources
   - Generates citation reports
   - Creates the authoritative `impact_data.json` file using comprehensive_citation_data.py

2. **DOI Format Fixing** (fix-dois.yml)
   - Runs the authoritative DOI fixing script auto_fix_dois.py
   - Fixes common DOI formatting issues across all metadata files

3. **Citation Data Validation** (citation-validation-direct.yml)
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
   - This file is used by the dashboard for visualization

## Recent Improvements

The citation data system has been refactored to:

1. **Remove synthetic data generation**:
   - Eliminated artificial citation distribution algorithms
   - Removed percentage-based "influential citations" estimates
   - Removed heuristic tool categorization

2. **Consolidate responsibilities**:
   - Made comprehensive_citation_data.py the single source of truth for impact_data.json
   - Consolidated DOI fixing to auto_fix_dois.py and its workflow
   - Removed redundant impact data generation from other workflows

3. **Simplify workflow structure**:
   - Removed deprecated workflows
   - Updated remaining workflows to have clear, non-overlapping responsibilities
   - Added clear naming to indicate authoritative sources

## Maintenance Guidelines

When making changes to the citation system:

1. **Follow the authoritative pattern**:
   - Use the designated scripts for their specific purposes
   - Don't create new scripts that duplicate functionality

2. **Avoid synthetic data**:
   - Never generate synthetic citation data
   - If data is missing, leave it empty rather than estimating

3. **Keep workflows focused**:
   - Each workflow should have a clear, non-overlapping responsibility
   - workflows.md contains details about each workflow's purpose and schedule

## Testing Before Pushing Changes

Before pushing changes to GitHub, run the test suite to verify that workflows function correctly and don't generate synthetic data:

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
1. Scripts don't generate synthetic data
2. DOI fixing works correctly
3. Workflow files don't have syntax errors
4. Workflows can run without errors
5. impact_data.json contains only real data

If any test fails, fix the issues before pushing to GitHub.