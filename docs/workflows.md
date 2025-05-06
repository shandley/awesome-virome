# GitHub Workflows

This document provides details about the GitHub Actions workflows used in the awesome-virome repository.

## Citation Data Workflows

### Citation Data Collection and Impact Data Generation

**File**: `.github/workflows/pubmed-citations.yml`  
**Schedule**: Daily at 04:00 UTC  
**Manual Trigger**: Yes  

**Purpose**:
- Collects citation data from PubMed using NCBI E-utilities API
- Generates citation reports and analysis
- Runs the authoritative script for generating `impact_data.json`

**Key Steps**:
1. Collects PubMed citation data using `pubmed_citations.py`
2. Generates citation reports using `citation_report.py`
3. Updates `impact_data.json` using the authoritative `comprehensive_citation_data.py`

### DOI Format Fixing

**File**: `.github/workflows/fix-dois.yml`  
**Schedule**: Weekly on Sunday at 08:00 UTC  
**Manual Trigger**: Yes  

**Purpose**:
- Automatically fixes common DOI formatting issues across metadata files
- Consolidates DOI fixing functionality in one place

**Key Steps**:
1. Runs the authoritative `auto_fix_dois.py` script
2. Commits and pushes any changes to DOI formats

### Citation Data Validation

**File**: `.github/workflows/citation-validation-direct.yml`  
**Schedule**: Weekly on Sunday at 12:00 UTC  
**Manual Trigger**: Yes  

**Purpose**:
- Validates DOI consistency across metadata files
- Generates validation reports
- Creates issues for critical validation problems

**Key Steps**:
1. Runs citation validation script (without DOI fixing)
2. Generates validation reports
3. Checks for critical issues
4. Creates GitHub issues for critical problems if found

### Update Citation Counts from CrossRef

**File**: `.github/workflows/update-citation-counts.yml`  
**Schedule**: Weekly on Sunday at 03:00 AM UTC  
**Manual Trigger**: Yes  

**Purpose**:
- Updates citation counts from CrossRef API
- Does NOT generate `impact_data.json` (delegated to the authoritative workflow)

**Key Steps**:
1. Fixes citation counts using CrossRef data
2. Updates citation count data in metadata files
3. Commits changes to metadata files only

## Disabled/Deprecated Workflows

The following workflows have been deprecated or disabled in favor of the authoritative workflows above:

### Citation Data Validation (DEPRECATED)

**File**: `.github/workflows/disabled/citation-validation.yml`  
**Status**: Moved to disabled directory  
**Replacement**: Use `citation-validation-direct.yml` instead

### Citation Data Consolidation (DEPRECATED)

**File**: `.github/workflows/update-impact-data.yml`  
**Status**: Deprecated but not removed (accessible via manual trigger only)  
**Replacement**: Use the authoritative Citation Data Collection workflow instead

## Other Workflows

### Site Health Check

**File**: `.github/workflows/site-health-check.yml`  
**Schedule**: Daily at 01:00 UTC  
**Purpose**: Checks the health of the website, verifying that all pages load correctly

### Repository Health Metrics

**File**: `.github/workflows/repository-health-metrics.yml`  
**Schedule**: Weekly on Monday  
**Purpose**: Collects and reports on repository health metrics

### Broken Link Checker

**File**: `.github/workflows/broken-link-checker.yml`  
**Schedule**: Weekly  
**Purpose**: Checks for broken links throughout the site

### Auto Changelog

**File**: `.github/workflows/automated-changelog.yml`  
**Trigger**: On pull request merge to main  
**Purpose**: Automatically updates the changelog when PRs are merged