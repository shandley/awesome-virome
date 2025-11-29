# Disabled Workflows

This directory contains workflows that have been disabled or replaced by more streamlined, authoritative workflows.

## Why Workflows Were Disabled

As part of the citation workflow refactoring and repository cleanup:

1. Consolidated responsibilities into clear authoritative workflows
2. Removed synthetic data generation
3. Eliminated redundant functionality
4. Archived deprecated Citation System V1 (premium APIs)
5. Removed references to missing/archived scripts

## Recently Disabled (November 2025 Cleanup)

### `update-citation-counts.yml`
**Reason**: References missing script `update_citation_counts.py` and archived script `fix_citation_counts.py`
**Note**: Workflow itself notes it's "NOT Authoritative for impact_data.json" - functionality consolidated

### `scopus-test.yml` & `wos-test.yml`
**Reason**: Test workflows for premium APIs (Scopus, Web of Science) from deprecated Citation System V1
**Note**: Project moved to free APIs only (GitHub API, CrossRef, etc.) in Citation System V2

### `test-citation-system.yml`
**Reason**: Tests the archived `scripts/citation_system/` package (Citation System V1)
**Note**: Tests Scopus, WOS, and iCite integration - replaced by simpler Phase 1 approach

## Previously Disabled Workflows

For the workflows in this directory, use these authoritative replacements instead:

| Disabled Workflow | Authoritative Replacement |
|-------------------|---------------------------|
| `auto-fix-dois.yml` | `fix-dois.yml` (DOI Format Fixing) |
| `citation-validation.yml` | `citation-validation-direct.yml` (Citation Data Validation) |
| `dashboard-update-workflow.yml` | `pubmed-citations.yml` (Citation Data Collection) |
| `update-dashboard-data.yml` | `pubmed-citations.yml` (Citation Data Collection) |
| `update-citation-data.yml` | `simplified-update-workflow.yml` (Main Update) |
| `update-dashboard-with-real-data.yml` | `pubmed-citations.yml` (Citation Data Collection) |

## Re-enabling Workflows

To re-enable a workflow:

1. **Fix the underlying issues** - Create missing scripts or restore archived code
2. **Move the file back**: `mv .github/workflows/disabled/workflow-name.yml .github/workflows/`
3. **Test locally** if possible before re-enabling
4. **Update documentation** to reflect the change

## Active Core Workflows

The main active workflows are:
- `simplified-update-workflow.yml` - Core update workflow (weekly/monthly)
- `github_metrics_workflow.py` - Citation System V2 Phase 1 (GitHub metrics)
- `cache-maintenance.yml` - Cache system maintenance
- `simplified-pages-deploy.yml` - GitHub Pages deployment

For more information on the citation workflow architecture, see `/CITATION_WORKFLOWS.md` in the repository root.