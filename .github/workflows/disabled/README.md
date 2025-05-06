# Disabled Workflows

This directory contains workflows that have been disabled or replaced by more streamlined, authoritative workflows.

## Why Workflows Were Disabled

As part of the citation workflow refactoring, we have:

1. Consolidated responsibilities into clear authoritative workflows
2. Removed synthetic data generation
3. Eliminated redundant functionality

## Authoritative Replacement Workflows

For the workflows in this directory, use these authoritative replacements instead:

| Disabled Workflow | Authoritative Replacement |
|-------------------|---------------------------|
| `auto-fix-dois.yml` | `fix-dois.yml` (DOI Format Fixing) |
| `citation-validation.yml` | `citation-validation-direct.yml` (Citation Data Validation) |
| `dashboard-update-workflow.yml` | `pubmed-citations.yml` (Citation Data Collection) |
| `update-dashboard-data.yml` | `pubmed-citations.yml` (Citation Data Collection) |
| `update-citation-data.yml` | `update-citation-counts.yml` (Update Citation Counts) |
| `update-dashboard-with-real-data.yml` | `pubmed-citations.yml` (Citation Data Collection) |

For more information on the citation workflow architecture, see `/CITATION_WORKFLOWS.md` in the repository root.