# GitHub Workflows

This document describes the GitHub Actions workflows used in the awesome-virome repository. The set was reduced in July 2026 to the four workflows that each do a distinct job; the citation, DOI, health-metrics, cache-maintenance, site-health, changelog, and duplicate pages-deploy workflows were removed.

## Active workflows

### Simplified Data Update Workflow

**File**: `.github/workflows/simplified-update-workflow.yml`
**Schedule**: Weekly (Sunday 08:00 UTC) and monthly (1st of the month, 08:00 UTC)
**Manual trigger**: Yes

This is the core data engine. It pulls repository metadata from the GitHub API (stars, forks, language, topics), regenerates `data.json`, and generates the REST API endpoints under `api/`. It commits the refreshed data back to `main`.

Key scripts: `update_check.py`, `github_metrics_workflow.py`, `enhance_metadata.py`, `bioinformatics_metadata.py`, `update_data_json.py`, `generate_api.py`, `verify_readme_content.py`.

### Unified GitHub Pages Deployment

**File**: `.github/workflows/unified-pages-deploy.yml`
**Trigger**: Push to `main` touching `data.json`, `*.html`, `metadata/**`, `reports/**`, `js/**`, `mkdocs/**`, or the landing page
**Manual trigger**: Yes

Builds and publishes the site to the `gh-pages` branch, which is the source GitHub Pages serves from. It deploys the MkDocs docs (via `mike`) and copies the dashboard, comparison, selection-guide, landing page, and data assets into the branch. This is the only site deployer; the former `simplified-pages-deploy.yml` used the Actions-artifact deploy method, which does not apply when Pages serves from a branch, so it was removed.

### Broken Link Checker

**File**: `.github/workflows/broken-link-checker.yml`
**Schedule**: Weekly (Monday 00:00 UTC), and after the data update workflow
**Manual trigger**: Yes (with an optional `check_all` input to also scan HTML and JSON)

Checks the links in `README.md`, `API.md`, and `CONTRIBUTING.md`, and opens or updates a single rolling issue labeled `broken-links` when it finds problems. It skips entries already tagged `[unavailable]` and web.archive.org snapshots, allowlists domains that block automated checkers but work in a browser (publishers, DOI, Zenodo, SourceForge), treats 401/403/405/429/503 as "blocked" rather than broken, and retries a connection failure once before declaring a link dead.

### Validate Tool Contribution

**File**: `.github/workflows/validate-contribution.yml`
**Trigger**: New or edited issues using the tool-submission template, and pull requests touching `README.md`, `data.json`, or `metadata/**`
**Manual trigger**: Yes (by issue or PR number)

Validates tool submissions and posts feedback. For issues it parses the submission template, runs `tool_validator.py`, and comments a quality report with a maintainer checklist, applying `validated` / `needs-changes` labels. For pull requests it validates the changed tool data or README (with duplicate checking) and comments the result. This replaced the separate `contributor-suggestion.yml`, which duplicated the PR-validation path.

## Disabled workflows

`.github/workflows/disabled/` holds older experimental workflows that are not run. They are kept for reference only and are not wired to any trigger.
