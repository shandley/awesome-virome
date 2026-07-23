# Simplified GitHub Workflows

This directory contains simplified GitHub Actions workflows for the Awesome-Virome repository. The previous complex workflow system has been replaced with a streamlined set of workflows.

## Active Workflows

1. **simplified-update-workflow.yml**
   - Pulls GitHub metadata (stars, forks, language, topics), regenerates data.json, and generates the API
   - Runs on schedule (weekly and monthly) or manual trigger
   - Commits refreshed data directly to main

2. **unified-pages-deploy.yml**
   - Builds and publishes the site to the gh-pages branch (the Pages source)
   - Deploys the MkDocs docs and copies dashboard/comparison/landing-page assets
   - Triggered on push to main that touches site content

3. **broken-link-checker.yml**
   - Checks links in README/API/CONTRIBUTING and maintains a single rolling broken-links issue
   - Runs weekly and after the data update

4. **validate-contribution.yml**
   - Validates tool submissions from issues (submission template) and pull requests
   - Comments quality feedback and applies validated / needs-changes labels

See `docs/workflows.md` for the full description of each.

## Disabled Workflows

All previous workflows have been moved to the `disabled/` directory for reference. These include:
- update-repos.yml
- update-citation-data.yml
- sync-json-from-readme.yml
- cache-related workflows
- and others

These disabled workflows are kept for historical reference but are no longer active.

## How the Simplified Process Works

1. **Data Management**:
   - data.json is the single source of truth
   - README.md content is synchronized with data.json
   - Metadata is stored in the metadata/ directory

2. **Update Process**:
   - Scheduled updates run weekly (basic) and monthly (comprehensive)
   - Updates are committed directly to main branch
   - GitHub Pages automatically deploy with new data

3. **Contribution Process**:
   - External contributors create PRs with changes to README or metadata
   - Validation workflow checks their changes
   - Maintainers review and merge valid contributions

## Implementation Notes

This simplified approach eliminates the complex PR-based workflow and branch management previously used. Instead, we rely on direct commits to main for automatic updates while still supporting external contributions through the standard PR process.

Benefits:
- Reduced complexity and fewer moving parts
- Fewer permissions issues with GitHub Actions
- More straightforward management of the repository
- Clearer separation between automated updates and human contributions