# Simplified GitHub Workflows

This directory contains simplified GitHub Actions workflows for the Awesome-Virome repository. The previous complex workflow system has been replaced with a streamlined set of workflows.

## Active Workflows

1. **simplified-update-workflow.yml**
   - Consolidates all update operations (repo info, citations, metadata)
   - Runs on schedule or manual trigger
   - Directly commits changes to main branch
   - No PR creation or approval process needed
   - Replaces multiple specialized workflows (update-repos.yml, update-citation-data.yml, etc.)

2. **simplified-pages-deploy.yml**
   - Deploys to GitHub Pages when data.json changes
   - Handles GitHub Pages configuration and deployment
   - Automatically triggered by changes to data.json
   - Replaces previous custom-pages-deploy.yml and simple-pages-deploy.yml

3. **contributor-suggestion.yml**
   - Processes contributions from external users
   - Validates suggested changes
   - Provides automatic feedback on PRs
   - Ensures quality of external contributions

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