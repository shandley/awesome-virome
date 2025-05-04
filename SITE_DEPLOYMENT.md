# GitHub Pages Deployment Structure

This document explains how the Awesome-Virome repository is deployed to GitHub Pages, with both documentation and an interactive dashboard.

## Site Structure

The site is deployed with the following structure:

- **Root (`/`)**: Landing page with links to both the documentation and dashboard
- **Dashboard (`/dashboard.html`)**: Interactive dashboard for exploring virome tools
- **Documentation (`/1.0.0/`)**: MkDocs documentation site (versioned)

## Deployment Workflow

We use a unified GitHub Pages deployment workflow (`/.github/workflows/unified-pages-deploy.yml`) that:

1. Builds the MkDocs site using the `mike` tool for versioning
2. Maintains the dashboard by copying all relevant HTML, JS, and data files
3. Creates a landing page that provides access to both components

## How It Works

The workflow:

1. Checks out the repository
2. Installs Python and MkDocs dependencies
3. Builds the versioned MkDocs site and deploys it to the `gh-pages` branch
4. Checks out the `gh-pages` branch (which now contains the MkDocs site)
5. Copies all dashboard files from the main branch to the `gh-pages` branch
6. Adds a landing page as the main `index.html`
7. Uploads and deploys the contents of the `gh-pages` branch to GitHub Pages

## Benefits

This approach:

- Keeps the MkDocs documentation properly versioned
- Maintains the interactive dashboard functionality
- Provides a clear entry point for users to choose what they need
- Prevents the two components from overwriting each other

## Manual Deployment

If needed, you can trigger a manual deployment by:

1. Going to the Actions tab in GitHub
2. Selecting the "Unified GitHub Pages Deployment" workflow
3. Clicking "Run workflow" and selecting the main branch

## Important Files

- `/.github/workflows/unified-pages-deploy.yml`: The deployment workflow
- `/landing-page.html`: The source for the landing page (copied to index.html during deployment)
- `/mkdocs/mkdocs.yml`: MkDocs configuration for the documentation site
- `/dashboard.html`: Main dashboard file (and related JavaScript files)