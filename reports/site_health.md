# Site Health Monitoring

This document provides information about the site health monitoring process for the Awesome Virome repository.

## Health Check Process

The site health is monitored daily through an automated workflow that performs the following checks:

1. **Site Accessibility**: Verifies that the GitHub Pages site is accessible
2. **Data File Availability**: Ensures that the data.json file is available
3. **Data Structure Validation**: Confirms that the data.json file has the expected format
4. **Content Verification**: Checks for essential content on the site

## Health Check Results

The results of the latest health check can be viewed in the [GitHub Actions](https://github.com/shandley/awesome-virome/actions/workflows/site-health-check.yml) tab.

## Dashboard Status Badge

The site health status is displayed as a badge in the repository README:

[![Site Health](https://github.com/shandley/awesome-virome/actions/workflows/site-health-check.yml/badge.svg)](https://github.com/shandley/awesome-virome/actions/workflows/site-health-check.yml)

This badge reflects the current health status of the site based on the latest automated check.

## Troubleshooting

If the health check fails, examine the workflow logs for specific errors. Common issues include:

- GitHub Pages deployment failures
- Malformed data.json file
- Missing content due to repository structure changes

For persistent issues, check the GitHub Pages settings in the repository options or review recent commits that may have affected the site structure.