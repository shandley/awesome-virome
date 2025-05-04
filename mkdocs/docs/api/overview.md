# Awesome-Virome API

Awesome-Virome provides a REST API that allows programmatic access to the curated database of virome analysis tools. This enables developers and researchers to build applications, workflows, or custom analyses on top of the dataset.

## API Features

- **Base URL:** `https://shandley.github.io/awesome-virome/api/v1/`
- **Format:** All endpoints return JSON
- **Authentication:** No authentication required
- **CORS-Enabled:** Accessible from browser applications
- **Cache-Friendly:** Responses include appropriate cache headers

## Getting Started

Here's a simple example of how to fetch the complete list of tools using JavaScript:

```javascript
fetch('https://shandley.github.io/awesome-virome/api/v1/tools.json')
  .then(response => response.json())
  .then(data => {
    console.log(`Found ${data.count} tools`);
    // Process the tools data
  });
```

Or using Python:

```python
import requests

response = requests.get('https://shandley.github.io/awesome-virome/api/v1/tools.json')
data = response.json()
print(f"Found {data['count']} tools")
# Process the tools data
```

## Available Endpoints

| Endpoint | Description | Example |
|----------|-------------|---------|
| [`/api/v1/tools.json`](endpoints.md#tools) | Complete list of all tools with metadata | [View](https://shandley.github.io/awesome-virome/api/v1/tools.json) |
| [`/api/v1/categories.json`](endpoints.md#categories) | List of all tool categories | [View](https://shandley.github.io/awesome-virome/api/v1/categories.json) |
| [`/api/v1/search.json`](endpoints.md#search) | Lightweight index for client-side filtering | [View](https://shandley.github.io/awesome-virome/api/v1/search.json) |
| [`/api/v1/stats.json`](endpoints.md#stats) | Aggregate statistics about the tools collection | [View](https://shandley.github.io/awesome-virome/api/v1/stats.json) |
| [`/api/v1/categories/{category_slug}.json`](endpoints.md#category-specific) | Tools filtered by category | [Example](https://shandley.github.io/awesome-virome/api/v1/categories/virus-identification.json) |

## Data Structure

The API returns data with a consistent structure across endpoints. Here's an example of a tool entry:

```json
{
  "name": "VirSorter2",
  "description": "Random forest classifier for virus detection",
  "url": "https://bitbucket.org/MAVERICLab/virsorter2/",
  "category": "Virus and Phage Identification",
  "subcategory": "Metagenome Analysis",
  "language": "Python",
  "github_stars": null,
  "package_manager": "conda",
  "latest_version": "v2.2.4",
  "latest_release_date": "2023-04-15",
  "license": "GPL-3.0",
  "citation_count": 342,
  "maintenance_status": "Active"
}
```

## Rate Limiting

The API is served from GitHub Pages and doesn't have explicit rate limits. However, we recommend:

- Implementing caching in your applications
- Limiting requests to a reasonable frequency
- Downloading the full dataset for heavy processing rather than making many small requests

## Use Cases

The API enables many interesting use cases:

- **Tool Discovery Applications**: Build custom interfaces for discovering tools
- **Scientific Workflows**: Integrate tool metadata into computational pipelines
- **Academic Research**: Analyze trends in virome analysis tool development
- **Recommendation Systems**: Create intelligent tool recommendation systems based on user needs
- **Custom Dashboards**: Build specialized visualizations for specific types of tools

## Versioning

The API uses a versioned URL structure (currently `v1`). This ensures that future changes to the API won't break existing applications. If we introduce breaking changes, we'll create a new version (e.g., `v2`) while maintaining the old version for backward compatibility.

## API Changelog

### v1 (Current)
- Initial release with core endpoints
- Support for all tool categories
- Basic statistics and search index
- Complete metadata for all tools

## Need Help?

If you have questions about the API or encounter any issues, please:

1. Check the complete [API documentation](endpoints.md)
2. Look at the [example code](examples.md)
3. [Open an issue](https://github.com/shandley/awesome-virome/issues) on GitHub if you find a bug or have a feature request