# API Endpoints

This page details all available endpoints in the Awesome-Virome API, including their parameters, response format, and example usage.

## Tools

Returns a complete list of all tools in the database with their full metadata.

- **URL**: `/api/v1/tools.json`
- **Method**: GET
- **Parameters**: None
- **Response Format**: JSON

**Response Structure**:
```json
{
  "count": 300,
  "last_updated": "2025-05-01",
  "tools": [
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
    },
    // Additional tools...
  ]
}
```

**Example Request**:
```javascript
fetch('https://shandley.github.io/awesome-virome/api/v1/tools.json')
  .then(response => response.json())
  .then(data => console.log(data));
```

## Categories

Returns a list of all tool categories and their subcategories.

- **URL**: `/api/v1/categories.json`
- **Method**: GET
- **Parameters**: None
- **Response Format**: JSON

**Response Structure**:
```json
{
  "count": 8,
  "categories": [
    {
      "name": "Virus and Phage Identification",
      "slug": "virus-identification",
      "description": "Tools for identifying viral sequences in metagenomic data",
      "subcategories": [
        "Metagenome Analysis",
        "Integrated Viruses",
        "RNA Virus Identification"
      ],
      "tool_count": 45
    },
    // Additional categories...
  ]
}
```

**Example Request**:
```python
import requests
response = requests.get('https://shandley.github.io/awesome-virome/api/v1/categories.json')
data = response.json()
print(f"Found {data['count']} categories")
```

## Search

Provides a lightweight index for client-side searching and filtering.

- **URL**: `/api/v1/search.json`
- **Method**: GET
- **Parameters**: None
- **Response Format**: JSON

**Response Structure**:
```json
{
  "count": 300,
  "last_updated": "2025-05-01",
  "tools": [
    {
      "name": "VirSorter2",
      "category": "Virus and Phage Identification",
      "subcategory": "Metagenome Analysis",
      "tags": ["virus detection", "metagenomics", "classifier"],
      "language": "Python"
    },
    // Additional simplified tool entries...
  ]
}
```

**Example Request**:
```javascript
// Client-side filtering example
fetch('https://shandley.github.io/awesome-virome/api/v1/search.json')
  .then(response => response.json())
  .then(data => {
    const pythonTools = data.tools.filter(tool => tool.language === 'Python');
    console.log(`Found ${pythonTools.length} Python tools`);
  });
```

## Stats

Provides aggregate statistics about the tool collection.

- **URL**: `/api/v1/stats.json`
- **Method**: GET
- **Parameters**: None
- **Response Format**: JSON

**Response Structure**:
```json
{
  "total_tools": 300,
  "tools_by_category": {
    "Virus and Phage Identification": 45,
    "Host Prediction": 20,
    "Genome Analysis": 35,
    // Additional categories...
  },
  "tools_by_language": {
    "Python": 150,
    "R": 50,
    "C++": 30,
    // Additional languages...
  },
  "tools_by_maintenance": {
    "Active": 200,
    "Inactive": 80,
    "Unknown": 20
  },
  "tools_by_package_manager": {
    "conda": 120,
    "pip": 80,
    "other": 100
  },
  "average_citation_count": 105.3,
  "last_updated": "2025-05-01"
}
```

**Example Request**:
```python
import requests
import matplotlib.pyplot as plt

response = requests.get('https://shandley.github.io/awesome-virome/api/v1/stats.json')
data = response.json()

# Visualize tools by category
categories = list(data['tools_by_category'].keys())
counts = list(data['tools_by_category'].values())

plt.figure(figsize=(10, 6))
plt.bar(categories, counts)
plt.title('Tools by Category')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.show()
```

## Category-Specific

Returns tools filtered by a specific category.

- **URL**: `/api/v1/categories/{category_slug}.json`
- **Method**: GET
- **Path Parameters**:
  - `category_slug`: The slug of the category to filter by (e.g., `virus-identification`)
- **Response Format**: JSON

**Response Structure**:
```json
{
  "category": "Virus and Phage Identification",
  "description": "Tools for identifying viral sequences in metagenomic data",
  "subcategories": [
    "Metagenome Analysis",
    "Integrated Viruses",
    "RNA Virus Identification"
  ],
  "count": 45,
  "tools": [
    {
      "name": "VirSorter2",
      "description": "Random forest classifier for virus detection",
      "url": "https://bitbucket.org/MAVERICLab/virsorter2/",
      "subcategory": "Metagenome Analysis",
      "language": "Python",
      "github_stars": null,
      "package_manager": "conda",
      "latest_version": "v2.2.4",
      "latest_release_date": "2023-04-15",
      "license": "GPL-3.0",
      "citation_count": 342,
      "maintenance_status": "Active"
    },
    // Additional tools in this category...
  ]
}
```

**Example Request**:
```javascript
// Get all host prediction tools
fetch('https://shandley.github.io/awesome-virome/api/v1/categories/host-prediction.json')
  .then(response => response.json())
  .then(data => {
    console.log(`Found ${data.count} host prediction tools`);
    // Display tools in a table
    const tableData = data.tools.map(tool => ({
      name: tool.name,
      description: tool.description,
      language: tool.language || 'Unknown',
      stars: tool.github_stars || 'N/A'
    }));
    console.table(tableData);
  });
```

## Metadata

Returns detailed information about the API itself.

- **URL**: `/api/v1/metadata.json`
- **Method**: GET
- **Parameters**: None
- **Response Format**: JSON

**Response Structure**:
```json
{
  "name": "Awesome-Virome API",
  "version": "v1",
  "documentation_url": "https://shandley.github.io/awesome-virome/docs/api",
  "description": "REST API for accessing the Awesome-Virome database of virome analysis tools",
  "last_updated": "2025-05-01",
  "contact": "https://github.com/shandley/awesome-virome/issues",
  "endpoints": [
    "/api/v1/tools.json",
    "/api/v1/categories.json",
    "/api/v1/search.json",
    "/api/v1/stats.json",
    "/api/v1/categories/{category_slug}.json",
    "/api/v1/metadata.json"
  ],
  "total_tools": 300,
  "repository_url": "https://github.com/shandley/awesome-virome"
}
```

**Example Request**:
```python
import requests

response = requests.get('https://shandley.github.io/awesome-virome/api/v1/metadata.json')
data = response.json()
print(f"API Version: {data['version']}")
print(f"Last Updated: {data['last_updated']}")
print(f"Available Endpoints: {', '.join(data['endpoints'])}")
```

## Error Handling

If an error occurs, the API will return a JSON object with an error message:

```json
{
  "error": "Category not found",
  "status": 404,
  "message": "The category 'invalid-category' does not exist"
}
```

Note that since the API is hosted on GitHub Pages, only HTTP status 200 (OK) and 404 (Not Found) are possible. Error details are included in the response body.

## Cross-Origin Resource Sharing

All API endpoints support CORS (Cross-Origin Resource Sharing), allowing them to be accessed from web applications hosted on different domains.

## Next Steps

- Check out the [API Examples](examples.md) page for more detailed usage examples
- View the [API Overview](overview.md) for general information about the API
- See the [Tools Overview](../tools/overview.md) to understand the data model better