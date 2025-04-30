# Awesome Virome API Documentation

Awesome Virome provides a simple JSON API that allows programmatic access to our curated database of virome analysis tools. This enables developers and researchers to build applications, workflows, or analyses on top of this dataset.

## API Overview

- **Base URL:** `https://shandley.github.io/awesome-virome/api/v1/`
- **Format:** All endpoints return JSON
- **Authentication:** No authentication required
- **Rate Limiting:** Please limit requests to a reasonable number
- **CORS:** Enabled for cross-origin requests

## API Endpoints

### List All Tools

Get a complete list of all virome analysis tools in the database.

- **URL:** `/api/v1/tools.json`
- **Method:** `GET`
- **Response Format:**
  ```json
  {
    "api_version": "1.0",
    "generated_at": "2025-04-30T12:34:56.789Z",
    "count": 250,
    "tools": [
      {
        "id": "VirSorter2",
        "name": "VirSorter2",
        "description": "Random forest classifier for virus detection",
        "url": "https://github.com/jiarong/VirSorter2",
        "category": "Virus and Phage Identification",
        "subcategory": "Metagenome Analysis",
        "language": "Python",
        "stars": 123,
        "forks": 45,
        "license": "GPL-3.0",
        "topics": ["virus", "identification", "metagenomics"],
        "updated_at": "2025-02-15T10:20:30Z",
        "created_at": "2020-01-01T00:00:00Z",
        "languages": {
          "Python": 12345,
          "Shell": 678
        },
        "doi": "10.1093/bioinformatics/btaa1234",
        "citation_count": 42,
        "package_manager": "conda"
      },
      // More tools...
    ]
  }
  ```

### List Categories

Get a list of all tool categories.

- **URL:** `/api/v1/categories.json`
- **Method:** `GET`
- **Response Format:**
  ```json
  {
    "api_version": "1.0",
    "generated_at": "2025-04-30T12:34:56.789Z",
    "count": 15,
    "categories": [
      {
        "slug": "virus-identification",
        "name": "Virus and Phage Identification",
        "count": 45,
        "endpoint": "/api/v1/categories/virus-identification.json"
      },
      {
        "slug": "host-prediction",
        "name": "Host Prediction",
        "count": 20,
        "endpoint": "/api/v1/categories/host-prediction.json"
      },
      // More categories...
    ]
  }
  ```

### Tools By Category

Get tools filtered by a specific category.

- **URL:** `/api/v1/categories/{category_slug}.json`
- **Method:** `GET`
- **URL Parameters:**
  - `category_slug`: The slug of the category (e.g., `virus-identification`, `host-prediction`)
- **Response Format:**
  ```json
  {
    "api_version": "1.0",
    "generated_at": "2025-04-30T12:34:56.789Z",
    "category": "Virus and Phage Identification",
    "count": 45,
    "tools": [
      // Tools in this category...
    ]
  }
  ```

### Search Index

Get a lightweight index optimized for client-side searching and filtering.

- **URL:** `/api/v1/search.json`
- **Method:** `GET`
- **Response Format:**
  ```json
  {
    "api_version": "1.0",
    "generated_at": "2025-04-30T12:34:56.789Z",
    "count": 250,
    "tools": [
      {
        "id": "VirSorter2",
        "name": "VirSorter2",
        "description": "Random forest classifier for virus detection",
        "category": "Virus and Phage Identification",
        "subcategory": "Metagenome Analysis",
        "language": "Python",
        "topics": ["virus", "identification", "metagenomics"],
        "stars": 123,
        "url": "https://github.com/jiarong/VirSorter2"
      },
      // More tools...
    ]
  }
  ```

### Statistics

Get aggregate statistics about the tools database.

- **URL:** `/api/v1/stats.json`
- **Method:** `GET`
- **Response Format:**
  ```json
  {
    "api_version": "1.0",
    "generated_at": "2025-04-30T12:34:56.789Z",
    "total_tools": 250,
    "total_stars": 12500,
    "average_stars": 50.0,
    "languages": {
      "Python": 120,
      "R": 45,
      "C++": 30,
      "Shell": 25,
      "Java": 20,
      "Other": 10
    },
    "package_managers": {
      "conda": 80,
      "pip": 60,
      "docker": 30,
      "npm": 5,
      "other": 15,
      "none": 60
    },
    "categories": {
      "Virus and Phage Identification": 45,
      "Host Prediction": 20,
      "Genome Analysis": 35
      // More categories...
    }
  }
  ```

### API Metadata

Get information about the API itself.

- **URL:** `/api/v1/metadata.json`
- **Method:** `GET`
- **Response Format:**
  ```json
  {
    "api_name": "Awesome Virome API",
    "api_version": "1.0",
    "generated_at": "2025-04-30T12:34:56.789Z",
    "description": "API for accessing the awesome-virome database of virus and phage analysis tools",
    "documentation_url": "https://github.com/shandley/awesome-virome/blob/main/API.md",
    "endpoints": [
      {
        "path": "/api/v1/tools.json",
        "description": "Complete list of tools",
        "methods": ["GET"]
      },
      // More endpoints...
    ],
    "license": "CC0",
    "repository": "https://github.com/shandley/awesome-virome"
  }
  ```

## Using the API

### JavaScript Example

```javascript
// Fetch all tools
fetch('https://shandley.github.io/awesome-virome/api/v1/tools.json')
  .then(response => response.json())
  .then(data => {
    console.log(`Found ${data.count} virome analysis tools`);
    
    // Filter tools with over 50 stars
    const popularTools = data.tools.filter(tool => tool.stars > 50);
    console.log(`${popularTools.length} tools have more than 50 stars`);
    
    // Find tools by category
    const phageTools = data.tools.filter(tool => 
      tool.category === "Virus and Phage Identification"
    );
    console.log(`${phageTools.length} tools are for phage identification`);
  })
  .catch(error => console.error('Error fetching tools:', error));
```

### Python Example

```python
import requests
import pandas as pd

# Fetch tools by category
response = requests.get('https://shandley.github.io/awesome-virome/api/v1/categories/host-prediction.json')
data = response.json()

print(f"Found {data['count']} host prediction tools")

# Convert to pandas DataFrame for analysis
tools_df = pd.DataFrame(data['tools'])

# Analyze by programming language
language_counts = tools_df['language'].value_counts()
print("Programming languages used:")
print(language_counts)

# Find most-cited tools
top_cited = tools_df.sort_values('citation_count', ascending=False).head(5)
print("\nTop 5 most cited tools:")
for _, tool in top_cited.iterrows():
    print(f"{tool['name']}: {tool['citation_count']} citations")
```

### R Example

```r
library(httr)
library(jsonlite)
library(dplyr)

# Get statistics
response <- GET("https://shandley.github.io/awesome-virome/api/v1/stats.json")
stats <- fromJSON(content(response, "text"))

# Plot distribution of tools by category
categories <- data.frame(
  category = names(stats$categories),
  count = unlist(stats$categories)
)

categories <- categories %>%
  arrange(desc(count)) %>%
  mutate(category = factor(category, levels = category))

barplot(
  categories$count, 
  names.arg = categories$category, 
  las = 2, 
  cex.names = 0.7,
  main = "Distribution of Virome Analysis Tools by Category"
)
```

## Notes and Best Practices

1. **Caching**: The API is updated whenever the main dataset changes. Consider caching responses to reduce load on the server.

2. **Error Handling**: Always include proper error handling in your code when making API requests.

3. **Large Responses**: Some endpoints (like `tools.json`) return relatively large responses. Consider using the category-specific endpoints if you only need tools from a specific category.

4. **Client-Side Filtering**: For complex filtering or search, fetch the search index and perform filtering on the client side.

5. **Cross-Origin Requests**: The API supports CORS, so you can use it directly from browser-based applications.

## Data Schema

### Tool Object

| Field | Type | Description |
|-------|------|-------------|
| `id` | String | Unique identifier for the tool |
| `name` | String | Display name of the tool |
| `description` | String | Brief description of the tool |
| `url` | String | URL to the tool's repository or website |
| `category` | String | Primary category |
| `subcategory` | String | Subcategory within the primary category |
| `language` | String | Primary programming language |
| `stars` | Number | GitHub/GitLab stars (if available) |
| `forks` | Number | Number of repository forks (if available) |
| `license` | String | License type (if available) |
| `topics` | Array | List of topic tags |
| `updated_at` | String | ISO 8601 timestamp of last update |
| `created_at` | String | ISO 8601 timestamp of creation date |
| `languages` | Object | Key-value pairs of languages and byte counts |
| `doi` | String | Digital Object Identifier for the publication |
| `citation_count` | Number | Number of citations (if available) |
| `package_manager` | String | Primary package manager (pip, conda, etc.) |

## Feedback and Contributions

If you have suggestions for improving the API or find any issues, please open an issue on the [GitHub repository](https://github.com/shandley/awesome-virome/issues).

To contribute new tools to the database, please follow the [contribution guidelines](https://github.com/shandley/awesome-virome/blob/main/CONTRIBUTING.md).

## License

The API and all data it provides are released under the [Creative Commons Zero v1.0 Universal](https://creativecommons.org/publicdomain/zero/1.0/) license.