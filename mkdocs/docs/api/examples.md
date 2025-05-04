# API Usage Examples

This page provides practical examples of how to use the Awesome-Virome API in different programming languages and for various use cases.

## Basic Examples

### JavaScript (Browser)

Fetch all tools and display them in a table:

```javascript
async function fetchAndDisplayTools() {
  try {
    const response = await fetch('https://shandley.github.io/awesome-virome/api/v1/tools.json');
    const data = await response.json();
    
    const toolsTable = document.getElementById('tools-table');
    const tbody = toolsTable.querySelector('tbody');
    
    data.tools.forEach(tool => {
      const row = document.createElement('tr');
      row.innerHTML = `
        <td><a href="${tool.url}" target="_blank">${tool.name}</a></td>
        <td>${tool.description}</td>
        <td>${tool.category}</td>
        <td>${tool.language || 'N/A'}</td>
        <td>${tool.maintenance_status || 'Unknown'}</td>
      `;
      tbody.appendChild(row);
    });
  } catch (error) {
    console.error('Error fetching tools:', error);
  }
}

// Call the function when the page loads
document.addEventListener('DOMContentLoaded', fetchAndDisplayTools);
```

### Python

Find tools by category and save to CSV:

```python
import requests
import pandas as pd

def get_tools_by_category(category_slug):
    url = f'https://shandley.github.io/awesome-virome/api/v1/categories/{category_slug}.json'
    response = requests.get(url)
    
    if response.status_code == 200:
        data = response.json()
        return data
    else:
        print(f"Error fetching category: {response.status_code}")
        return None

# Get all host prediction tools
host_prediction_data = get_tools_by_category('host-prediction')

if host_prediction_data:
    # Convert to DataFrame
    tools_df = pd.DataFrame(host_prediction_data['tools'])
    
    # Save to CSV
    tools_df.to_csv('host_prediction_tools.csv', index=False)
    
    print(f"Saved {len(tools_df)} host prediction tools to CSV")
    
    # Basic analysis
    if 'language' in tools_df.columns:
        print("\nTools by Language:")
        print(tools_df['language'].value_counts())
    
    if 'maintenance_status' in tools_df.columns:
        print("\nTools by Maintenance Status:")
        print(tools_df['maintenance_status'].value_counts())
```

### R

Create a visualization of tools by category:

```r
library(httr)
library(jsonlite)
library(ggplot2)
library(dplyr)

# Fetch statistics
stats_url <- "https://shandley.github.io/awesome-virome/api/v1/stats.json"
stats_response <- GET(stats_url)
stats_data <- fromJSON(content(stats_response, "text", encoding = "UTF-8"))

# Extract category data for visualization
categories_df <- data.frame(
  Category = names(stats_data$tools_by_category),
  Count = unlist(stats_data$tools_by_category),
  stringsAsFactors = FALSE
)

# Sort by count descending
categories_df <- categories_df %>% arrange(desc(Count))

# Create bar chart
ggplot(categories_df, aes(x = reorder(Category, -Count), y = Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = "Virome Analysis Tools by Category",
       x = "Category",
       y = "Number of Tools") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the plot
ggsave("tools_by_category.png", width = 10, height = 6)
```

## Advanced Examples

### Tool Recommendation System

This JavaScript example creates a simple tool recommendation system based on user requirements:

```javascript
async function recommendTools(requirements) {
  try {
    // Fetch all tools
    const response = await fetch('https://shandley.github.io/awesome-virome/api/v1/tools.json');
    const data = await response.json();
    
    // Score each tool based on requirements
    const scoredTools = data.tools.map(tool => {
      let score = 0;
      
      // Category match
      if (requirements.category && tool.category === requirements.category) {
        score += 3;
      }
      
      // Language match
      if (requirements.language && tool.language === requirements.language) {
        score += 2;
      }
      
      // Maintenance status (prefer active tools)
      if (tool.maintenance_status === 'Active') {
        score += 1;
      }
      
      // Citation count (normalized 0-2 points)
      if (tool.citation_count) {
        const normalizedCitations = Math.min(tool.citation_count / 500, 1) * 2;
        score += normalizedCitations;
      }
      
      // GitHub stars (normalized 0-2 points)
      if (tool.github_stars) {
        const normalizedStars = Math.min(tool.github_stars / 1000, 1) * 2;
        score += normalizedStars;
      }
      
      return {
        ...tool,
        score
      };
    });
    
    // Return top 5 recommendations
    const recommendations = scoredTools
      .sort((a, b) => b.score - a.score)
      .slice(0, 5);
    
    return recommendations;
  } catch (error) {
    console.error('Error fetching tools for recommendation:', error);
    return [];
  }
}

// Example usage
const userRequirements = {
  category: 'Host Prediction',
  language: 'Python'
};

recommendTools(userRequirements)
  .then(recommendations => {
    console.log('Top recommendations:');
    recommendations.forEach((tool, index) => {
      console.log(`${index + 1}. ${tool.name} (Score: ${tool.score.toFixed(1)})`);
      console.log(`   Description: ${tool.description}`);
      console.log(`   URL: ${tool.url}`);
      console.log('');
    });
  });
```

### Python Data Visualization

This Python script creates simple visualizations for tool data:

```python
import requests
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Fetch all tools
def fetch_tools():
    url = 'https://shandley.github.io/awesome-virome/api/v1/tools.json'
    response = requests.get(url)
    if response.status_code == 200:
        return response.json()['tools']
    else:
        print(f"Error fetching tools: {response.status_code}")
        return []

tools = fetch_tools()
df = pd.DataFrame(tools)

# Set up the plotting environment
plt.style.use('ggplot')
plt.figure(figsize=(12, 8))

# Create a bar chart of tools by category
category_counts = df['category'].value_counts()
plt.subplot(2, 1, 1)
sns.barplot(x=category_counts.index, y=category_counts.values)
plt.title('Tools by Category')
plt.xticks(rotation=45, ha='right')
plt.ylabel('Number of Tools')
plt.tight_layout()

# Create a pie chart of tool programming languages
plt.subplot(2, 1, 2)
language_counts = df['language'].value_counts()
plt.pie(language_counts, labels=language_counts.index, autopct='%1.1f%%')
plt.title('Tools by Programming Language')
plt.axis('equal')
plt.tight_layout()

# Save the figure
plt.savefig('awesome_virome_stats.png')
plt.close()

print(f"Visualization saved as awesome_virome_stats.png")

# Print some basic statistics
print("\nBasic Statistics:")
print(f"Total number of tools: {len(df)}")

if 'maintenance_status' in df.columns:
    print("\nTools by Maintenance Status:")
    print(df['maintenance_status'].value_counts())

if 'category' in df.columns:
    print("\nTop Categories:")
    print(category_counts.head(5))

if 'citation_count' in df.columns:
    print("\nMost Cited Tools:")
    most_cited = df.sort_values('citation_count', ascending=False).head(5)
    for _, tool in most_cited.iterrows():
        print(f"- {tool['name']}: {tool['citation_count']} citations")
```

## Command-Line Tool Examples

### Shell Script for Category Summary

This Bash script summarizes tools by category:

```bash
#!/bin/bash

# Function to fetch and display category information
fetch_category() {
  local category_slug=$1
  echo "Fetching information for category: $category_slug"
  
  # Fetch category data
  response=$(curl -s "https://shandley.github.io/awesome-virome/api/v1/categories/$category_slug.json")
  
  # Extract and display information
  category=$(echo $response | jq -r '.category')
  description=$(echo $response | jq -r '.description')
  count=$(echo $response | jq -r '.count')
  
  echo "==================================="
  echo "Category: $category"
  echo "Description: $description"
  echo "Number of tools: $count"
  echo "==================================="
  
  # Show top 5 tools
  echo "Top tools in this category:"
  echo $response | jq -r '.tools[0:5] | .[] | "- \(.name): \(.description)"'
  echo ""
}

# Fetch categories
echo "Fetching available categories..."
categories=$(curl -s "https://shandley.github.io/awesome-virome/api/v1/categories.json")

# Display categories
echo "Available categories:"
echo $categories | jq -r '.categories[] | .slug'
echo ""

# Process each category
echo $categories | jq -r '.categories[] | .slug' | while read -r category_slug; do
  fetch_category "$category_slug"
done
```

## Next Steps

- Explore the [API Reference](endpoints.md) for detailed endpoint documentation
- Check out the [API Overview](overview.md) for general information
- Visit the [Tools](../tools/overview.md) section to learn about the available tools