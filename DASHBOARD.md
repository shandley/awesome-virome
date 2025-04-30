# Enhanced Visualization Dashboard

The Awesome-Virome dashboard provides interactive visualizations to explore relationships between viral analysis tools, track citation growth, and analyze adoption trends.

## Features

### Tool Relationship Network
An interactive force-directed graph showing connections between tools, categories, and subcategories in the virome analysis ecosystem. The visualization helps users understand how different tools relate to each other and identify clusters of related functionality.

### Citation Growth Visualization
Dynamic charts showing:
- Citation trends over time for all tools
- Cumulative citation growth
- Most cited tools in the collection
- Detailed citation metrics for individual tools

### Tool Adoption Timeline
A combined bar and line chart showing:
- New tools created each year
- Cumulative growth of the tool ecosystem
- Adoption trends by category

### Ecosystem Analytics
Interactive charts showing:
- Distribution of tools by category
- Programming language prevalence
- Maintenance status overview
- Key statistics about the collection

## Implementation Details

### Data Sources
The visualization dashboard utilizes two main data sources:

1. **data.json**: Primary data file containing the complete collection of tools, categories, and relationships

2. **impact_data.json**: Aggregated academic impact data including citation trends, tool adoption timelines, and relationship metrics

### Update Process
The impact data can be updated automatically using the provided script:

```bash
python scripts/update_impact_data.py
```

This script processes the main data.json file and generates updated impact metrics for the visualization dashboard.

### Libraries Used
- **D3.js**: For network visualization and custom charts
- **Chart.js**: For standard charts (line, bar, pie)
- **vis-network.js**: For the interactive relationship network
- **DataTables**: For searchable and sortable tool listings

## Future Enhancements

- **Citation Prediction**: Machine learning model to predict future citation trends
- **Comparative Tool Analysis**: Side-by-side comparison of similar tools
- **Research Impact Metrics**: Integration with additional academic impact metrics
- **Community Adoption Indicators**: GitHub stars growth and community engagement metrics
- **Geospatial Usage Data**: World map showing tool usage by region

## Usage Examples

### Finding Related Tools
1. Navigate to the Network section
2. Click on a tool node to see its connections
3. Explore the connected nodes to discover related tools

### Identifying Trending Tools
1. Visit the Citations section
2. View the Citation Trends chart to identify tools with rapid citation growth
3. Select specific years to focus on recent trends

### Analyzing Tool Ecosystem Evolution
1. Go to the Analytics section
2. Examine the Creation Timeline chart to see how the ecosystem has evolved
3. Filter by category to see growth patterns in specific areas