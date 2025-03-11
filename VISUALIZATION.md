# Awesome-Virome Visualization

This is an interactive visualization for the Awesome-Virome repository that displays relationships between virome analysis tools across different categories.

## Features

- Interactive network/graph visualization using D3.js
- Tools represented as nodes, color-coded by update status
- Node size corresponds to GitHub stars (larger = more popular)
- Tool categories and subcategories shown in the graph structure
- Interconnected tools that are commonly used together in workflows
- Comprehensive filtering options:
  - Filter by category
  - Filter by programming language
  - Filter by last update date
  - Filter by minimum star count
  - Search for specific tools
- Detailed information display for selected tools
- Responsive design for both desktop and mobile viewing
- Interactive zooming and panning

## How to Use

1. View the visualization by opening `index.html` in a browser or visiting the GitHub Pages site
2. Use the filters on the left to narrow down the tools shown
3. Hover over nodes to see basic information in a tooltip
4. Click on a tool node to see detailed information
5. Drag nodes to rearrange the visualization
6. Zoom in and out using mouse wheel or touchpad gestures
7. Use the reset button to clear all filters

## Data Generation

The visualization uses data generated from the main README.md file. If you want to update the visualization with the latest data:

1. Make sure you have Node.js installed
2. Run the data generation script:
   ```
   node generate_data.js > data.json
   ```
3. This will parse the README.md and repo_updates.json files and generate an updated data.json

## File Structure

- `index.html` - Main HTML file for the visualization
- `styles.css` - CSS styles for the visualization
- `visualization.js` - D3.js code for creating the interactive visualization
- `data.json` - Data file containing the tool information
- `generate_data.js` - Script to generate the data file from README.md
- `VISUALIZATION.md` - This documentation file

## Implementation Details

The visualization uses the following technologies:

- D3.js (v7) for the force-directed graph visualization
- Bootstrap 5 for responsive layout and UI components
- JavaScript for data processing and interaction handling

Nodes are sized based on GitHub stars and colored based on update time:
- Green: Updated within 6 months
- Yellow: Updated within 1 year
- Red: Not updated for over 1 year
- Gray: No update data available

## Future Improvements

Potential future enhancements for this visualization:

- Add ability to compare tools side by side
- Implement better clustering of related tools
- Add statistics section showing distribution of tools by language, update frequency, etc.
- Improve mobile interaction for touch devices
- Add integration with GitHub API for real-time stars and update data

## Contributing

If you'd like to improve this visualization:

1. Fork the repository
2. Make your changes
3. Submit a pull request

Focus areas for contributions:
- Performance optimizations for handling larger datasets
- More intuitive interactions for filtering and exploration
- Additional visualization options or views
- Improved mobile experience

## License

This visualization is released under the same license as the main Awesome-Virome repository.

---

Return to the [main repository](https://github.com/yourusername/awesome-virome)