# Awesome-Virome Dashboard

This interactive dashboard provides a visual exploration of the various tools in the Awesome-Virome repository.

## Features

- Filter tools by category, language, update recency, and star count
- Interactive visualizations of tool distribution across categories
- Maintenance status overview showing recently updated vs. stale tools
- Language usage breakdown
- Stars distribution histogram
- Update timeline showing when tools were last updated
- Top tools ranking by popularity
- Detailed tool information panel
- Complete searchable and sortable data table

## Running the Dashboard

### Option 1: Using a Simple Server (Recommended)

To avoid CORS issues when loading the data file, it's best to run the dashboard using a local server.

1. Make sure you have Node.js installed
2. Run the included server script:
   ```
   node serve.js
   ```
3. Open your browser and navigate to:
   ```
   http://localhost:8080
   ```

### Option 2: Direct File Opening

You can also open the `index.html` file directly in your browser, but you may encounter CORS issues when loading the data file. If you see an error loading data, click the "Use Sample Data Instead" button to use generated sample data.

## Troubleshooting

### Error Loading Data

If you see "Error loading data" when using the dashboard, it could be due to:

1. The `data.json` file is missing or in the wrong location
2. There's a syntax error in the JSON file
3. You're opening the file directly in a browser (CORS policy prevents local file access)

Solutions:
- Use the "Use Sample Data Instead" button to view the dashboard with generated data
- Run the included server script to serve the files properly
- Check that `data.json` exists and is valid JSON

### Generating Data from README.md

If you want to generate real data from your README.md file:

1. Make sure you have Node.js installed
2. Run:
   ```
   node generate_data.js > data.json
   ```
3. This will parse the README.md file and create a properly formatted data.json file

## Browser Compatibility

The dashboard works best in modern browsers:
- Chrome 80+
- Firefox 75+
- Safari 13+
- Edge 80+

## Technical Details

This dashboard uses:
- D3.js v7 for visualizations
- Bootstrap 5 for layout and styling
- JavaScript ES6 features for data processing and interaction

The code is structured into modular functions that handle specific aspects of the dashboard:
- Data loading and processing
- Filter application
- Chart creation and updates
- Tool details display