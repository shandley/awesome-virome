/**
 * Tool Categorization Treemap Visualization
 * Implements a hierarchical visualization to display tools grouped by categories
 */

class TreemapVisualization {
    constructor() {
        this.data = null;
        this.chartInstance = null;
        this.containerId = 'treemapContainer';
        this.loadingSpinner = null; // Reference to the loading spinner
    }

    // Initialize treemap with data
    initialize(data) {
        console.log('TreemapVisualization.initialize called with data');
        
        if (!data || !data.nodes || !data.links) {
            console.error('Invalid data format provided to treemap visualization');
            this.showError('Invalid data format provided to treemap visualization');
            return;
        }
        
        this.data = data;
        
        // Get container and loading spinner
        const container = document.getElementById(this.containerId);
        if (!container) {
            console.error('Treemap container not found');
            return;
        }
        
        // Find the loading spinner parent
        const chartContainer = container.closest('.chart-container');
        if (chartContainer) {
            this.loadingSpinner = chartContainer.querySelector('.loading-container');
        }
        
        // Create the visualization
        this.createTreemap();
    }
    
    // Show error message in container
    showError(message) {
        const container = document.getElementById(this.containerId);
        if (container) {
            container.innerHTML = `
                <div class="alert alert-danger">
                    <i class="bi bi-exclamation-triangle-fill me-2"></i>
                    ${message}
                </div>
            `;
        }
        
        // Hide the loading spinner if it exists
        if (this.loadingSpinner) {
            this.loadingSpinner.style.display = 'none';
        }
    }
    
    // Create treemap visualization
    createTreemap() {
        try {
            // Process data into hierarchical structure
            const hierarchyData = this.processData();
            if (!hierarchyData) {
                return;
            }
            
            // Get the container
            const container = document.getElementById(this.containerId);
            if (!container) {
                console.error('Treemap container not found');
                return;
            }
            
            // Reset container
            container.innerHTML = '';
            
            // Create table-based treemap (simple fallback if Chart.js doesn't work)
            this.createTableTreemap(container, hierarchyData);
            
            // Hide loading spinner
            if (this.loadingSpinner) {
                this.loadingSpinner.style.display = 'none';
            }
            
        } catch (error) {
            console.error('Error creating treemap:', error);
            this.showError(`Error creating treemap: ${error.message}`);
        }
    }
    
    // Process data into hierarchical structure for treemap
    processData() {
        if (!this.data || !this.data.nodes) {
            this.showError('No data available for treemap visualization');
            return null;
        }
        
        console.log('Processing data for treemap...');
        
        // Extract categories and tools
        const categories = this.data.nodes.filter(node => node.type === 'category');
        const tools = this.data.nodes.filter(node => node.type === 'tool');
        
        // Create a map of tools by category
        const toolsByCategory = {};
        
        // Initialize categories
        categories.forEach(category => {
            toolsByCategory[category.name] = [];
        });
        
        // Find which category each tool belongs to
        tools.forEach(tool => {
            let categoryName = tool.category || 'Uncategorized';
            
            // If categoryName is not found in our map, add it to Uncategorized
            if (!toolsByCategory[categoryName]) {
                if (!toolsByCategory['Uncategorized']) {
                    toolsByCategory['Uncategorized'] = [];
                }
                categoryName = 'Uncategorized';
            }
            
            toolsByCategory[categoryName].push({
                name: tool.name,
                stars: tool.stars || 0,
                url: tool.url || '',
                description: tool.description || ''
            });
        });
        
        // Filter out empty categories and sort by tool count
        const result = Object.entries(toolsByCategory)
            .filter(([_, tools]) => tools.length > 0)
            .sort((a, b) => b[1].length - a[1].length)
            .map(([category, tools]) => ({
                name: category,
                tools: tools.sort((a, b) => b.stars - a.stars)
            }));
            
        console.log(`Processed ${result.length} categories with tools`);
        return result;
    }
    
    // Create a table-based treemap (simple fallback)
    createTableTreemap(container, data) {
        // Create a div for the treemap
        const treemapDiv = document.createElement('div');
        treemapDiv.className = 'table-treemap';
        treemapDiv.style.cssText = 'width: 100%; height: 100%; overflow: auto; padding: 10px;';
        
        // Add title
        const title = document.createElement('h4');
        title.textContent = 'Tool Categorization by Function';
        title.style.cssText = 'text-align: center; margin-bottom: 20px;';
        treemapDiv.appendChild(title);
        
        // Create grid container for categories
        const grid = document.createElement('div');
        grid.style.cssText = 'display: grid; grid-template-columns: repeat(auto-fill, minmax(280px, 1fr)); gap: 15px;';
        treemapDiv.appendChild(grid);
        
        // Color palette for categories
        const colors = [
            '#4682B4', // Steel Blue
            '#20B2AA', // Light Sea Green
            '#6A5ACD', // Slate Blue
            '#DB7093', // Pale Violet Red
            '#3CB371', // Medium Sea Green
            '#FF6347', // Tomato
            '#9370DB', // Medium Purple
            '#1E90FF', // Dodger Blue
            '#FA8072', // Salmon
            '#43CD80'  // Sea Green 3
        ];
        
        // Create category cards
        data.forEach((category, index) => {
            // Pick a color from the palette
            const colorIndex = index % colors.length;
            const color = colors[colorIndex];
            
            // Create category card
            const card = document.createElement('div');
            card.className = 'category-card';
            card.style.cssText = `
                border-radius: 8px;
                border: 1px solid rgba(0,0,0,0.1);
                box-shadow: 0 2px 4px rgba(0,0,0,0.05);
                overflow: hidden;
                height: 100%;
                display: flex;
                flex-direction: column;
            `;
            
            // Category header
            const header = document.createElement('div');
            header.style.cssText = `
                background-color: ${color};
                color: white;
                padding: 10px;
                font-weight: 600;
                font-size: 14px;
                display: flex;
                justify-content: space-between;
                align-items: center;
            `;
            header.innerHTML = `
                <span>${category.name}</span>
                <span class="badge bg-light text-dark">${category.tools.length}</span>
            `;
            card.appendChild(header);
            
            // Tools list
            const toolsList = document.createElement('div');
            toolsList.style.cssText = `
                padding: 10px;
                flex: 1;
                overflow-y: auto;
                max-height: 200px;
                background-color: #f8f9fa;
            `;
            
            // Add top tools (max 5)
            const topTools = category.tools.slice(0, 5);
            topTools.forEach(tool => {
                const toolItem = document.createElement('div');
                toolItem.style.cssText = `
                    padding: 6px 10px;
                    border-bottom: 1px solid rgba(0,0,0,0.05);
                    font-size: 13px;
                    display: flex;
                    justify-content: space-between;
                `;
                toolItem.innerHTML = `
                    <span title="${tool.description || ''}">${tool.name}</span>
                    <span class="stars">
                        <i class="bi bi-star-fill" style="color: #ffc107; font-size: 10px;"></i>
                        ${tool.stars || 0}
                    </span>
                `;
                toolsList.appendChild(toolItem);
            });
            
            // If there are more tools, add a "more" indicator
            if (category.tools.length > 5) {
                const moreItem = document.createElement('div');
                moreItem.style.cssText = `
                    padding: 6px 10px;
                    font-size: 12px;
                    color: #6c757d;
                    text-align: center;
                    font-style: italic;
                `;
                moreItem.textContent = `+ ${category.tools.length - 5} more tools`;
                toolsList.appendChild(moreItem);
            }
            
            card.appendChild(toolsList);
            grid.appendChild(card);
        });
        
        container.appendChild(treemapDiv);
    }
}

// Create and expose treemap visualizer
const treemapVisualizer = new TreemapVisualization();
window.treemapVisualizer = treemapVisualizer;

// Initialize when document is loaded
document.addEventListener('DOMContentLoaded', function() {
    console.log('Treemap script loaded, waiting for data...');
    
    // Check if data is already loaded
    if (window.visualizationData) {
        console.log('Data already available, initializing treemap...');
        treemapVisualizer.initialize(window.visualizationData);
    }
});