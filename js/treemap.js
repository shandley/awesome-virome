/**
 * Tool Categorization Treemap Visualization
 * Implements a hierarchical treemap visualization to display tools grouped by categories and subcategories
 */

class TreemapVisualization {
    constructor() {
        this.data = null;
        this.chartInstance = null;
        this.containerId = 'treemapContainer';
        this.canvasId = 'treemapChart';
    }

    // Initialize treemap with data
    initialize(data) {
        console.log('TreemapVisualization.initialize called with data:', data);
        if (!data) {
            console.error('Treemap initialization failed: No data provided');
            return;
        }
        this.data = data;
        this.buildTreemap();
    }

    // Process data into hierarchical structure for treemap
    processData() {
        if (!this.data || !this.data.nodes) {
            console.error('No data available for treemap visualization');
            return null;
        }

        // Extract categories and subcategories
        const categories = this.data.nodes.filter(node => node.type === 'category');
        const subcategories = this.data.nodes.filter(node => node.type === 'subcategory');
        const tools = this.data.nodes.filter(node => node.type === 'tool');

        // Create hierarchical structure
        const hierarchyData = {
            name: 'Tools',
            children: []
        };

        // Process categories
        categories.forEach(category => {
            const categoryNode = {
                name: category.name,
                id: category.id,
                children: []
            };

            // Find subcategories for this category
            const categorySubcategories = this.findSubcategoriesForCategory(category, subcategories);
            
            // Process tools directly in this category (without subcategory)
            const directTools = this.findToolsForCategory(category, tools, true);
            
            // Add subcategories
            categorySubcategories.forEach(subcategory => {
                const subcategoryTools = this.findToolsForSubcategory(subcategory, tools);
                
                if (subcategoryTools.length > 0) {
                    const subcategoryNode = {
                        name: subcategory.name,
                        id: subcategory.id,
                        children: subcategoryTools.map(tool => ({
                            name: tool.name,
                            id: tool.id,
                            value: this.calculateToolValue(tool),
                            category: category.name,
                            subcategory: subcategory.name,
                            url: tool.url,
                            description: tool.description
                        }))
                    };
                    categoryNode.children.push(subcategoryNode);
                }
            });
            
            // Add direct tools
            if (directTools.length > 0) {
                // Add tools directly under category (without subcategory)
                directTools.forEach(tool => {
                    categoryNode.children.push({
                        name: tool.name,
                        id: tool.id,
                        value: this.calculateToolValue(tool),
                        category: category.name,
                        subcategory: null,
                        url: tool.url,
                        description: tool.description
                    });
                });
            }
            
            // Only add categories that have tools
            if (categoryNode.children.length > 0) {
                hierarchyData.children.push(categoryNode);
            }
        });
        
        // Check for uncategorized tools
        const uncategorizedTools = tools.filter(tool => !tool.category);
        if (uncategorizedTools.length > 0) {
            hierarchyData.children.push({
                name: 'Uncategorized',
                id: 'category-uncategorized',
                children: uncategorizedTools.map(tool => ({
                    name: tool.name,
                    id: tool.id,
                    value: this.calculateToolValue(tool),
                    category: 'Uncategorized',
                    subcategory: null,
                    url: tool.url,
                    description: tool.description
                }))
            });
        }

        return hierarchyData;
    }

    // Find subcategories linked to a category
    findSubcategoriesForCategory(category, subcategories) {
        // Use links to find subcategories connected to this category
        const categoryId = category.id;
        const links = this.data.links || [];
        
        const connectedSubcategoryIds = links
            .filter(link => 
                (link.source === categoryId && this.isSubcategoryId(link.target)) ||
                (link.target === categoryId && this.isSubcategoryId(link.source))
            )
            .map(link => 
                link.source === categoryId ? link.target : link.source
            );
        
        return subcategories.filter(subcategory => 
            connectedSubcategoryIds.includes(subcategory.id)
        );
    }

    // Helper to check if an ID is a subcategory ID
    isSubcategoryId(id) {
        return typeof id === 'string' && id.startsWith('subcategory-');
    }

    // Find tools for a category
    findToolsForCategory(category, tools, directOnly = false) {
        const categoryId = category.id;
        const categoryName = category.name;
        const links = this.data.links || [];
        
        // Find tools directly connected to this category
        const connectedToolIds = links
            .filter(link => 
                (link.source === categoryId && this.isToolId(link.target)) ||
                (link.target === categoryId && this.isToolId(link.source))
            )
            .map(link => 
                link.source === categoryId ? link.target : link.source
            );
        
        // If directOnly is true, return only tools directly connected
        if (directOnly) {
            return tools.filter(tool => 
                connectedToolIds.includes(tool.id) ||
                tool.category === categoryName
            );
        }
        
        // Otherwise, include tools listing this category name
        return tools.filter(tool => 
            connectedToolIds.includes(tool.id) ||
            tool.category === categoryName
        );
    }

    // Find tools for a subcategory
    findToolsForSubcategory(subcategory, tools) {
        const subcategoryId = subcategory.id;
        const subcategoryName = subcategory.name;
        const links = this.data.links || [];
        
        // Find tools connected to this subcategory
        const connectedToolIds = links
            .filter(link => 
                (link.source === subcategoryId && this.isToolId(link.target)) ||
                (link.target === subcategoryId && this.isToolId(link.source))
            )
            .map(link => 
                link.source === subcategoryId ? link.target : link.source
            );
        
        return tools.filter(tool => 
            connectedToolIds.includes(tool.id) ||
            tool.subcategory === subcategoryName
        );
    }

    // Helper to check if an ID is a tool ID
    isToolId(id) {
        return typeof id === 'string' && !id.startsWith('category-') && !id.startsWith('subcategory-');
    }

    // Calculate a size value for the tool based on stars, citations, etc.
    calculateToolValue(tool) {
        // Base value - all tools have a minimum size
        let value = 20;
        
        // Add value based on stars
        if (tool.stars) {
            value += Math.min(tool.stars / 10, 80);
        }
        
        // Add value based on citation count if available
        if (tool.citation_count) {
            value += Math.min(tool.citation_count / 5, 60);
        } else if (tool.citations_by_year) {
            const totalCitations = Object.values(tool.citations_by_year)
                .reduce((sum, count) => sum + count, 0);
            value += Math.min(totalCitations / 5, 60);
        }
        
        return value;
    }

    // Build the treemap visualization
    buildTreemap() {
        console.log('Building treemap visualization...');
        
        // Process data for treemap
        const treeData = this.processData();
        if (!treeData) {
            console.error('Failed to process data for treemap');
            this.showError('Failed to process data for treemap visualization');
            return;
        }
        
        console.log('Treemap data processed successfully:', treeData);
        
        // Create container if it doesn't exist
        let container = document.getElementById(this.containerId);
        if (!container) {
            console.error(`Treemap container #${this.containerId} not found`);
            return;
        }
        
        console.log('Found treemap container:', container);
        
        // Clear previous chart if it exists
        if (this.chartInstance) {
            this.chartInstance.destroy();
        }
        
        // Create canvas if it doesn't exist
        let canvas = document.getElementById(this.canvasId);
        if (!canvas) {
            console.log('Canvas element not found, creating a new one...');
            // Clear container and create canvas
            container.innerHTML = '';
            canvas = document.createElement('canvas');
            canvas.id = this.canvasId;
            canvas.style.width = '100%';
            canvas.style.height = '100%';
            container.appendChild(canvas);
            console.log('Created new canvas element:', canvas);
        } else {
            console.log('Found existing canvas element:', canvas);
        }
        
        // Get the chart context
        const ctx = canvas.getContext('2d');
        if (!ctx) {
            console.error('Failed to get 2D context from canvas');
            this.showError('Failed to initialize canvas for chart rendering');
            return;
        }
        console.log('Canvas 2D context obtained successfully');
        
        // Check if Chart.js is available
        if (typeof Chart === 'undefined') {
            console.error('Chart.js library is required for treemap visualization');
            container.innerHTML = `
                <div class="alert alert-warning">
                    Chart.js library is required for treemap visualization.
                </div>
            `;
            return;
        }
        
        // Generate colors for categories
        const categoryColors = this.generateCategoryColors(treeData.children);
        
        // Check if treemap plugin is loaded
        if (!Chart.registry.controllers.get('treemap')) {
            console.error('Chart.js treemap plugin is not registered');
            this.showError('Chart.js treemap plugin is not available. Please ensure the treemap plugin is properly loaded.');
            return;
        }
        
        console.log('Chart.js and treemap plugin are ready for chart creation');
        
        // Create Treemap using Chart.js
        try {
            console.log('Creating treemap chart instance...');
            this.chartInstance = new Chart(ctx, {
                type: 'treemap',
                data: {
                    datasets: [{
                        tree: treeData,
                        key: 'value',
                        groups: ['name'],
                        spacing: 2,
                        borderWidth: 1,
                        borderColor: 'rgba(0, 0, 0, 0.2)',
                        backgroundColor: (ctx) => {
                            const item = ctx.raw;
                            if (!item || !item._data) return 'rgba(0, 0, 0, 0.1)';
                            
                            // Get category name
                            const categoryName = item._data.category || item._data.name;
                            
                            // Return appropriate color
                            if (categoryColors[categoryName]) {
                                // Leaf nodes (tools) get a lighter version of the category color
                                if (item._data.value) {
                                    return this.lightenColor(categoryColors[categoryName], 0.2);
                                }
                                return categoryColors[categoryName];
                            }
                            
                            // Default color for anything missing
                            return 'rgba(100, 100, 100, 0.5)';
                        },
                        labels: {
                            display: true,
                            align: 'center',
                            position: 'middle',
                            color: (ctx) => {
                                // For smaller rectangles, use white text for readability
                                if (ctx.raw && ctx.raw.h < 30) {
                                    return 'rgba(255, 255, 255, 0.8)';
                                }
                                return 'rgba(0, 0, 0, 0.7)';
                            },
                            font: {
                                size: (ctx) => {
                                    // Adjust font size based on rectangle height
                                    if (ctx.raw) {
                                        if (ctx.raw.h < 25) return 10;
                                        if (ctx.raw.h < 50) return 12;
                                        return 14;
                                    }
                                    return 12;
                                }
                            },
                            formatter: (ctx) => {
                                if (ctx.raw._data.name) {
                                    return ctx.raw._data.name;
                                }
                                return '';
                            }
                        }
                    }]
                },
                options: {
                    responsive: true,
                    maintainAspectRatio: false,
                    plugins: {
                        title: {
                            display: true,
                            text: 'Tool Categorization by Function',
                            font: {
                                size: 16,
                                weight: 'bold'
                            },
                            padding: {
                                top: 10,
                                bottom: 20
                            }
                        },
                        legend: {
                            display: false
                        },
                        tooltip: {
                            callbacks: {
                                title: (tooltipItems) => {
                                    const item = tooltipItems[0];
                                    if (item && item.raw && item.raw._data) {
                                        return item.raw._data.name;
                                    }
                                    return '';
                                },
                                label: (tooltipItem) => {
                                    const item = tooltipItem.raw._data;
                                    
                                    // For leaf nodes (tools)
                                    if (item.value) {
                                        let label = [];
                                        if (item.category) {
                                            label.push(`Category: ${item.category}`);
                                        }
                                        if (item.subcategory) {
                                            label.push(`Subcategory: ${item.subcategory}`);
                                        }
                                        if (item.description) {
                                            // Truncate description if too long
                                            const desc = item.description.length > 100 
                                                ? item.description.substring(0, 100) + '...' 
                                                : item.description;
                                            label.push(`Description: ${desc}`);
                                        }
                                        return label;
                                    }
                                    
                                    // For group nodes, show count of children
                                    if (item.children) {
                                        return `Contains ${item.children.length} items`;
                                    }
                                    
                                    return '';
                                }
                            }
                        }
                    },
                    onClick: (event, elements) => {
                        if (elements && elements.length > 0) {
                            const clickedItem = elements[0];
                            if (clickedItem.raw && clickedItem.raw._data) {
                                const data = clickedItem.raw._data;
                                
                                // Handle click on tool
                                if (data.url) {
                                    // Find the tool in the original data
                                    const toolId = data.id;
                                    const tool = this.data.nodes.find(node => node.id === toolId);
                                    
                                    // Show tool details if available
                                    if (tool && typeof window.showToolDetail === 'function') {
                                        window.showToolDetail(tool);
                                    }
                                }
                            }
                        }
                    }
                }
            });
        } catch (error) {
            console.error('Error creating treemap chart:', error);
            console.log('Error details:', {
                errorType: error.name,
                message: error.message,
                stack: error.stack,
                chartContext: ctx ? 'valid' : 'invalid',
                containerStatus: container ? 'found' : 'missing',
                chartJsStatus: typeof Chart !== 'undefined' ? 'loaded' : 'missing',
                treemapPluginStatus: (typeof Chart !== 'undefined' && Chart.registry && Chart.registry.controllers) 
                    ? (Chart.registry.controllers.get('treemap') ? 'loaded' : 'missing') 
                    : 'unknown'
            });
            
            container.innerHTML = `
                <div class="alert alert-danger">
                    <strong>Error creating treemap visualization:</strong> ${error.message}
                    <br>
                    <small>Make sure Chart.js and the Treemap plugin are properly loaded.</small>
                    <br>
                    <small>Check browser console for detailed error information.</small>
                </div>
            `;
        }
    }

    // Generate consistent colors for categories
    generateCategoryColors(categories) {
        const colors = {};
        const rootColorMap = {
            'Virus and Phage Identification': 'var(--color-identification, #4682B4)',
            'Host Prediction': 'var(--color-host, #6A5ACD)',
            'Genome Analysis': 'var(--color-genome, #20B2AA)',
            'Taxonomy': 'var(--color-taxonomy, #DB7093)',
            'Functional Analysis': 'var(--color-functional, #9370DB)',
            'CRISPR Analysis': 'var(--color-crispr, #FF6347)',
            'Databases': 'var(--color-database, #3CB371)',
            'Sequence Analysis': 'var(--color-sequence, #1E90FF)',
            'Uncategorized': 'var(--secondary-color, #2C4B7C)'
        };
        
        // Assign pre-defined colors where available, generate others
        categories.forEach((category, index) => {
            const name = category.name;
            if (rootColorMap[name]) {
                colors[name] = rootColorMap[name];
            } else {
                // Generate colors for categories without predefined colors
                const hue = (index * 137.5) % 360; // Golden angle approximation for good distribution
                colors[name] = `hsl(${hue}, 70%, 50%)`;
            }
        });
        
        return colors;
    }
    
    // Helper function to lighten a color
    lightenColor(color, amount) {
        // Handle CSS variables
        if (color.startsWith('var(')) {
            // Extract the var value and convert to a standard color
            const varName = color.match(/var\((.*?)(?:,|\))/)[1].trim();
            const style = getComputedStyle(document.documentElement);
            color = style.getPropertyValue(varName).trim() || '#6A5ACD'; // Default if var not found
        }
        
        // Convert hex to rgb
        let hex = color;
        if (color.startsWith('#')) {
            hex = color.substring(1);
        } else if (color.startsWith('hsl')) {
            // For HSL, just adjust the lightness
            return color.replace(/(\d+)%\)$/, (match, lightness) => {
                const newLightness = Math.min(parseInt(lightness) + amount * 100, 90);
                return `${newLightness}%)`;
            });
        } else if (color.startsWith('rgb')) {
            // Extract RGB values and convert to hex
            const rgbMatch = color.match(/rgba?\((\d+),\s*(\d+),\s*(\d+)(?:,\s*[\d.]+)?\)/);
            if (rgbMatch) {
                const r = parseInt(rgbMatch[1]);
                const g = parseInt(rgbMatch[2]);
                const b = parseInt(rgbMatch[3]);
                hex = ((1 << 24) + (r << 16) + (g << 8) + b).toString(16).slice(1);
            }
        }
        
        // Convert hex to RGB, lighten, and convert back
        let r = parseInt(hex.substring(0, 2), 16);
        let g = parseInt(hex.substring(2, 4), 16);
        let b = parseInt(hex.substring(4, 6), 16);
        
        r = Math.min(255, r + Math.round(amount * 255));
        g = Math.min(255, g + Math.round(amount * 255));
        b = Math.min(255, b + Math.round(amount * 255));
        
        const newHex = ((1 << 24) + (r << 16) + (g << 8) + b).toString(16).slice(1);
        return `#${newHex}`;
    }
    
    // Helper method to show error messages in treemap container
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
    }
}

// Load treemap plugin for Chart.js if needed
function loadTreemapPlugin() {
    return new Promise((resolve, reject) => {
        console.log('Checking if Chart.js treemap plugin needs to be loaded...');
        
        if (typeof Chart === 'undefined') {
            console.error('Chart.js is not loaded. Please load Chart.js before the treemap plugin.');
            reject(new Error('Chart.js is not loaded'));
            return;
        }
        
        // Check if Chart registry exists
        if (!Chart.registry || !Chart.registry.controllers) {
            console.error('Chart.js registry not found. Your Chart.js version might be incompatible.');
            reject(new Error('Chart.js registry not found'));
            return;
        }
        
        // Check if treemap plugin is already registered
        if (Chart.registry.controllers.get('treemap')) {
            console.log('Treemap plugin is already registered with Chart.js');
            resolve();
            return;
        }
        
        console.log('Treemap plugin not found, attempting to load it dynamically...');
        
        // Create a script element to load the treemap plugin
        const script = document.createElement('script');
        script.src = 'https://cdn.jsdelivr.net/npm/chartjs-chart-treemap@2.3.0/dist/chartjs-chart-treemap.min.js';
        script.onload = () => {
            console.log('Treemap plugin loaded successfully');
            
            // Verify it's been registered correctly
            if (Chart.registry.controllers.get('treemap')) {
                console.log('Treemap plugin registered with Chart.js successfully');
                resolve();
            } else {
                console.error('Treemap plugin loaded but not registered with Chart.js');
                reject(new Error('Treemap plugin loaded but not registered'));
            }
        };
        script.onerror = (err) => {
            console.error('Failed to load treemap plugin:', err);
            reject(new Error('Failed to load treemap plugin'));
        };
        document.head.appendChild(script);
    });
}

// Create treemap visualizer instance
const treemapVisualizer = new TreemapVisualization();

// Initialize when document is loaded and Chart.js is available
document.addEventListener('DOMContentLoaded', function() {
    // Check if Chart.js is loaded
    if (typeof Chart === 'undefined') {
        console.warn('Chart.js not found, loading dynamically...');
        const chartScript = document.createElement('script');
        chartScript.src = 'https://cdn.jsdelivr.net/npm/chart.js@4.4.0/dist/chart.umd.min.js';
        document.head.appendChild(chartScript);
        
        chartScript.onload = loadTreemapAndInitialize;
    } else {
        loadTreemapAndInitialize();
    }

    function loadTreemapAndInitialize() {
        loadTreemapPlugin()
            .then(() => {
                console.log('Treemap plugin loaded successfully');
                
                // Check if data is already loaded (window.visualizationData)
                if (window.visualizationData) {
                    treemapVisualizer.initialize(window.visualizationData);
                } else {
                    // If not loaded, wait for data to be loaded
                    const originalInitializeDashboard = window.initializeDashboard;
                    if (typeof originalInitializeDashboard === 'function') {
                        window.initializeDashboard = function(data) {
                            // Call the original dashboard initialization
                            originalInitializeDashboard(data);
                            
                            // Initialize treemap with the data
                            treemapVisualizer.initialize(data);
                        };
                    }
                }
            })
            .catch(error => {
                console.error('Failed to initialize treemap:', error);
                // Display error message in treemap container
                const container = document.getElementById('treemapContainer');
                if (container) {
                    container.innerHTML = `
                        <div class="alert alert-warning">
                            <i class="bi bi-exclamation-triangle me-2"></i>
                            Failed to load treemap visualization: ${error.message}
                        </div>
                    `;
                }
            });
    }
});

// Expose treemap visualizer globally
window.treemapVisualizer = treemapVisualizer;