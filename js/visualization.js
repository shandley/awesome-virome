/**
 * Enhanced Visualization Dashboard
 * Implements interactive visualizations for tool adoption trends,
 * citation growth, and relationships between tools.
 */

// Main visualization controller
class VisualizationManager {
    constructor() {
        this.data = null;
        this.networkInstance = null;
        this.chartsInitialized = false;
    }

    // Initialize all visualizations with data
    initialize(data) {
        this.data = data;
        this.initializeToolRelationships();
        this.initializeCitationTrends();
        this.initializeTopCitedTools();
        this.initializeAdoptionTimeline();
        this.initializeCategoryDistribution();
        this.initializeLanguageDistribution();
        this.chartsInitialized = true;
    }

    // Force-directed network graph showing tool relationships
    initializeToolRelationships() {
        const container = document.getElementById('networkGraph');
        if (!container) return;

        // Filter and format data for network visualization
        const nodes = this.data.nodes.map(node => ({
            id: node.id,
            label: node.name,
            title: node.description || node.name,
            group: node.type,
            value: node.type === 'tool' ? (node.citation_count || 1) : (node.size || 5),
            color: this.getNodeColor(node)
        }));

        const edges = this.data.links.map(link => ({
            from: link.source,
            to: link.target,
            width: link.value || 1,
            title: link.type || 'related'
        }));

        // Create network with vis.js
        const networkData = { nodes, edges };
        const options = {
            physics: {
                stabilization: true,
                barnesHut: {
                    gravitationalConstant: -10000,
                    springConstant: 0.01,
                    springLength: 200
                }
            },
            nodes: {
                shape: 'dot',
                scaling: {
                    min: 5,
                    max: 30,
                    label: {
                        enabled: true,
                        min: 14,
                        max: 24
                    }
                },
                font: {
                    size: 14
                }
            },
            edges: {
                width: 0.5,
                smooth: {
                    type: 'continuous'
                },
                hoverWidth: 3,
                selectionWidth: 2
            },
            interaction: {
                hover: true,
                tooltipDelay: 300,
                zoomView: true,
                dragView: true
            },
            groups: {
                tool: {
                    color: { background: '#5470C6', border: '#3C5CAA' },
                    shape: 'dot'
                },
                category: {
                    color: { background: '#91CC75', border: '#6E994E' },
                    shape: 'diamond'
                },
                subcategory: {
                    color: { background: '#FAC858', border: '#D1A246' },
                    shape: 'triangle'
                }
            }
        };

        // Create a local reference to vis network if visualizer exists
        if (typeof vis !== 'undefined') {
            this.networkInstance = new vis.Network(container, networkData, options);
            
            // Add event listeners for interactions
            this.networkInstance.on('selectNode', (params) => {
                const nodeId = params.nodes[0];
                const node = this.data.nodes.find(n => n.id === nodeId);
                if (node && node.type === 'tool') {
                    // Find the tool object and show its details
                    if (typeof showToolDetail === 'function') {
                        showToolDetail(node);
                    }
                }
            });
            
            // Set up zoom controls
            document.getElementById('zoomIn')?.addEventListener('click', () => {
                this.networkInstance.moveTo({
                    scale: this.networkInstance.getScale() * 1.2
                });
            });
            
            document.getElementById('zoomOut')?.addEventListener('click', () => {
                this.networkInstance.moveTo({
                    scale: this.networkInstance.getScale() / 1.2
                });
            });
            
            document.getElementById('resetNetwork')?.addEventListener('click', () => {
                this.networkInstance.fit({animation: true});
            });
        } else {
            // Show placeholder if vis.js is not available
            container.innerHTML = `
                <div class="d-flex justify-content-center align-items-center h-100">
                    <div class="text-center">
                        <i class="bi bi-diagram-3 mb-3" style="font-size: 3rem;"></i>
                        <h4>Network Visualization</h4>
                        <p class="text-muted">Include vis-network.js to enable interactive tool relationship visualization.</p>
                    </div>
                </div>
            `;
        }
    }

    // Line chart showing citation trends over time
    initializeCitationTrends() {
        const container = document.getElementById('citationTrendsChart');
        if (!container) return;

        // Use impact_data.json for citation trends instead of extracting from node data
        if (!window.impactData || !window.impactData.citations || !window.impactData.citations.by_year) {
            this.showEmptyChart(container, 'Citation Trends', 'No citation trend data available.');
            return;
        }
        
        const citationsByYear = window.impactData.citations.by_year;
        
        if (Object.keys(citationsByYear).length === 0) {
            this.showEmptyChart(container, 'Citation Trends', 'No citation trend data available.');
            return;
        }

        // Get the years in sequence
        const years = Object.keys(citationsByYear).sort();
        
        // Calculate cumulative citations
        let cumulativeCount = 0;
        const cumulativeData = years.map(year => {
            cumulativeCount += citationsByYear[year];
            return { year, count: cumulativeCount };
        });

        // Create Chart.js visualization
        if (typeof Chart !== 'undefined') {
            const ctx = document.createElement('canvas');
            ctx.width = container.clientWidth;
            ctx.height = container.clientHeight;
            container.innerHTML = '';
            container.appendChild(ctx);
            
            new Chart(ctx, {
                type: 'line',
                data: {
                    labels: years,
                    datasets: [
                        {
                            label: 'Citations per Year',
                            data: years.map(year => citationsByYear[year]),
                            backgroundColor: 'rgba(75, 192, 192, 0.2)',
                            borderColor: 'rgba(75, 192, 192, 1)',
                            borderWidth: 2,
                            pointRadius: 4,
                            tension: 0.1
                        },
                        {
                            label: 'Cumulative Citations',
                            data: cumulativeData.map(d => d.count),
                            backgroundColor: 'rgba(54, 162, 235, 0.2)',
                            borderColor: 'rgba(54, 162, 235, 1)',
                            borderWidth: 2,
                            pointRadius: 4,
                            tension: 0.1
                        }
                    ]
                },
                options: {
                    responsive: true,
                    maintainAspectRatio: false,
                    interaction: {
                        mode: 'index',
                        intersect: false
                    },
                    scales: {
                        x: {
                            display: true,
                            title: {
                                display: true,
                                text: 'Year'
                            }
                        },
                        y: {
                            display: true,
                            title: {
                                display: true,
                                text: 'Citation Count'
                            },
                            beginAtZero: true
                        }
                    },
                    plugins: {
                        title: {
                            display: true,
                            text: 'Citation Growth Over Time',
                            font: {
                                size: 16
                            }
                        },
                        tooltip: {
                            callbacks: {
                                label: function(context) {
                                    const label = context.dataset.label || '';
                                    const value = context.parsed.y;
                                    return `${label}: ${value}`;
                                }
                            }
                        }
                    }
                }
            });
        } else {
            this.showEmptyChart(container, 'Citation Trends', 'Chart.js library required for citation trend visualization.');
        }
    }

    // Bar chart showing most cited tools
    initializeTopCitedTools() {
        const container = document.getElementById('topCitedToolsChart');
        if (!container) return;

        // Use impact_data.json for top cited tools instead of extracting from node data
        if (!window.impactData || !window.impactData.tools || !window.impactData.tools.length) {
            this.showEmptyChart(container, 'Most Cited Tools', 'No citation data available.');
            return;
        }
        
        // Get top tools from impact_data.json
        const toolsWithCitations = window.impactData.tools
            .map(tool => ({
                name: tool.name,
                citation_count: Object.values(tool.citations_by_year || {}).reduce((sum, count) => sum + count, 0)
            }))
            .sort((a, b) => (b.citation_count || 0) - (a.citation_count || 0))
            .slice(0, 10); // Top 10
        
        if (toolsWithCitations.length === 0) {
            this.showEmptyChart(container, 'Most Cited Tools', 'No citation data available.');
            return;
        }

        // Create Chart.js visualization
        if (typeof Chart !== 'undefined') {
            const ctx = document.createElement('canvas');
            ctx.width = container.clientWidth;
            ctx.height = container.clientHeight;
            container.innerHTML = '';
            container.appendChild(ctx);
            
            new Chart(ctx, {
                type: 'bar',
                data: {
                    labels: toolsWithCitations.map(tool => tool.name),
                    datasets: [{
                        label: 'Citations',
                        data: toolsWithCitations.map(tool => tool.citation_count),
                        backgroundColor: toolsWithCitations.map((_, i) => {
                            const hue = (i * 25) % 360;
                            return `hsla(${hue}, 70%, 60%, 0.8)`;
                        }),
                        borderWidth: 1
                    }]
                },
                options: {
                    responsive: true,
                    maintainAspectRatio: false,
                    indexAxis: 'y', // Horizontal bar chart
                    scales: {
                        x: {
                            beginAtZero: true,
                            title: {
                                display: true,
                                text: 'Citations'
                            }
                        }
                    },
                    plugins: {
                        title: {
                            display: true,
                            text: 'Most Cited Tools',
                            font: {
                                size: 16
                            }
                        },
                        legend: {
                            display: false
                        }
                    }
                }
            });
        } else {
            this.showEmptyChart(container, 'Most Cited Tools', 'Chart.js library required for citation visualization.');
        }
    }

    // Timeline showing tool adoption over time
    initializeAdoptionTimeline() {
        const container = document.getElementById('creationTimelineChart');
        if (!container) return;

        // Extract tool creation dates
        const toolsByYear = this.extractToolsByCreationYear();
        
        if (Object.keys(toolsByYear).length === 0) {
            this.showEmptyChart(container, 'Tool Adoption Timeline', 'No tool creation date data available.');
            return;
        }

        // Get years in sequence
        const years = Object.keys(toolsByYear).sort();
        
        // Calculate cumulative tools
        let cumulativeCount = 0;
        const cumulativeData = years.map(year => {
            cumulativeCount += toolsByYear[year];
            return { year, count: cumulativeCount };
        });

        // Create Chart.js visualization
        if (typeof Chart !== 'undefined') {
            const ctx = document.createElement('canvas');
            ctx.width = container.clientWidth;
            ctx.height = container.clientHeight;
            container.innerHTML = '';
            container.appendChild(ctx);
            
            new Chart(ctx, {
                type: 'line',
                data: {
                    labels: years,
                    datasets: [
                        {
                            label: 'New Tools per Year',
                            data: years.map(year => toolsByYear[year]),
                            backgroundColor: 'rgba(153, 102, 255, 0.2)',
                            borderColor: 'rgba(153, 102, 255, 1)',
                            borderWidth: 2,
                            pointRadius: 4,
                            tension: 0.1,
                            type: 'bar'
                        },
                        {
                            label: 'Cumulative Tools',
                            data: cumulativeData.map(d => d.count),
                            backgroundColor: 'rgba(255, 159, 64, 0.2)',
                            borderColor: 'rgba(255, 159, 64, 1)',
                            borderWidth: 2,
                            pointRadius: 4,
                            tension: 0.1,
                            type: 'line'
                        }
                    ]
                },
                options: {
                    responsive: true,
                    maintainAspectRatio: false,
                    scales: {
                        x: {
                            display: true,
                            title: {
                                display: true,
                                text: 'Year'
                            }
                        },
                        y: {
                            display: true,
                            title: {
                                display: true,
                                text: 'Tool Count'
                            },
                            beginAtZero: true
                        }
                    },
                    plugins: {
                        title: {
                            display: true,
                            text: 'Tool Adoption Timeline',
                            font: {
                                size: 16
                            }
                        }
                    }
                }
            });
        } else {
            this.showEmptyChart(container, 'Tool Adoption Timeline', 'Chart.js library required for timeline visualization.');
        }
    }

    // Pie chart showing tools by category
    initializeCategoryDistribution() {
        const container = document.getElementById('categoryDistributionChart');
        if (!container) return;

        // Count tools by category
        const categories = {};
        this.data.nodes
            .filter(node => node.type === 'tool' && node.category)
            .forEach(tool => {
                const category = tool.category;
                categories[category] = (categories[category] || 0) + 1;
            });
        
        // Prepare data for chart
        const categoryNames = Object.keys(categories).sort();
        const categoryCounts = categoryNames.map(cat => categories[cat]);
        
        // Generate colors
        const categoryColors = categoryNames.map((_, i) => {
            const hue = (i * 35) % 360;
            return `hsla(${hue}, 70%, 60%, 0.8)`;
        });

        // Create Chart.js visualization
        if (typeof Chart !== 'undefined') {
            const ctx = document.createElement('canvas');
            ctx.width = container.clientWidth;
            ctx.height = container.clientHeight;
            container.innerHTML = '';
            container.appendChild(ctx);
            
            new Chart(ctx, {
                type: 'pie',
                data: {
                    labels: categoryNames,
                    datasets: [{
                        data: categoryCounts,
                        backgroundColor: categoryColors,
                        borderWidth: 1
                    }]
                },
                options: {
                    responsive: true,
                    maintainAspectRatio: false,
                    plugins: {
                        title: {
                            display: true,
                            text: 'Tools by Category',
                            font: {
                                size: 16
                            }
                        },
                        tooltip: {
                            callbacks: {
                                label: function(context) {
                                    const label = context.label || '';
                                    const value = context.parsed;
                                    const total = context.dataset.data.reduce((a, b) => a + b, 0);
                                    const percentage = Math.round((value / total) * 100);
                                    return `${label}: ${value} (${percentage}%)`;
                                }
                            }
                        }
                    }
                }
            });
        } else {
            this.showEmptyChart(container, 'Category Distribution', 'Chart.js library required for category visualization.');
        }
    }

    // Bar chart showing tools by programming language
    initializeLanguageDistribution() {
        const container = document.getElementById('languageDistributionChart');
        if (!container) return;

        // Count tools by primary language
        const languages = {};
        this.data.nodes
            .filter(node => node.type === 'tool' && node.languages && Object.keys(node.languages).length > 0)
            .forEach(tool => {
                const primaryLanguage = this.getPrimaryLanguage(tool.languages);
                languages[primaryLanguage] = (languages[primaryLanguage] || 0) + 1;
            });
        
        // Sort languages by count (descending)
        const sortedLanguages = Object.entries(languages)
            .sort((a, b) => b[1] - a[1])
            .slice(0, 10); // Top 10
        
        const languageNames = sortedLanguages.map(entry => entry[0]);
        const languageCounts = sortedLanguages.map(entry => entry[1]);
        
        // Generate colors
        const languageColors = languageNames.map((_, i) => {
            const hue = (i * 35) % 360;
            return `hsla(${hue}, 70%, 60%, 0.8)`;
        });

        // Create Chart.js visualization
        if (typeof Chart !== 'undefined') {
            const ctx = document.createElement('canvas');
            ctx.width = container.clientWidth;
            ctx.height = container.clientHeight;
            container.innerHTML = '';
            container.appendChild(ctx);
            
            new Chart(ctx, {
                type: 'bar',
                data: {
                    labels: languageNames,
                    datasets: [{
                        label: 'Number of Tools',
                        data: languageCounts,
                        backgroundColor: languageColors,
                        borderWidth: 1
                    }]
                },
                options: {
                    responsive: true,
                    maintainAspectRatio: false,
                    scales: {
                        y: {
                            beginAtZero: true,
                            title: {
                                display: true,
                                text: 'Number of Tools'
                            }
                        }
                    },
                    plugins: {
                        title: {
                            display: true,
                            text: 'Top Programming Languages',
                            font: {
                                size: 16
                            }
                        },
                        legend: {
                            display: false
                        }
                    }
                }
            });
        } else {
            this.showEmptyChart(container, 'Language Distribution', 'Chart.js library required for language visualization.');
        }
    }

    // Helper function to get primary language from language object
    getPrimaryLanguage(languages) {
        if (!languages || typeof languages !== 'object') return 'Unknown';
        
        // Find language with highest percentage
        let primaryLanguage = 'Unknown';
        let maxPercentage = 0;
        
        for (const [language, percentage] of Object.entries(languages)) {
            if (percentage > maxPercentage) {
                primaryLanguage = language;
                maxPercentage = percentage;
            }
        }
        
        return primaryLanguage;
    }

    // Helper function to extract citation data by year
    extractCitationsByYear() {
        const citationsByYear = {};
        
        this.data.nodes
            .filter(node => node.type === 'tool' && node.citations_by_year)
            .forEach(tool => {
                for (const [year, count] of Object.entries(tool.citations_by_year || {})) {
                    citationsByYear[year] = (citationsByYear[year] || 0) + count;
                }
            });
        
        return citationsByYear;
    }

    // Helper function to extract tool creation dates by year
    extractToolsByCreationYear() {
        const toolsByYear = {};
        
        this.data.nodes
            .filter(node => node.type === 'tool' && node.createdAt)
            .forEach(tool => {
                const year = new Date(tool.createdAt).getFullYear().toString();
                if (year && !isNaN(parseInt(year))) {
                    toolsByYear[year] = (toolsByYear[year] || 0) + 1;
                }
            });
        
        return toolsByYear;
    }

    // Helper function to get color for network nodes
    getNodeColor(node) {
        if (node.type === 'category') {
            return '#91CC75';
        } else if (node.type === 'subcategory') {
            return '#FAC858';
        } else if (node.type === 'tool') {
            // Color tools by category
            if (node.category === 'Virus and Phage Identification') {
                return '#5470C6';
            } else if (node.category === 'Host Prediction') {
                return '#6A5ACD';
            } else if (node.category === 'Genome Analysis') {
                return '#20B2AA';
            } else {
                // Default tool color
                return '#5470C6';
            }
        }
        return '#bbb'; // Default gray
    }

    // Helper function to show placeholder for empty charts
    showEmptyChart(container, title, message) {
        container.innerHTML = `
            <div class="d-flex justify-content-center align-items-center h-100">
                <div class="text-center">
                    <i class="bi bi-bar-chart mb-3" style="font-size: 3rem;"></i>
                    <h4>${title}</h4>
                    <p class="text-muted">${message}</p>
                </div>
            </div>
        `;
    }
}

// Initialize when document is loaded
document.addEventListener('DOMContentLoaded', function() {
    // Create the visualization manager
    window.visualizationManager = new VisualizationManager();
    
    // Modified data loading to initialize visualizations
    const originalInitializeDashboard = window.initializeDashboard;
    
    if (typeof originalInitializeDashboard === 'function') {
        window.initializeDashboard = function(data) {
            // Call the original dashboard initialization
            originalInitializeDashboard(data);
            
            // Load impact data first, then initialize visualizations
            fetch('impact_data.json')
                .then(response => response.json())
                .then(impactData => {
                    window.impactData = impactData;
                    console.log('Impact data loaded successfully');
                    
                    // Initialize enhanced visualizations with both data sources
                    window.visualizationManager.initialize(data);
                    console.log('Enhanced visualizations initialized');
                })
                .catch(error => {
                    console.warn('Could not load impact data:', error);
                    // Initialize visualizations anyway, they'll show placeholders for missing data
                    window.visualizationManager.initialize(data);
                });
        };
    }
});