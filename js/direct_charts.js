/**
 * Direct Charts - Standalone visualization script
 * 
 * This file provides a completely independent implementation
 * of the dashboard visualization that directly manipulates the DOM 
 * without relying on the existing framework.
 */

// Execute when DOM is fully loaded
document.addEventListener('DOMContentLoaded', function() {
    console.log('Direct Charts: Initializing...');
    
    // Check if the dashboard has already loaded the impact data
    if (window.directChartsLoaded) {
        console.log('Direct Charts: Data already loaded by dashboard, skipping fetch');
        return;
    }
    
    if (window.impactData) {
        console.log('Direct Charts: Using pre-loaded impact data from window.impactData');
        renderAllCharts(window.impactData);
        return;
    }
    
    // If not already loaded, load impact data
    console.log('Direct Charts: Loading impact data directly');
    fetch('impact_data.json')
        .then(response => response.json())
        .then(data => {
            console.log('Direct Charts: Successfully loaded impact_data.json');
            window.impactData = data;
            renderAllCharts(data);
        })
        .catch(error => {
            console.error('Direct Charts: Error loading impact data:', error);
        });
});

// Render all charts
function renderAllCharts(data) {
    console.log('Direct Charts: Starting chart rendering');
    
    // Check if we have citation data
    if (data.citations && data.citations.by_year && Object.keys(data.citations.by_year).length > 0) {
        console.log('Direct Charts: Citation data found, years:', Object.keys(data.citations.by_year).length);
        
        // Replace the entire citation trends chart container
        const citationTrendsContainer = document.getElementById('citationTrendsChart');
        if (citationTrendsContainer) {
            renderCitationTrendsChart(citationTrendsContainer, data);
        } else {
            console.error('Direct Charts: Citation trends container not found!');
        }
        
        // Replace the entire top cited tools chart container
        const topCitedContainer = document.getElementById('topCitedToolsChart');
        if (topCitedContainer) {
            renderTopCitedToolsChart(topCitedContainer, data);
        } else {
            console.error('Direct Charts: Top cited tools container not found!');
        }
    } else {
        console.error('Direct Charts: No citation data found in impact_data.json');
    }
    
    // Render network visualization
    const networkContainer = document.getElementById('networkGraph');
    if (networkContainer) {
        renderNetworkGraph(networkContainer, data);
    }
    
    // Render category distribution
    const categoryContainer = document.getElementById('categoryDistributionChart');
    if (categoryContainer) {
        renderCategoryDistributionChart(categoryContainer, data);
    }
    
    // Update stats in the sidebar
    updateStats(data);
}

// Render citation trends chart
function renderCitationTrendsChart(container, data) {
    console.log('Direct Charts: Rendering citation trends chart');
    
    // Make sure Chart.js is available
    if (typeof Chart === 'undefined') {
        container.innerHTML = '<div class="alert alert-warning">Chart.js library not found</div>';
        return;
    }
    
    // Get citation data by year
    const citationsByYear = data.citations.by_year;
    const years = Object.keys(citationsByYear).sort();
    
    // Calculate cumulative citations
    let cumulativeCount = 0;
    const cumulativeData = years.map(year => {
        cumulativeCount += citationsByYear[year];
        return cumulativeCount;
    });
    
    // Create canvas for the chart
    container.innerHTML = '';
    const canvas = document.createElement('canvas');
    canvas.id = 'citationTrendsChartCanvas';
    canvas.style.width = '100%';
    canvas.style.height = '100%';
    container.appendChild(canvas);
    
    // Create the chart
    new Chart(canvas, {
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
                    data: cumulativeData,
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
    
    console.log('Direct Charts: Citation trends chart rendered successfully');
}

// Render top cited tools chart
function renderTopCitedToolsChart(container, data) {
    console.log('Direct Charts: Rendering top cited tools chart');
    
    // Make sure Chart.js is available
    if (typeof Chart === 'undefined') {
        container.innerHTML = '<div class="alert alert-warning">Chart.js library not found</div>';
        return;
    }
    
    // Calculate total citations for each tool
    const tools = data.tools || [];
    const toolsWithCitations = tools.map(tool => {
        const totalCitations = Object.values(tool.citations_by_year || {}).reduce((sum, count) => sum + count, 0);
        return {
            name: tool.name,
            citations: totalCitations
        };
    });
    
    // Sort by citation count and take top 10
    const topTools = toolsWithCitations
        .sort((a, b) => b.citations - a.citations)
        .slice(0, 10);
    
    // Create canvas for the chart
    container.innerHTML = '';
    const canvas = document.createElement('canvas');
    canvas.id = 'topCitedToolsChartCanvas';
    canvas.style.width = '100%';
    canvas.style.height = '100%';
    container.appendChild(canvas);
    
    // Generate colors
    const colors = topTools.map((_, i) => {
        const hue = (i * 36) % 360;
        return `hsla(${hue}, 80%, 60%, 0.8)`;
    });
    
    // Create the chart
    new Chart(canvas, {
        type: 'bar',
        data: {
            labels: topTools.map(tool => tool.name),
            datasets: [{
                label: 'Citations',
                data: topTools.map(tool => tool.citations),
                backgroundColor: colors,
                borderWidth: 1
            }]
        },
        options: {
            indexAxis: 'y',
            responsive: true,
            maintainAspectRatio: false,
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
    
    console.log('Direct Charts: Top cited tools chart rendered successfully');
}

// Render category distribution chart
function renderCategoryDistributionChart(container, data) {
    console.log('Direct Charts: Rendering category distribution chart');
    
    // Make sure Chart.js is available
    if (typeof Chart === 'undefined') {
        container.innerHTML = '<div class="alert alert-warning">Chart.js library not found</div>';
        return;
    }
    
    // Prepare data
    const categories = data.categories || {};
    const categoryNames = Object.keys(categories);
    const categoryCounts = categoryNames.map(name => categories[name]);
    
    // Generate colors
    const colors = categoryNames.map((_, i) => {
        const hue = (i * 25) % 360;
        return `hsla(${hue}, 70%, 60%, 0.8)`;
    });
    
    // Create canvas for the chart
    container.innerHTML = '';
    const canvas = document.createElement('canvas');
    canvas.id = 'categoryDistributionChartCanvas';
    canvas.style.width = '100%';
    canvas.style.height = '100%';
    container.appendChild(canvas);
    
    // Create the chart
    new Chart(canvas, {
        type: 'pie',
        data: {
            labels: categoryNames,
            datasets: [{
                data: categoryCounts,
                backgroundColor: colors,
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
    
    console.log('Direct Charts: Category distribution chart rendered successfully');
}

// Render network graph
function renderNetworkGraph(container, data) {
    console.log('Direct Charts: Rendering network graph');
    
    // Check if vis-network.js is available
    if (typeof vis === 'undefined') {
        container.innerHTML = '<div class="alert alert-warning">vis-network.js library not found</div>';
        return;
    }
    
    // Prepare nodes and edges
    const nodes = data.nodes.map(node => ({
        id: node.id,
        label: node.name,
        title: node.description || node.name,
        group: node.type,
        value: node.type === 'tool' ? (node.citation_count || 5) : (node.size || 10),
        // Different colors for different node types
        color: node.type === 'category' ? '#91CC75' : 
               node.type === 'subcategory' ? '#FAC858' : '#5470C6'
    }));
    
    const edges = data.links.map(link => ({
        from: link.source,
        to: link.target,
        width: link.value || 1
    }));
    
    // Create network
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
        }
    };
    
    try {
        // Create network
        const network = new vis.Network(container, networkData, options);
        
        // Set up zoom controls
        document.getElementById('zoomIn')?.addEventListener('click', () => {
            network.moveTo({
                scale: network.getScale() * 1.2
            });
        });
        
        document.getElementById('zoomOut')?.addEventListener('click', () => {
            network.moveTo({
                scale: network.getScale() / 1.2
            });
        });
        
        document.getElementById('resetNetwork')?.addEventListener('click', () => {
            network.fit({animation: true});
        });
        
        console.log('Direct Charts: Network graph rendered successfully');
    } catch (error) {
        console.error('Direct Charts: Error rendering network graph:', error);
        container.innerHTML = `<div class="alert alert-danger">Error rendering network: ${error.message}</div>`;
    }
}

// Update dashboard statistics
function updateStats(data) {
    console.log('Direct Charts: Updating dashboard statistics');
    
    // Total citations
    const totalCitationsElement = document.getElementById('totalCitations');
    if (totalCitationsElement && data.citations) {
        totalCitationsElement.textContent = data.citations.total.toLocaleString();
    }
    
    // Total tools
    const totalToolsElement = document.getElementById('totalTools');
    if (totalToolsElement && data.nodes) {
        const toolCount = data.nodes.filter(node => node.type === 'tool').length;
        totalToolsElement.textContent = toolCount.toLocaleString();
    }
    
    // Average citations
    const avgCitationsElement = document.getElementById('avgCitations');
    if (avgCitationsElement && data.citations && data.nodes) {
        const toolCount = data.nodes.filter(node => node.type === 'tool').length;
        const avgCitations = toolCount > 0 ? data.citations.total / toolCount : 0;
        avgCitationsElement.textContent = avgCitations.toFixed(1);
    }
}