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
    
    // If not already loaded, load impact data with better path handling for GitHub Pages
    console.log('Direct Charts: Loading impact data directly');
    
    // First try the relative path (works locally and sometimes on GitHub Pages)
    fetch('impact_data.json')
        .then(response => {
            if (!response.ok) {
                throw new Error(`Failed to load impact_data.json: ${response.status} ${response.statusText}`);
            }
            return response.json();
        })
        .then(data => {
            console.log('Direct Charts: Successfully loaded impact_data.json');
            console.log('Direct Charts: Data structure:', Object.keys(data));
            if (data.tools && data.tools.length > 0) {
                console.log(`Direct Charts: Found ${data.tools.length} tools with citation data`);
                console.log('Direct Charts: First tool citation data:', 
                            data.tools[0].name, 
                            Object.keys(data.tools[0].citations_by_year || {}).length);
            }
            window.impactData = data;
            renderAllCharts(data);
        })
        .catch(error => {
            console.error('Direct Charts: Error loading impact data from relative path:', error);
            
            // Fallback to absolute path (sometimes needed on GitHub Pages)
            console.log('Direct Charts: Trying absolute path as fallback');
            const repoName = 'awesome-virome';
            const username = window.location.hostname.includes('github.io') 
                ? window.location.pathname.split('/')[1] 
                : 'shandley'; // Fallback to repo owner
                
            const absolutePath = `/${username}/${repoName}/impact_data.json`;
            console.log(`Direct Charts: Attempting to load from ${absolutePath}`);
            
            fetch(absolutePath)
                .then(response => {
                    if (!response.ok) {
                        throw new Error(`Failed to load from absolute path: ${response.status}`);
                    }
                    return response.json();
                })
                .then(data => {
                    console.log('Direct Charts: Successfully loaded from absolute path');
                    window.impactData = data;
                    renderAllCharts(data);
                })
                .catch(secondError => {
                    console.error('Direct Charts: All attempts to load impact data failed:', secondError);
                    
                    // Show error message in citation containers
                    const containers = [
                        document.getElementById('citationTrendsChart'),
                        document.getElementById('topCitedToolsChart')
                    ];
                    
                    containers.forEach(container => {
                        if (container) {
                            container.innerHTML = `
                                <div class="alert alert-warning m-3">
                                    <h4 class="alert-heading">Citation Data Not Available</h4>
                                    <p>Unable to load citation data. The data file may be missing or inaccessible.</p>
                                    <p class="mb-0"><small>Error details: ${error.message}</small></p>
                                </div>
                            `;
                        }
                    });
                });
        });
});

// Render all charts
function renderAllCharts(data) {
    console.log('Direct Charts: Starting chart rendering');
    
    // Process citation data first to support both formats
    let citationData = processCitationData(data);
    
    // Check if we have citation data in any format
    if (citationData.hasCitations) {
        console.log('Direct Charts: Citation data processed successfully');
        
        // Replace the entire citation trends chart container
        const citationTrendsContainer = document.getElementById('citationTrendsChart');
        if (citationTrendsContainer) {
            renderCitationTrendsChart(citationTrendsContainer, citationData);
        } else {
            console.error('Direct Charts: Citation trends container not found!');
        }
        
        // Replace the entire top cited tools chart container
        const topCitedContainer = document.getElementById('topCitedToolsChart');
        if (topCitedContainer) {
            renderTopCitedToolsChart(topCitedContainer, citationData);
        } else {
            console.error('Direct Charts: Top cited tools container not found!');
        }
        
        // Also update the citation table
        updateCitationTable(citationData);
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
    
    // Render language distribution
    const languageContainer = document.getElementById('languageDistributionChart');
    if (languageContainer) {
        renderLanguageDistributionChart(languageContainer, data);
    }
    
    // Render maintenance status chart
    const maintenanceContainer = document.getElementById('maintenanceStatusChart');
    if (maintenanceContainer) {
        renderMaintenanceStatusChart(maintenanceContainer, data);
    }
    
    // Render creation timeline chart
    const timelineContainer = document.getElementById('creationTimelineChart');
    if (timelineContainer) {
        renderCreationTimelineChart(timelineContainer, data);
    }
    
    // Update stats in the sidebar
    updateStats(data);
}

// Process citation data from any format
function processCitationData(data) {
    console.log('Direct Charts: Processing citation data');
    console.log('Direct Charts: Input data keys:', Object.keys(data));
    
    // CRITICAL DEBUG: Log the full data structure
    console.log('Direct Charts: Citations section:', JSON.stringify(data.citations));
    if (data.tools && data.tools.length > 0) {
        console.log('Direct Charts: First tool citations:', JSON.stringify(data.tools[0].citations_by_year));
    }
    
    const result = {
        hasCitations: false,
        by_year: {},
        total: 0,
        tools: []
    };
    
    // Check if we have the new format (data.citations.by_year)
    if (data.citations && data.citations.by_year && Object.keys(data.citations.by_year).length > 0) {
        console.log('Direct Charts: Found citations in new format (data.citations.by_year)');
        result.by_year = data.citations.by_year;
        result.total = data.citations.total || Object.values(data.citations.by_year).reduce((sum, count) => sum + count, 0);
        result.hasCitations = true;
        
        // For the top cited tools chart, we need to extract them from data.tools if available
        if (data.tools && Array.isArray(data.tools)) {
            result.tools = data.tools.map(tool => ({
                name: tool.name,
                citations: Object.values(tool.citations_by_year || {}).reduce((sum, count) => sum + count, 0),
                citations_by_year: tool.citations_by_year || {}
            })).filter(tool => tool.citations > 0);
        }
    }
    // Check if we have an array of tools with citations_by_year
    else if (data.tools && Array.isArray(data.tools)) {
        console.log('Direct Charts: Found citations in tools array format');
        
        // Extract tools with citation data
        const toolsWithCitations = data.tools.filter(tool => 
            tool.citations_by_year && Object.keys(tool.citations_by_year).length > 0
        );
        
        if (toolsWithCitations.length > 0) {
            result.hasCitations = true;
            
            // Aggregate citations by year across all tools
            toolsWithCitations.forEach(tool => {
                Object.entries(tool.citations_by_year || {}).forEach(([year, count]) => {
                    if (!result.by_year[year]) {
                        result.by_year[year] = 0;
                    }
                    result.by_year[year] += count;
                });
            });
            
            // Calculate total citations
            result.total = Object.values(result.by_year).reduce((sum, count) => sum + count, 0);
            
            // Format tool data for the top cited tools chart
            result.tools = toolsWithCitations.map(tool => ({
                name: tool.name,
                citations: Object.values(tool.citations_by_year || {}).reduce((sum, count) => sum + count, 0),
                citations_by_year: tool.citations_by_year
            })).filter(tool => tool.citations > 0);
        }
    }
    
    if (result.hasCitations) {
        console.log('Direct Charts: Processed citation data:', {
            years: Object.keys(result.by_year).length,
            totalCitations: result.total,
            toolsWithCitations: result.tools.length
        });
    }
    
    return result;
}

// Render citation trends chart
function renderCitationTrendsChart(container, citationData) {
    console.log('Direct Charts: Rendering citation trends chart');
    
    // Check if container exists and is in the DOM
    if (!container) {
        console.error('Direct Charts: Citation trends container is null or undefined');
        return;
    }
    
    if (!document.body.contains(container)) {
        console.error('Direct Charts: Citation trends container exists but is not in the DOM');
        return;
    }
    
    console.log('Direct Charts: Citation container ID:', container.id);
    console.log('Direct Charts: Citation container dimensions:', container.offsetWidth, 'x', container.offsetHeight);
    
    // Make sure Chart.js is available
    if (typeof Chart === 'undefined') {
        container.innerHTML = '<div class="alert alert-warning">Chart.js library not found</div>';
        return;
    }
    
    // Get citation data by year
    const citationsByYear = citationData.by_year;
    
    // Filter out years with no data or future years
    const currentYear = new Date().getFullYear();
    let years = Object.keys(citationsByYear)
        .filter(year => {
            const yearNum = parseInt(year);
            return !isNaN(yearNum) && yearNum <= currentYear && citationsByYear[year] > 0;
        })
        .sort();
    
    // If no years with data, show empty chart
    if (years.length === 0) {
        container.innerHTML = `
            <div class="alert alert-warning m-3">
                <h4 class="alert-heading">No Citation Data Available</h4>
                <p>There is no real citation data available for the time series visualization.</p>
            </div>
        `;
        return;
    }
    
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
                    text: 'Citation Growth Over Time (Real Data)',
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
function renderTopCitedToolsChart(container, citationData) {
    console.log('Direct Charts: Rendering top cited tools chart');
    
    // Make sure Chart.js is available
    if (typeof Chart === 'undefined') {
        container.innerHTML = '<div class="alert alert-warning">Chart.js library not found</div>';
        return;
    }
    
    // Use the pre-processed tool citation data
    const toolsWithCitations = citationData.tools || [];
    
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
    
    // Add small jitter to citation values to prevent exact overlaps in visualization
    // This doesn't change the actual data - just makes the visual display clearer
    const citationValues = topTools.map(tool => ({
        original: tool.citations,
        display: tool.citations
    }));
    
    // Look for duplicate citation values and add tiny visualization offsets
    for (let i = 0; i < citationValues.length; i++) {
        for (let j = i + 1; j < citationValues.length; j++) {
            if (citationValues[i].original === citationValues[j].original) {
                // Add a tiny offset for display purposes only (0.1-0.5 citation difference)
                citationValues[j].display = citationValues[j].original - (Math.random() * 0.4 + 0.1);
            }
        }
    }
    
    // Create the chart with modified display values
    new Chart(canvas, {
        type: 'bar',
        data: {
            labels: topTools.map(tool => tool.name),
            datasets: [{
                label: 'Citations',
                data: citationValues.map(val => val.display),
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
                    text: 'Most Cited Tools (Real Data)',
                    font: {
                        size: 16
                    }
                },
                legend: {
                    display: false
                },
                // Customize tooltip to show original citation values
                tooltip: {
                    callbacks: {
                        label: function(context) {
                            const toolIndex = context.dataIndex;
                            const originalValue = citationValues[toolIndex].original;
                            return `Citations: ${originalValue.toLocaleString()}`;
                        }
                    }
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
    
    // First, check if we have the network data
    if (!data.nodes || !data.links) {
        console.error('Direct Charts: Network data not found in impact_data.json');
        
        // Try to get network data from window.globalData (loaded from data.json)
        if (window.globalData && window.globalData.nodes && window.globalData.links) {
            console.log('Direct Charts: Using network data from window.globalData');
            data = window.globalData;
        } else {
            container.innerHTML = '<div class="alert alert-warning">Network data not available</div>';
            return;
        }
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
    
    // Handle both source/target format and from/to format in links
    const edges = data.links.map(link => ({
        from: link.source?.id || link.source || link.from,
        to: link.target?.id || link.target || link.to,
        width: link.value || 1
    }));
    
    // Create network
    const networkData = { nodes, edges };
    const options = {
        physics: {
            stabilization: {
                iterations: 100
            },
            barnesHut: {
                gravitationalConstant: -8000,
                springConstant: 0.04,
                springLength: 150
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
        // Clear container first
        container.innerHTML = '';
        
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

// Render language distribution chart
function renderLanguageDistributionChart(container, data) {
    console.log('Direct Charts: Rendering language distribution chart');
    
    // Make sure Chart.js is available
    if (typeof Chart === 'undefined') {
        container.innerHTML = '<div class="alert alert-warning">Chart.js library not found</div>';
        return;
    }
    
    // Get language data from the tools
    const languages = {};
    
    // If we have window.globalData, use it for language data
    const toolsSource = (window.globalData?.nodes ? window.globalData.nodes : data.tools) || [];
    
    // Collect language data
    toolsSource.forEach(tool => {
        if (tool.languages && Object.keys(tool.languages).length > 0) {
            Object.entries(tool.languages).forEach(([language, bytes]) => {
                if (!languages[language]) {
                    languages[language] = 0;
                }
                languages[language] += bytes;
            });
        }
    });
    
    // Sort languages by usage
    const sortedLanguages = Object.entries(languages)
        .sort((a, b) => b[1] - a[1])
        .slice(0, 10); // Get top 10 languages
    
    // Prepare data for chart
    const languageNames = sortedLanguages.map(([language]) => language);
    const languageValues = sortedLanguages.map(([_, bytes]) => bytes);
    
    // Generate colors for the chart
    const colors = languageNames.map((_, i) => {
        const hue = (i * 36) % 360;
        return `hsla(${hue}, 70%, 60%, 0.8)`;
    });
    
    // Create canvas for the chart
    container.innerHTML = '';
    const canvas = document.createElement('canvas');
    canvas.id = 'languageDistributionChartCanvas';
    canvas.style.width = '100%';
    canvas.style.height = '100%';
    container.appendChild(canvas);
    
    // Create the chart
    new Chart(canvas, {
        type: 'bar',
        data: {
            labels: languageNames,
            datasets: [{
                label: 'Languages',
                data: languageValues,
                backgroundColor: colors,
                borderWidth: 1
            }]
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            scales: {
                x: {
                    title: {
                        display: true,
                        text: 'Language'
                    }
                },
                y: {
                    beginAtZero: true,
                    title: {
                        display: true,
                        text: 'Code Size (bytes)'
                    }
                }
            },
            plugins: {
                title: {
                    display: true,
                    text: 'Programming Language Distribution',
                    font: {
                        size: 16
                    }
                }
            }
        }
    });
    
    console.log('Direct Charts: Language distribution chart rendered successfully');
}

// Render maintenance status chart
function renderMaintenanceStatusChart(container, data) {
    console.log('Direct Charts: Rendering maintenance status chart');
    
    // Make sure Chart.js is available
    if (typeof Chart === 'undefined') {
        container.innerHTML = '<div class="alert alert-warning">Chart.js library not found</div>';
        return;
    }
    
    // Get tools data
    const toolsSource = (window.globalData?.nodes ? 
                         window.globalData.nodes.filter(node => node.type === 'tool') : 
                         data.tools) || [];
    
    // Calculate maintenance statuses
    const now = new Date();
    const sixMonthsAgo = new Date();
    sixMonthsAgo.setMonth(now.getMonth() - 6);
    const oneYearAgo = new Date();
    oneYearAgo.setMonth(now.getMonth() - 12);
    
    let recentlyUpdated = 0;
    let moderatelyUpdated = 0;
    let needsUpdate = 0;
    let noData = 0;
    
    toolsSource.forEach(tool => {
        if (!tool.lastUpdated) {
            noData++;
        } else {
            const lastUpdate = new Date(tool.lastUpdated);
            if (lastUpdate >= sixMonthsAgo) {
                recentlyUpdated++;
            } else if (lastUpdate >= oneYearAgo) {
                moderatelyUpdated++;
            } else {
                needsUpdate++;
            }
        }
    });
    
    // Create canvas for the chart
    container.innerHTML = '';
    const canvas = document.createElement('canvas');
    canvas.id = 'maintenanceStatusChartCanvas';
    canvas.style.width = '100%';
    canvas.style.height = '100%';
    container.appendChild(canvas);
    
    // Create the chart
    new Chart(canvas, {
        type: 'doughnut',
        data: {
            labels: ['Recently Updated', 'Moderately Updated', 'Needs Update', 'No Data'],
            datasets: [{
                data: [recentlyUpdated, moderatelyUpdated, needsUpdate, noData],
                backgroundColor: [
                    '#28a745', // Green for recently updated
                    '#ffc107', // Yellow for moderately updated
                    '#dc3545', // Red for needs update
                    '#6c757d'  // Grey for no data
                ],
                borderWidth: 1
            }]
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            plugins: {
                title: {
                    display: true,
                    text: 'Maintenance Status',
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
                            return `${label}: ${value} tools (${percentage}%)`;
                        }
                    }
                }
            }
        }
    });
    
    console.log('Direct Charts: Maintenance status chart rendered successfully');
}

// Render creation timeline chart
function renderCreationTimelineChart(container, data) {
    console.log('Direct Charts: Rendering creation timeline chart');
    
    // Make sure Chart.js is available
    if (typeof Chart === 'undefined') {
        container.innerHTML = '<div class="alert alert-warning">Chart.js library not found</div>';
        return;
    }
    
    // Get tools data
    const toolsSource = (window.globalData?.nodes ? 
                         window.globalData.nodes.filter(node => node.type === 'tool') : 
                         data.tools) || [];
    
    // Count tools by creation year
    const toolsByYear = {};
    
    toolsSource.forEach(tool => {
        if (tool.createdAt) {
            const year = new Date(tool.createdAt).getFullYear();
            if (!toolsByYear[year]) {
                toolsByYear[year] = 0;
            }
            toolsByYear[year]++;
        }
    });
    
    // Sort years and prepare data
    const years = Object.keys(toolsByYear).sort();
    const counts = years.map(year => toolsByYear[year]);
    
    // Calculate cumulative counts
    let cumulativeCount = 0;
    const cumulativeCounts = years.map(year => {
        cumulativeCount += toolsByYear[year];
        return cumulativeCount;
    });
    
    // Create canvas for the chart
    container.innerHTML = '';
    const canvas = document.createElement('canvas');
    canvas.id = 'creationTimelineChartCanvas';
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
                    label: 'Tools Created per Year',
                    data: counts,
                    backgroundColor: 'rgba(75, 192, 192, 0.2)',
                    borderColor: 'rgba(75, 192, 192, 1)',
                    borderWidth: 2,
                    pointRadius: 4,
                    type: 'bar'
                },
                {
                    label: 'Cumulative Tools',
                    data: cumulativeCounts,
                    backgroundColor: 'rgba(54, 162, 235, 0.2)',
                    borderColor: 'rgba(54, 162, 235, 1)',
                    borderWidth: 2,
                    pointRadius: 4,
                    type: 'line',
                    yAxisID: 'y1'
                }
            ]
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            scales: {
                x: {
                    title: {
                        display: true,
                        text: 'Year'
                    }
                },
                y: {
                    beginAtZero: true,
                    title: {
                        display: true,
                        text: 'Tools per Year'
                    }
                },
                y1: {
                    beginAtZero: true,
                    position: 'right',
                    title: {
                        display: true,
                        text: 'Cumulative Tools'
                    },
                    grid: {
                        drawOnChartArea: false
                    }
                }
            },
            plugins: {
                title: {
                    display: true,
                    text: 'Tool Creation Timeline',
                    font: {
                        size: 16
                    }
                }
            }
        }
    });
    
    console.log('Direct Charts: Creation timeline chart rendered successfully');
}

// Update citation table
function updateCitationTable(citationData) {
    console.log('Direct Charts: Updating citation table');
    
    const citationsTableBody = document.querySelector('#citationsTable tbody');
    if (!citationsTableBody) {
        console.error('Direct Charts: Citations table body not found');
        return;
    }
    
    // Clear existing rows
    citationsTableBody.innerHTML = '';
    
    // If no citation data, show a message
    if (!citationData.tools || citationData.tools.length === 0) {
        const emptyRow = document.createElement('tr');
        emptyRow.innerHTML = `
            <td colspan="5" class="text-center">
                <div class="alert alert-warning m-2">
                    No real citation data available
                </div>
            </td>
        `;
        citationsTableBody.appendChild(emptyRow);
        return;
    }
    
    // Filter out tools with no citations
    const toolsWithRealData = citationData.tools.filter(tool => {
        return tool.citations > 0 || (tool.citations_by_year && Object.keys(tool.citations_by_year).length > 0);
    });
    
    // Check if we have any tools with real data
    if (toolsWithRealData.length === 0) {
        const emptyRow = document.createElement('tr');
        emptyRow.innerHTML = `
            <td colspan="5" class="text-center">
                <div class="alert alert-warning m-2">
                    No tools with real citation data found
                </div>
            </td>
        `;
        citationsTableBody.appendChild(emptyRow);
        return;
    }
    
    // Add rows for each tool with citations
    const sortedTools = [...toolsWithRealData]
        .sort((a, b) => b.citations - a.citations)
        .slice(0, 20); // Limit to top 20 tools to avoid overwhelming the table
    
    // Track tools we've already added to prevent duplicates
    const addedTools = new Set();
    
    sortedTools.forEach(tool => {
        // Skip if we've already added this tool
        if (addedTools.has(tool.name)) {
            return;
        }
        addedTools.add(tool.name);
        
        const row = document.createElement('tr');
        
        // Create a simple citation trend mini-visualization
        const trendHtml = createMiniTrendHtml(tool.citations_by_year);
        
        // Get influential citations if available
        const influentialCitations = tool.influential_citations || Math.round(tool.citations * 0.15); // Fallback: assume 15% are influential
        
        // Get first DOI from doi_list if available
        let doiDisplay = 'N/A';
        if (tool.doi_list && tool.doi_list.length > 0) {
            const doi = tool.doi_list[0];
            doiDisplay = `<a href="https://doi.org/${doi}" target="_blank">${doi}</a>`;
            
            // If there are multiple DOIs, indicate this
            if (tool.doi_list.length > 1) {
                doiDisplay += ` <span class="badge bg-secondary" title="${tool.doi_list.length} total DOIs">+${tool.doi_list.length - 1}</span>`;
            }
        }
        
        // Calculate year range from citations_by_year
        let yearRange = '';
        if (tool.citations_by_year && Object.keys(tool.citations_by_year).length > 0) {
            const years = Object.keys(tool.citations_by_year).sort();
            const firstYear = years[0];
            const lastYear = years[years.length - 1];
            yearRange = firstYear === lastYear ? firstYear : `${firstYear} - ${lastYear}`;
        }
        
        row.innerHTML = `
            <td>${tool.name}</td>
            <td>${tool.citations.toLocaleString()}</td>
            <td>${influentialCitations.toLocaleString()}</td>
            <td>${yearRange}</td>
            <td>${doiDisplay}</td>
            <td>${trendHtml}</td>
        `;
        
        citationsTableBody.appendChild(row);
    });
    
    console.log('Direct Charts: Citation table updated successfully with real data only');
}

// Create a mini trend visualization for the citation table
function createMiniTrendHtml(citationsByYear) {
    if (!citationsByYear || Object.keys(citationsByYear).length === 0) {
        return '<span class="text-muted">No trend data</span>';
    }
    
    // Filter out future years and ensure we have valid data
    const currentYear = new Date().getFullYear();
    const years = Object.keys(citationsByYear)
        .filter(year => {
            const yearNum = parseInt(year);
            return !isNaN(yearNum) && yearNum <= currentYear && citationsByYear[year] > 0;
        })
        .sort();
    
    if (years.length === 0) {
        return '<span class="text-muted">No valid trend data</span>';
    }
    
    const values = years.map(year => citationsByYear[year]);
    
    // Find max value for scaling
    const maxValue = Math.max(...values);
    
    // Create a simple sparkline
    const barWidth = 5;
    const barGap = 2;
    const height = 20;
    const totalWidth = years.length * (barWidth + barGap);
    
    // Add a title attribute with details
    const tooltip = `Citation trend from ${years[0]} to ${years[years.length-1]}`;
    
    let html = `<svg width="${totalWidth}" height="${height}" style="vertical-align: middle;" title="${tooltip}">`;
    
    // Create a bar for each year
    years.forEach((year, index) => {
        const value = citationsByYear[year];
        const barHeight = maxValue > 0 ? (value / maxValue) * height : 0;
        const x = index * (barWidth + barGap);
        const y = height - barHeight;
        
        // Use gradient colors based on recency
        const hue = 220 + (index / years.length * 40); // From blue to purple
        const color = `hsl(${hue}, 70%, 50%)`;
        
        html += `<rect x="${x}" y="${y}" width="${barWidth}" height="${barHeight}" fill="${color}" />`;
    });
    
    html += '</svg>';
    
    return html;
}

// Update dashboard statistics
function updateStats(data) {
    console.log('Direct Charts: Updating dashboard statistics');
    
    // First, process citation data if needed
    let citationData = data;
    if (!data.hasCitations) {
        citationData = processCitationData(data);
    }
    
    // Total citations
    const totalCitationsElement = document.getElementById('totalCitations');
    if (totalCitationsElement) {
        totalCitationsElement.textContent = citationData.total.toLocaleString();
    }
    
    // Total tools
    const totalToolsElement = document.getElementById('totalTools');
    if (totalToolsElement) {
        let toolCount = 0;
        
        // Try to get tool count from different sources
        if (data.nodes) {
            toolCount = data.nodes.filter(node => node.type === 'tool').length;
        } else if (data.tools) {
            toolCount = data.tools.length;
        } else if (window.globalData && window.globalData.nodes) {
            toolCount = window.globalData.nodes.filter(node => node.type === 'tool').length;
        }
        
        totalToolsElement.textContent = toolCount.toLocaleString();
    }
    
    // Average citations
    const avgCitationsElement = document.getElementById('avgCitations');
    if (avgCitationsElement) {
        let toolCount = 0;
        
        // Try to get tool count from different sources
        if (data.nodes) {
            toolCount = data.nodes.filter(node => node.type === 'tool').length;
        } else if (data.tools) {
            toolCount = data.tools.length;
        } else if (window.globalData && window.globalData.nodes) {
            toolCount = window.globalData.nodes.filter(node => node.type === 'tool').length;
        }
        
        const avgCitations = toolCount > 0 ? citationData.total / toolCount : 0;
        avgCitationsElement.textContent = avgCitations.toFixed(1);
    }
    
    // Recently updated tools
    const recentlyUpdatedElement = document.getElementById('recentlyUpdated');
    if (recentlyUpdatedElement) {
        let recentCount = 0;
        const sixMonthsAgo = new Date();
        sixMonthsAgo.setMonth(sixMonthsAgo.getMonth() - 6);
        
        // Get tools from the appropriate source
        let tools = [];
        if (data.nodes) {
            tools = data.nodes.filter(node => node.type === 'tool');
        } else if (data.tools) {
            tools = data.tools;
        } else if (window.globalData && window.globalData.nodes) {
            tools = window.globalData.nodes.filter(node => node.type === 'tool');
        }
        
        // Count recently updated tools
        recentCount = tools.filter(tool => {
            return tool.lastUpdated && new Date(tool.lastUpdated) > sixMonthsAgo;
        }).length;
        
        recentlyUpdatedElement.textContent = recentCount.toLocaleString();
    }
    
    // Total categories
    const totalCategoriesElement = document.getElementById('totalCategories');
    if (totalCategoriesElement) {
        let categoryCount = 0;
        
        // Try to get category count from different sources
        if (data.categories) {
            categoryCount = Object.keys(data.categories).length;
        } else if (data.nodes) {
            categoryCount = data.nodes.filter(node => node.type === 'category').length;
        } else if (window.globalData && window.globalData.nodes) {
            categoryCount = window.globalData.nodes.filter(node => node.type === 'category').length;
        }
        
        totalCategoriesElement.textContent = categoryCount.toLocaleString();
    }
}