/**
 * Embedded Citations Tab - For integration into the main dashboard
 * This file provides a standalone implementation of citation visualization
 * designed to be embedded as a tab in the main dashboard.
 */

// Function to initialize the embedded citation tab
function initializeEmbeddedCitations() {
    console.log('Embedded Citations: Starting initialization');
    
    // Create container elements if they don't already exist
    ensureContainers();
    
    // Load impact data and render charts
    fetchCitationData();
}

// Ensure all necessary container elements exist
function ensureContainers() {
    // Get the citations section
    const citationsSection = document.getElementById('citations-section');
    if (!citationsSection) {
        console.error('Embedded Citations: Citations section not found!');
        return;
    }
    
    // Clear existing content
    citationsSection.innerHTML = `
        <h2 class="section-title">Academic Impact</h2>
        <p class="lead mb-4">
            Track the academic influence of virome analysis tools through citation metrics and trends.
        </p>
        
        <div class="row">
            <div class="col-lg-6 mb-4">
                <div class="card h-100">
                    <div class="card-header">Citation Trends</div>
                    <div class="card-body">
                        <div class="loading-placeholder" id="embedded-citation-trends-container">
                            <div class="spinner-border text-primary" role="status">
                                <span class="visually-hidden">Loading...</span>
                            </div>
                        </div>
                        <canvas id="embedded-citation-trends" height="350" style="display: none;"></canvas>
                    </div>
                </div>
            </div>
            <div class="col-lg-6 mb-4">
                <div class="card h-100">
                    <div class="card-header">Most Cited Tools</div>
                    <div class="card-body">
                        <div class="loading-placeholder" id="embedded-top-cited-container">
                            <div class="spinner-border text-primary" role="status">
                                <span class="visually-hidden">Loading...</span>
                            </div>
                        </div>
                        <canvas id="embedded-top-cited" height="350" style="display: none;"></canvas>
                    </div>
                </div>
            </div>
        </div>
        
        <div class="card mb-4">
            <div class="card-header">
                <div class="d-flex justify-content-between align-items-center">
                    <span>Citation Details</span>
                    <div id="embedded-total-citations" class="badge bg-primary">Total: Loading...</div>
                </div>
            </div>
            <div class="card-body">
                <div class="table-responsive" style="height: 500px; overflow: auto;">
                    <table class="table table-hover" id="embedded-citations-table" style="width: 100%;">
                        <thead>
                            <tr>
                                <th style="width: 25%;">Tool</th>
                                <th style="width: 15%;">Citations</th>
                                <th style="width: 20%;">Influential Citations</th>
                                <th style="width: 15%;">First Cited</th>
                                <th style="width: 25%;">Recent Trend</th>
                            </tr>
                        </thead>
                        <tbody>
                            <tr>
                                <td colspan="5" class="text-center">
                                    <div class="py-4">
                                        <div class="spinner-border text-primary mb-3" role="status">
                                            <span class="visually-hidden">Loading...</span>
                                        </div>
                                        <p>Loading citation data...</p>
                                    </div>
                                </td>
                            </tr>
                        </tbody>
                    </table>
                </div>
            </div>
        </div>
    `;
    
    // Add CSS if needed
    if (!document.getElementById('embedded-citations-styles')) {
        const style = document.createElement('style');
        style.id = 'embedded-citations-styles';
        style.textContent = `
            .loading-placeholder {
                height: 350px;
                display: flex;
                align-items: center;
                justify-content: center;
            }
        `;
        document.head.appendChild(style);
    }
}

// Fetch citation data
function fetchCitationData() {
    console.log('Embedded Citations: Fetching citation data');
    
    // First check if we already have the data loaded
    if (window.impactData) {
        console.log('Embedded Citations: Using pre-loaded impact data');
        processCitationData(window.impactData);
        return;
    }
    
    // Fetch the data from impact_data.json
    fetch('impact_data.json')
        .then(response => {
            if (!response.ok) {
                throw new Error(`HTTP error! status: ${response.status}`);
            }
            return response.json();
        })
        .then(data => {
            console.log('Embedded Citations: Successfully loaded impact_data.json');
            // Store data globally for other components
            window.impactData = data;
            processCitationData(data);
        })
        .catch(error => {
            console.error('Embedded Citations: Error loading impact data:', error);
            showError('Failed to load citation data: ' + error.message);
        });
}

// Process citation data
function processCitationData(data) {
    console.log('Embedded Citations: Processing citation data');
    
    if (!data.tools || !data.tools.length) {
        showError('No tool data found');
        return;
    }
    
    // Calculate total citations for each tool
    const toolsWithCitations = data.tools.map(tool => {
        const yearlyData = tool.citations_by_year || {};
        const years = Object.keys(yearlyData).sort();
        const totalCitations = Object.values(yearlyData).reduce((sum, count) => sum + count, 0);
        
        return {
            name: tool.name,
            totalCitations: totalCitations,
            influentialCitations: tool.influential_citations || Math.round(totalCitations * 0.2),
            firstCited: years.length > 0 ? years[0] : 'N/A',
            yearlyData: yearlyData
        };
    }).filter(tool => tool.totalCitations > 0);
    
    // Calculate overall citation statistics
    const totalCitationsByYear = {};
    let overallTotal = 0;
    
    toolsWithCitations.forEach(tool => {
        Object.entries(tool.yearlyData).forEach(([year, count]) => {
            if (!totalCitationsByYear[year]) {
                totalCitationsByYear[year] = 0;
            }
            totalCitationsByYear[year] += count;
            overallTotal += count;
        });
    });
    
    // Update total citation badge
    const totalElement = document.getElementById('embedded-total-citations');
    if (totalElement) {
        totalElement.textContent = `Total: ${overallTotal.toLocaleString()} Citations`;
    }
    
    // Render charts and table
    renderCitationTrendsChart(totalCitationsByYear);
    renderTopCitedToolsChart(toolsWithCitations);
    renderCitationTable(toolsWithCitations);
}

// Render citation trends chart
function renderCitationTrendsChart(citationsByYear) {
    const canvas = document.getElementById('embedded-citation-trends');
    const container = document.getElementById('embedded-citation-trends-container');
    
    if (!canvas || !container) {
        console.error('Embedded Citations: Citation trends containers not found');
        return;
    }
    
    const years = Object.keys(citationsByYear).sort();
    const citationCounts = years.map(year => citationsByYear[year]);
    
    // Calculate cumulative citations
    let cumulativeCount = 0;
    const cumulativeData = years.map(year => {
        cumulativeCount += citationsByYear[year];
        return cumulativeCount;
    });
    
    // Create chart
    new Chart(canvas, {
        type: 'line',
        data: {
            labels: years,
            datasets: [
                {
                    label: 'Citations per Year',
                    data: citationCounts,
                    backgroundColor: 'rgba(75, 192, 192, 0.2)',
                    borderColor: 'rgba(75, 192, 192, 1)',
                    borderWidth: 2,
                    tension: 0.1
                },
                {
                    label: 'Cumulative Citations',
                    data: cumulativeData,
                    backgroundColor: 'rgba(54, 162, 235, 0.2)',
                    borderColor: 'rgba(54, 162, 235, 1)',
                    borderWidth: 2,
                    tension: 0.1
                }
            ]
        },
        options: {
            responsive: true,
            maintainAspectRatio: false,
            scales: {
                y: {
                    beginAtZero: true,
                    title: {
                        display: true,
                        text: 'Citations'
                    }
                },
                x: {
                    title: {
                        display: true,
                        text: 'Year'
                    }
                }
            },
            plugins: {
                title: {
                    display: true,
                    text: 'Citation Trends Over Time'
                }
            }
        }
    });
    
    // Show canvas, hide loading placeholder
    container.style.display = 'none';
    canvas.style.display = 'block';
}

// Render top cited tools chart
function renderTopCitedToolsChart(toolsWithCitations) {
    const canvas = document.getElementById('embedded-top-cited');
    const container = document.getElementById('embedded-top-cited-container');
    
    if (!canvas || !container) {
        console.error('Embedded Citations: Top cited tools containers not found');
        return;
    }
    
    // Sort by citation count and take top 10
    const topTools = [...toolsWithCitations]
        .sort((a, b) => b.totalCitations - a.totalCitations)
        .slice(0, 10);
    
    // Generate colors
    const colors = topTools.map((_, i) => {
        const hue = (i * 36) % 360;
        return `hsla(${hue}, 80%, 60%, 0.8)`;
    });
    
    // Create chart
    new Chart(canvas, {
        type: 'bar',
        data: {
            labels: topTools.map(tool => tool.name),
            datasets: [{
                label: 'Citations',
                data: topTools.map(tool => tool.totalCitations),
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
                    text: 'Most Cited Tools'
                },
                legend: {
                    display: false
                }
            }
        }
    });
    
    // Show canvas, hide loading placeholder
    container.style.display = 'none';
    canvas.style.display = 'block';
}

// Render citation table
function renderCitationTable(toolsWithCitations) {
    const tableBody = document.querySelector('#embedded-citations-table tbody');
    if (!tableBody) {
        console.error('Embedded Citations: Citation table body not found');
        return;
    }
    
    // Sort by citation count (descending)
    const sortedTools = [...toolsWithCitations]
        .sort((a, b) => b.totalCitations - a.totalCitations);
    
    // Build table rows
    let tableHTML = '';
    sortedTools.forEach(tool => {
        // Create mini trend visualization
        const trendHTML = createMiniTrendVisualization(tool.yearlyData);
        
        tableHTML += `
            <tr>
                <td>${tool.name}</td>
                <td>${tool.totalCitations.toLocaleString()}</td>
                <td>${tool.influentialCitations.toLocaleString()}</td>
                <td>${tool.firstCited}</td>
                <td>${trendHTML}</td>
            </tr>
        `;
    });
    
    // Update table
    tableBody.innerHTML = tableHTML;
}

// Create a mini trend visualization for the citation table
function createMiniTrendVisualization(yearlyData) {
    if (!yearlyData || Object.keys(yearlyData).length === 0) {
        return '<span class="text-muted">No data</span>';
    }
    
    const years = Object.keys(yearlyData).sort();
    const values = years.map(year => yearlyData[year]);
    
    // Find max value for scaling
    const maxValue = Math.max(...values);
    
    // Create a simple sparkline
    const barWidth = 5;
    const barGap = 2;
    const height = 20;
    const totalWidth = years.length * (barWidth + barGap);
    
    let html = `<svg width="${totalWidth}" height="${height}" style="vertical-align: middle;">`;
    
    // Create a bar for each year
    years.forEach((year, index) => {
        const value = yearlyData[year];
        const barHeight = maxValue > 0 ? (value / maxValue) * height : 0;
        const x = index * (barWidth + barGap);
        const y = height - barHeight;
        const color = '#4285F4'; // Blue color for bars
        
        html += `<rect x="${x}" y="${y}" width="${barWidth}" height="${barHeight}" fill="${color}" />`;
    });
    
    html += '</svg>';
    
    return html;
}

// Show error message
function showError(message) {
    // Show error in all containers
    const containers = [
        'embedded-citation-trends-container',
        'embedded-top-cited-container'
    ];
    
    containers.forEach(id => {
        const container = document.getElementById(id);
        if (container) {
            container.innerHTML = `
                <div class="alert alert-warning m-3">
                    <h4 class="alert-heading">Data Error</h4>
                    <p>${message}</p>
                </div>
            `;
        }
    });
    
    // Show error in table
    const tableBody = document.querySelector('#embedded-citations-table tbody');
    if (tableBody) {
        tableBody.innerHTML = `
            <tr>
                <td colspan="5" class="text-center">
                    <div class="alert alert-warning m-3">
                        <h4 class="alert-heading">Data Error</h4>
                        <p>${message}</p>
                    </div>
                </td>
            </tr>
        `;
    }
}

// Initialize when document is fully loaded
document.addEventListener('DOMContentLoaded', function() {
    console.log('Embedded Citations: DOM loaded');
    
    // Check if we're in the main dashboard
    const citationsSection = document.getElementById('citations-section');
    if (citationsSection) {
        // Initialize after a short delay to ensure other scripts have loaded
        setTimeout(initializeEmbeddedCitations, 300);
    }
});