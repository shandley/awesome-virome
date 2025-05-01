/**
 * Citation Charts - Direct visualization for dashboard citations
 * This file provides a standalone implementation of citation charts
 * that doesn't depend on the main dashboard initialization.
 */

// Function to load impact data and render charts
function loadAndRenderCitationCharts() {
    console.log('Citation charts: Starting initialization');
    
    // Check if Chart.js is available
    if (typeof Chart === 'undefined') {
        console.error('Citation charts: Chart.js library is not loaded');
        document.getElementById('citationTrendsChart').innerHTML = `
            <div class="d-flex justify-content-center align-items-center h-100">
                <div class="text-center">
                    <i class="bi bi-exclamation-circle mb-3" style="font-size: 3rem;"></i>
                    <h4>Chart Library Missing</h4>
                    <p class="text-muted">Required Chart.js library not loaded.</p>
                </div>
            </div>
        `;
        return;
    }
    
    // Load impact data
    fetch('impact_data.json')
        .then(response => {
            if (!response.ok) {
                throw new Error(`HTTP error! status: ${response.status}`);
            }
            return response.json();
        })
        .then(impactData => {
            console.log('Citation charts: Successfully loaded impact_data.json');
            console.log('Citation charts: Found', impactData.tools ? impactData.tools.length : 0, 'tools with citation data');
            
            // Render citation trends chart
            renderCitationTrendsChart(impactData);
            
            // Render top cited tools chart
            renderTopCitedToolsChart(impactData);
        })
        .catch(error => {
            console.error('Citation charts: Error loading impact data:', error);
            
            // Show error message in charts
            document.getElementById('citationTrendsChart').innerHTML = `
                <div class="d-flex justify-content-center align-items-center h-100">
                    <div class="text-center">
                        <i class="bi bi-exclamation-circle mb-3" style="font-size: 3rem;"></i>
                        <h4>Data Loading Error</h4>
                        <p class="text-muted">Could not load citation data: ${error.message}</p>
                    </div>
                </div>
            `;
            
            document.getElementById('topCitedToolsChart').innerHTML = `
                <div class="d-flex justify-content-center align-items-center h-100">
                    <div class="text-center">
                        <i class="bi bi-exclamation-circle mb-3" style="font-size: 3rem;"></i>
                        <h4>Data Loading Error</h4>
                        <p class="text-muted">Could not load citation data: ${error.message}</p>
                    </div>
                </div>
            `;
        });
}

// Render citation trends chart
function renderCitationTrendsChart(impactData) {
    const container = document.getElementById('citationTrendsChart');
    if (!container) {
        console.error('Citation charts: Container #citationTrendsChart not found');
        return;
    }
    
    console.log('Citation charts: Rendering citation trends chart');
    
    // Get citation data by year
    const citationsByYear = impactData.citations && impactData.citations.by_year ? 
        impactData.citations.by_year : {};
    
    if (Object.keys(citationsByYear).length === 0) {
        container.innerHTML = `
            <div class="d-flex justify-content-center align-items-center h-100">
                <div class="text-center">
                    <i class="bi bi-graph-up mb-3" style="font-size: 3rem;"></i>
                    <h4>Citation Trends</h4>
                    <p class="text-muted">No citation trend data available.</p>
                </div>
            </div>
        `;
        return;
    }
    
    // Get years in sequence
    const years = Object.keys(citationsByYear).sort();
    
    // Calculate cumulative citations
    let cumulativeCount = 0;
    const cumulativeData = years.map(year => {
        cumulativeCount += citationsByYear[year];
        return cumulativeCount;
    });
    
    // Create canvas for Chart.js
    const ctx = document.createElement('canvas');
    ctx.width = container.clientWidth;
    ctx.height = container.clientHeight;
    container.innerHTML = '';
    container.appendChild(ctx);
    
    // Create chart
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
    
    console.log('Citation charts: Citation trends chart rendered successfully');
}

// Render top cited tools chart
function renderTopCitedToolsChart(impactData) {
    const container = document.getElementById('topCitedToolsChart');
    if (!container) {
        console.error('Citation charts: Container #topCitedToolsChart not found');
        return;
    }
    
    console.log('Citation charts: Rendering top cited tools chart');
    
    // Get tools with citation data
    const tools = impactData.tools || [];
    
    if (tools.length === 0) {
        container.innerHTML = `
            <div class="d-flex justify-content-center align-items-center h-100">
                <div class="text-center">
                    <i class="bi bi-bar-chart mb-3" style="font-size: 3rem;"></i>
                    <h4>Most Cited Tools</h4>
                    <p class="text-muted">No citation data available.</p>
                </div>
            </div>
        `;
        return;
    }
    
    // Calculate total citations for each tool
    const toolsWithTotals = tools.map(tool => ({
        name: tool.name,
        total: Object.values(tool.citations_by_year || {}).reduce((sum, count) => sum + count, 0)
    }));
    
    // Sort by total citations and take top 10
    const topTools = toolsWithTotals
        .sort((a, b) => b.total - a.total)
        .slice(0, 10);
    
    // Create canvas for Chart.js
    const ctx = document.createElement('canvas');
    ctx.width = container.clientWidth;
    ctx.height = container.clientHeight;
    container.innerHTML = '';
    container.appendChild(ctx);
    
    // Generate colors
    const colors = topTools.map((_, i) => {
        const hue = (i * 25) % 360;
        return `hsla(${hue}, 70%, 60%, 0.8)`;
    });
    
    // Create chart
    new Chart(ctx, {
        type: 'bar',
        data: {
            labels: topTools.map(tool => tool.name),
            datasets: [{
                label: 'Citations',
                data: topTools.map(tool => tool.total),
                backgroundColor: colors,
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
    
    console.log('Citation charts: Top cited tools chart rendered successfully');
}

// Initialize when document is fully loaded
document.addEventListener('DOMContentLoaded', function() {
    console.log('Citation charts: DOM loaded, initializing charts');
    
    // Let the page finish loading other elements first
    setTimeout(loadAndRenderCitationCharts, 500);
});