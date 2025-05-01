/**
 * Direct Embedded Citations - The simplest possible implementation
 * 
 * This file directly loads impact_data.json and renders the citation data
 * without any complex integration with the main dashboard code.
 */

// Wait for DOM to be ready
document.addEventListener('DOMContentLoaded', function() {
    console.log('Direct Embedded Citations: DOM loaded');
    
    // Get the citation section
    const citationSection = document.getElementById('citations-section');
    if (!citationSection) {
        console.error('Direct Embedded Citations: Citations section not found!');
        return;
    }
    
    // Replace the citation section content with our stub
    citationSection.innerHTML = `
        <h2 class="section-title">Academic Impact</h2>
        <p class="lead mb-4">
            Track the academic influence of virome analysis tools through citation metrics and trends.
        </p>
        
        <div class="row">
            <div class="col-12 mb-4">
                <div class="card">
                    <div class="card-header d-flex justify-content-between align-items-center">
                        <span>Top 10 Most Cited Tools</span>
                        <div id="direct-total-citations" class="badge bg-primary">Loading...</div>
                    </div>
                    <div class="card-body">
                        <div id="direct-loading" class="py-4 text-center">
                            <div class="spinner-border text-primary mb-3" role="status">
                                <span class="visually-hidden">Loading...</span>
                            </div>
                            <p>Loading citation data...</p>
                        </div>
                        <div id="direct-error" class="alert alert-warning" style="display: none;"></div>
                        <table id="direct-citations-table" class="table table-striped" style="display: none;">
                            <thead>
                                <tr>
                                    <th>Rank</th>
                                    <th>Tool</th>
                                    <th>Total Citations</th>
                                    <th>Years Cited</th>
                                </tr>
                            </thead>
                            <tbody>
                                <!-- Will be populated by JavaScript -->
                            </tbody>
                        </table>
                    </div>
                </div>
            </div>
        </div>
    `;
    
    // Now load and display the data
    loadDirectCitationData();
});

// Load citation data directly from impact_data.json
function loadDirectCitationData() {
    console.log('Direct Embedded Citations: Loading data');
    
    fetch('impact_data.json')
        .then(response => response.json())
        .then(data => {
            console.log('Direct Embedded Citations: Data loaded');
            processAndDisplayCitations(data);
        })
        .catch(error => {
            console.error('Direct Embedded Citations: Error loading data:', error);
            document.getElementById('direct-loading').style.display = 'none';
            
            const errorElem = document.getElementById('direct-error');
            errorElem.style.display = 'block';
            errorElem.innerHTML = `<p>Error loading citation data: ${error.message}</p>`;
        });
}

// Process and display citation data
function processAndDisplayCitations(data) {
    console.log('Direct Embedded Citations: Processing data');
    
    // First, hide loading indicator
    document.getElementById('direct-loading').style.display = 'none';
    
    // Show the table
    document.getElementById('direct-citations-table').style.display = 'table';
    
    // Process tools data
    if (!data.tools || !Array.isArray(data.tools) || data.tools.length === 0) {
        showError('No tools data found in impact_data.json');
        return;
    }
    
    // Calculate citation totals for each tool
    const toolsWithCitations = data.tools.map(tool => {
        // Sum up citations from all years
        let totalCitations = 0;
        const yearlyData = tool.citations_by_year || {};
        
        // Print out the yearly data for debugging
        console.log(`${tool.name} yearly citations:`, yearlyData);
        
        Object.values(yearlyData).forEach(count => {
            totalCitations += Number(count);
        });
        
        return {
            name: tool.name,
            totalCitations: totalCitations,
            yearSpan: Object.keys(yearlyData).length > 0 ? 
                `${Math.min(...Object.keys(yearlyData).map(Number))}-${Math.max(...Object.keys(yearlyData).map(Number))}` : 
                'N/A'
        };
    });
    
    // Sort by citation count (descending)
    toolsWithCitations.sort((a, b) => b.totalCitations - a.totalCitations);
    
    // Take top 10
    const top10Tools = toolsWithCitations.slice(0, 10);
    
    // Calculate total citations
    const totalCitations = toolsWithCitations.reduce((sum, tool) => sum + tool.totalCitations, 0);
    
    // Update total citations badge
    document.getElementById('direct-total-citations').textContent = `Total: ${totalCitations.toLocaleString()} Citations`;
    
    // Update the table
    const tableBody = document.querySelector('#direct-citations-table tbody');
    
    let tableHTML = '';
    top10Tools.forEach((tool, index) => {
        tableHTML += `
            <tr>
                <td>${index + 1}</td>
                <td>${tool.name}</td>
                <td>${tool.totalCitations.toLocaleString()}</td>
                <td>${tool.yearSpan}</td>
            </tr>
        `;
    });
    
    tableBody.innerHTML = tableHTML;
    
    console.log('Direct Embedded Citations: Display complete');
    console.log('Top 10 tools:', top10Tools);
}

// Show error message
function showError(message) {
    const errorElem = document.getElementById('direct-error');
    errorElem.style.display = 'block';
    errorElem.innerHTML = `<p>${message}</p>`;
}