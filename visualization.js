// Awesome-Virome Dashboard Visualization
document.addEventListener('DOMContentLoaded', function() {
    // Global variables
    let allData = null;
    let filteredData = null;
    let selectedTool = null;
    
    // Status colors
    const statusColors = {
        recent: '#198754',    // green - updated within 6 months
        moderate: '#ffc107',  // yellow - updated within 12 months
        stale: '#dc3545',     // red - older than 12 months
        unknown: '#6c757d'    // gray - no update data
    };
    
    // Color scale for categories
    const categoryColorScale = d3.scaleOrdinal()
        .range(d3.schemeCategory10);
    
    // Initialize the dashboard
    init();
    
    // Main initialization function
    function init() {
        // Load data
        loadData();
        
        // Set up event listeners
        setupEventListeners();
    }
    
    // Load data from JSON file
    function loadData() {
        d3.json('data.json')
            .then(data => {
                console.log('Data loaded successfully:', data);
                
                // Process the data
                processData(data);
                
                // Show the dashboard content
                showDashboard();
            })
            .catch(error => {
                console.error('Error loading data:', error);
                showError();
            });
    }
    
    // Process the raw data
    function processData(data) {
        // Extract tool nodes from the data
        if (data.nodes) {
            allData = data.nodes.filter(node => node.type === 'tool');
        } else {
            // If data is already in the right format
            allData = data;
        }
        
        // Set the filtered data initially to all data
        filteredData = [...allData];
        
        // Populate filter dropdowns
        populateFilters();
        
        // Initialize summary statistics
        updateSummaryStats();
        
        // Create charts
        createCharts();
        
        // Populate top tools
        updateTopTools();
        
        // Populate data table
        updateToolsTable();
    }
    
    // Populate filter dropdowns
    function populateFilters() {
        // Get unique categories
        const categories = [...new Set(allData.map(tool => tool.category))].sort();
        const categorySelect = document.getElementById('category-filter');
        
        categories.forEach(category => {
            if (category) {
                const option = document.createElement('option');
                option.value = category;
                option.textContent = category;
                categorySelect.appendChild(option);
            }
        });
        
        // Get unique languages
        const languages = [...new Set(allData.map(tool => tool.language))].sort();
        const languageSelect = document.getElementById('language-filter');
        
        languages.forEach(language => {
            if (language) {
                const option = document.createElement('option');
                option.value = language;
                option.textContent = language;
                languageSelect.appendChild(option);
            }
        });
        
        // Update the max value for stars filter based on data
        const maxStars = Math.max(...allData.map(tool => tool.stars || 0));
        const starsFilter = document.getElementById('stars-filter');
        starsFilter.max = Math.min(Math.ceil(maxStars / 100) * 100, 1000); // Round up to nearest 100, max 1000
    }
    
    // Set up event listeners
    function setupEventListeners() {
        // Category filter
        document.getElementById('category-filter').addEventListener('change', applyFilters);
        
        // Language filter
        document.getElementById('language-filter').addEventListener('change', applyFilters);
        
        // Update time filter
        document.getElementById('update-filter').addEventListener('change', applyFilters);
        
        // Stars filter
        const starsFilter = document.getElementById('stars-filter');
        const starsValue = document.getElementById('stars-value');
        
        starsFilter.addEventListener('input', function() {
            starsValue.textContent = this.value;
        });
        
        starsFilter.addEventListener('change', applyFilters);
        
        // Search input
        document.getElementById('search-input').addEventListener('input', debounce(applyFilters, 300));
        
        // Reset filters
        document.getElementById('reset-filters').addEventListener('click', resetFilters);
    }
    
    // Apply all filters
    function applyFilters() {
        // Get filter values
        const categoryFilter = document.getElementById('category-filter').value;
        const languageFilter = document.getElementById('language-filter').value;
        const updateFilter = document.getElementById('update-filter').value;
        const starsFilter = parseInt(document.getElementById('stars-filter').value);
        const searchTerm = document.getElementById('search-input').value.toLowerCase();
        
        // Filter the data
        filteredData = allData.filter(tool => {
            // Category filter
            if (categoryFilter !== 'all' && tool.category !== categoryFilter) {
                return false;
            }
            
            // Language filter
            if (languageFilter !== 'all' && tool.language !== languageFilter) {
                return false;
            }
            
            // Update time filter
            if (updateFilter !== 'all') {
                if (tool.updateTime === null || tool.updateTime === undefined) {
                    return false;
                }
                if (parseInt(updateFilter) < tool.updateTime) {
                    return false;
                }
            }
            
            // Stars filter
            if (tool.stars < starsFilter) {
                return false;
            }
            
            // Search term
            if (searchTerm) {
                const searchFields = [
                    tool.name,
                    tool.description,
                    tool.category,
                    tool.subcategory,
                    tool.language
                ].filter(Boolean).map(field => field.toLowerCase());
                
                return searchFields.some(field => field.includes(searchTerm));
            }
            
            return true;
        });
        
        // Update the dashboard with filtered data
        updateDashboard();
    }
    
    // Reset all filters
    function resetFilters() {
        // Reset filter inputs
        document.getElementById('category-filter').value = 'all';
        document.getElementById('language-filter').value = 'all';
        document.getElementById('update-filter').value = 'all';
        document.getElementById('stars-filter').value = 0;
        document.getElementById('stars-value').textContent = '0';
        document.getElementById('search-input').value = '';
        
        // Reset filtered data
        filteredData = [...allData];
        
        // Update the dashboard
        updateDashboard();
    }
    
    // Update the entire dashboard
    function updateDashboard() {
        // Update summary statistics
        updateSummaryStats();
        
        // Update charts
        updateCharts();
        
        // Update top tools
        updateTopTools();
        
        // Update data table
        updateToolsTable();
    }
    
    // Update summary statistics
    function updateSummaryStats() {
        // Total tools
        document.getElementById('total-tools').textContent = filteredData.length;
        
        // Total categories
        const categories = [...new Set(filteredData.map(tool => tool.category))];
        document.getElementById('total-categories').textContent = categories.length;
        
        // Average stars
        const totalStars = filteredData.reduce((sum, tool) => sum + (tool.stars || 0), 0);
        const avgStars = filteredData.length > 0 ? Math.round(totalStars / filteredData.length) : 0;
        document.getElementById('avg-stars').textContent = avgStars;
        
        // Recently updated (within 6 months)
        const recentlyUpdated = filteredData.filter(tool => 
            tool.updateTime !== null && tool.updateTime !== undefined && tool.updateTime <= 6
        ).length;
        document.getElementById('recently-updated').textContent = recentlyUpdated;
        
        // Update filtered count
        document.getElementById('filtered-count').textContent = `${filteredData.length} of ${allData.length} tools`;
    }
    
    // Create all charts
    function createCharts() {
        createCategoryChart();
        createMaintenanceChart();
        createLanguageChart();
        createStarsChart();
        createTimelineChart();
    }
    
    // Update all charts
    function updateCharts() {
        updateCategoryChart();
        updateMaintenanceChart();
        updateLanguageChart();
        updateStarsChart();
        updateTimelineChart();
    }
    
    // Create category distribution chart
    function createCategoryChart() {
        const chartContainer = document.getElementById('category-chart');
        const width = chartContainer.clientWidth;
        const height = chartContainer.clientHeight;
        const margin = { top: 30, right: 30, bottom: 70, left: 60 };
        
        // Create SVG
        const svg = d3.select(chartContainer)
            .append('svg')
            .attr('width', width)
            .attr('height', height)
            .append('g')
            .attr('transform', `translate(${margin.left},${margin.top})`);
        
        // Create scales (will be updated in updateCategoryChart)
        const x = d3.scaleBand()
            .range([0, width - margin.left - margin.right])
            .padding(0.2);
        
        const y = d3.scaleLinear()
            .range([height - margin.top - margin.bottom, 0]);
        
        // Create axes
        svg.append('g')
            .attr('class', 'x-axis')
            .attr('transform', `translate(0,${height - margin.top - margin.bottom})`);
        
        svg.append('g')
            .attr('class', 'y-axis');
        
        // Create title
        svg.append('text')
            .attr('class', 'chart-title')
            .attr('x', (width - margin.left - margin.right) / 2)
            .attr('y', -10)
            .attr('text-anchor', 'middle')
            .style('font-size', '14px')
            .style('font-weight', 'bold');
        
        // Initial update
        updateCategoryChart();
    }
    
    // Update category distribution chart
    function updateCategoryChart() {
        const chartContainer = document.getElementById('category-chart');
        const width = chartContainer.clientWidth;
        const height = chartContainer.clientHeight;
        const margin = { top: 30, right: 30, bottom: 70, left: 60 };
        
        // Calculate category counts
        const categoryCounts = {};
        filteredData.forEach(tool => {
            if (tool.category) {
                categoryCounts[tool.category] = (categoryCounts[tool.category] || 0) + 1;
            }
        });
        
        // Convert to array and sort
        const data = Object.entries(categoryCounts)
            .map(([category, count]) => ({ category, count }))
            .sort((a, b) => b.count - a.count);
        
        // Update scales
        const x = d3.scaleBand()
            .domain(data.map(d => d.category))
            .range([0, width - margin.left - margin.right])
            .padding(0.2);
        
        const y = d3.scaleLinear()
            .domain([0, d3.max(data, d => d.count)])
            .nice()
            .range([height - margin.top - margin.bottom, 0]);
        
        // Update axes
        const svg = d3.select(chartContainer).select('svg').select('g');
        
        svg.select('.x-axis')
            .transition()
            .duration(500)
            .call(d3.axisBottom(x)
                .tickSizeOuter(0))
            .selectAll('text')
            .style('text-anchor', 'end')
            .attr('dx', '-.8em')
            .attr('dy', '.15em')
            .attr('transform', 'rotate(-45)');
        
        svg.select('.y-axis')
            .transition()
            .duration(500)
            .call(d3.axisLeft(y)
                .ticks(5)
                .tickSizeOuter(0));
        
        // Update bars
        const bars = svg.selectAll('.category-bar')
            .data(data, d => d.category);
        
        // Remove old bars
        bars.exit()
            .transition()
            .duration(300)
            .attr('y', height - margin.top - margin.bottom)
            .attr('height', 0)
            .remove();
        
        // Add new bars
        bars.enter()
            .append('rect')
            .attr('class', 'category-bar bar')
            .attr('x', d => x(d.category))
            .attr('width', x.bandwidth())
            .attr('y', height - margin.top - margin.bottom)
            .attr('height', 0)
            .attr('fill', d => categoryColorScale(d.category))
            .on('mouseover', function(event, d) {
                showTooltip(event, `<strong>${d.category}</strong><br>${d.count} tools`);
            })
            .on('mouseout', hideTooltip)
            .transition()
            .duration(500)
            .attr('y', d => y(d.count))
            .attr('height', d => height - margin.top - margin.bottom - y(d.count));
        
        // Update existing bars
        bars.transition()
            .duration(500)
            .attr('x', d => x(d.category))
            .attr('width', x.bandwidth())
            .attr('y', d => y(d.count))
            .attr('height', d => height - margin.top - margin.bottom - y(d.count))
            .attr('fill', d => categoryColorScale(d.category));
    }
    
    // Create maintenance status chart (donut chart)
    function createMaintenanceChart() {
        const chartContainer = document.getElementById('maintenance-chart');
        const width = chartContainer.clientWidth;
        const height = chartContainer.clientHeight;
        
        // Create SVG
        const svg = d3.select(chartContainer)
            .append('svg')
            .attr('width', width)
            .attr('height', height)
            .append('g')
            .attr('transform', `translate(${width / 2},${height / 2})`);
        
        // Create donut chart group
        svg.append('g')
            .attr('class', 'donut-chart');
        
        // Create center text
        svg.append('text')
            .attr('class', 'donut-center-text')
            .attr('text-anchor', 'middle')
            .attr('dy', '0.35em')
            .style('font-size', '16px')
            .style('font-weight', 'bold');
        
        // Create legend
        svg.append('g')
            .attr('class', 'donut-legend')
            .attr('transform', `translate(0,${height / 2 - 40})`);
        
        // Initial update
        updateMaintenanceChart();
    }
    
    // Update maintenance status chart
    function updateMaintenanceChart() {
        const chartContainer = document.getElementById('maintenance-chart');
        const width = chartContainer.clientWidth;
        const height = chartContainer.clientHeight;
        const radius = Math.min(width, height) / 2 - 40;
        
        // Calculate status counts
        const statusCounts = {
            recent: 0,
            moderate: 0,
            stale: 0,
            unknown: 0
        };
        
        filteredData.forEach(tool => {
            if (tool.updateTime === null || tool.updateTime === undefined) {
                statusCounts.unknown++;
            } else if (tool.updateTime <= 6) {
                statusCounts.recent++;
            } else if (tool.updateTime <= 12) {
                statusCounts.moderate++;
            } else {
                statusCounts.stale++;
            }
        });
        
        // Convert to array
        const data = [
            { status: 'Recent (≤6mo)', count: statusCounts.recent, color: statusColors.recent },
            { status: 'Moderate (≤12mo)', count: statusCounts.moderate, color: statusColors.moderate },
            { status: 'Stale (>12mo)', count: statusCounts.stale, color: statusColors.stale },
            { status: 'Unknown', count: statusCounts.unknown, color: statusColors.unknown }
        ].filter(d => d.count > 0);
        
        // Create pie layout
        const pie = d3.pie()
            .value(d => d.count)
            .sort(null);
        
        // Create arc generator
        const arc = d3.arc()
            .innerRadius(radius * 0.6)
            .outerRadius(radius);
        
        // Create SVG and update donut chart
        const svg = d3.select(chartContainer).select('svg').select('g');
        
        // Update donut segments
        const segments = svg.select('.donut-chart')
            .selectAll('path')
            .data(pie(data), d => d.data.status);
        
        // Remove old segments
        segments.exit().remove();
        
        // Add new segments
        segments.enter()
            .append('path')
            .attr('class', 'donut-segment')
            .attr('fill', d => d.data.color)
            .attr('d', arc)
            .on('mouseover', function(event, d) {
                showTooltip(event, `<strong>${d.data.status}</strong><br>${d.data.count} tools (${Math.round(d.data.count / filteredData.length * 100)}%)`);
                
                // Update center text
                svg.select('.donut-center-text')
                    .text(`${d.data.count} (${Math.round(d.data.count / filteredData.length * 100)}%)`)
                    .attr('fill', d.data.color);
            })
            .on('mouseout', function() {
                hideTooltip();
                
                // Reset center text
                svg.select('.donut-center-text')
                    .text(`${filteredData.length} tools`)
                    .attr('fill', '#212529');
            })
            .transition()
            .duration(500)
            .attrTween('d', function(d) {
                const interpolate = d3.interpolate({ startAngle: d.startAngle, endAngle: d.startAngle }, d);
                return function(t) {
                    return arc(interpolate(t));
                };
            });
        
        // Update existing segments
        segments.transition()
            .duration(500)
            .attrTween('d', function(d) {
                const current = this._current || { startAngle: 0, endAngle: 0 };
                const interpolate = d3.interpolate(current, d);
                this._current = interpolate(1);
                return function(t) {
                    return arc(interpolate(t));
                };
            })
            .attr('fill', d => d.data.color);
        
        // Update center text
        svg.select('.donut-center-text')
            .text(`${filteredData.length} tools`);
        
        // Update legend
        const legend = svg.select('.donut-legend')
            .selectAll('.legend-item')
            .data(data);
        
        // Remove old legend items
        legend.exit().remove();
        
        // Create legend item groups
        const legendEnter = legend.enter()
            .append('g')
            .attr('class', 'legend-item')
            .attr('transform', (d, i) => `translate(0,${i * 20})`);
        
        // Add legend color boxes
        legendEnter.append('rect')
            .attr('width', 12)
            .attr('height', 12)
            .attr('fill', d => d.color);
        
        // Add legend text
        legendEnter.append('text')
            .attr('x', 20)
            .attr('y', 10)
            .style('font-size', '12px')
            .text(d => `${d.status} (${d.count})`);
        
        // Update existing legend items
        legend.select('rect')
            .attr('fill', d => d.color);
        
        legend.select('text')
            .text(d => `${d.status} (${d.count})`);
        
        // Position legend
        svg.select('.donut-legend')
            .attr('transform', `translate(${-(width / 2 - 20)},${-(height / 4)})`);
    }
    
    // Create language usage chart
    function createLanguageChart() {
        const chartContainer = document.getElementById('language-chart');
        const width = chartContainer.clientWidth;
        const height = chartContainer.clientHeight;
        const margin = { top: 30, right: 30, bottom: 40, left: 100 };
        
        // Create SVG
        const svg = d3.select(chartContainer)
            .append('svg')
            .attr('width', width)
            .attr('height', height)
            .append('g')
            .attr('transform', `translate(${margin.left},${margin.top})`);
        
        // Create scales (will be updated in updateLanguageChart)
        const x = d3.scaleLinear()
            .range([0, width - margin.left - margin.right]);
        
        const y = d3.scaleBand()
            .range([0, height - margin.top - margin.bottom])
            .padding(0.2);
        
        // Create axes
        svg.append('g')
            .attr('class', 'x-axis')
            .attr('transform', `translate(0,${height - margin.top - margin.bottom})`);
        
        svg.append('g')
            .attr('class', 'y-axis');
        
        // Initial update
        updateLanguageChart();
    }
    
    // Update language usage chart
    function updateLanguageChart() {
        const chartContainer = document.getElementById('language-chart');
        const width = chartContainer.clientWidth;
        const height = chartContainer.clientHeight;
        const margin = { top: 30, right: 30, bottom: 40, left: 100 };
        
        // Calculate language counts
        const languageCounts = {};
        filteredData.forEach(tool => {
            if (tool.language) {
                languageCounts[tool.language] = (languageCounts[tool.language] || 0) + 1;
            }
        });
        
        // Convert to array and sort
        const data = Object.entries(languageCounts)
            .map(([language, count]) => ({ language, count }))
            .sort((a, b) => b.count - a.count)
            .slice(0, 10); // Show top 10 languages
        
        // Update scales
        const x = d3.scaleLinear()
            .domain([0, d3.max(data, d => d.count)])
            .nice()
            .range([0, width - margin.left - margin.right]);
        
        const y = d3.scaleBand()
            .domain(data.map(d => d.language))
            .range([0, height - margin.top - margin.bottom])
            .padding(0.2);
        
        // Color scale for languages
        const languageColor = d3.scaleOrdinal()
            .domain(data.map(d => d.language))
            .range(d3.schemeTableau10);
        
        // Update axes
        const svg = d3.select(chartContainer).select('svg').select('g');
        
        svg.select('.x-axis')
            .transition()
            .duration(500)
            .call(d3.axisBottom(x)
                .ticks(5)
                .tickSizeOuter(0));
        
        svg.select('.y-axis')
            .transition()
            .duration(500)
            .call(d3.axisLeft(y)
                .tickSizeOuter(0));
        
        // Update bars
        const bars = svg.selectAll('.language-bar')
            .data(data, d => d.language);
        
        // Remove old bars
        bars.exit()
            .transition()
            .duration(300)
            .attr('width', 0)
            .remove();
        
        // Add new bars
        bars.enter()
            .append('rect')
            .attr('class', 'language-bar bar')
            .attr('x', 0)
            .attr('y', d => y(d.language))
            .attr('height', y.bandwidth())
            .attr('width', 0)
            .attr('fill', d => languageColor(d.language))
            .on('mouseover', function(event, d) {
                showTooltip(event, `<strong>${d.language}</strong><br>${d.count} tools (${Math.round(d.count / filteredData.length * 100)}%)`);
            })
            .on('mouseout', hideTooltip)
            .transition()
            .duration(500)
            .attr('width', d => x(d.count));
        
        // Update existing bars
        bars.transition()
            .duration(500)
            .attr('y', d => y(d.language))
            .attr('height', y.bandwidth())
            .attr('width', d => x(d.count))
            .attr('fill', d => languageColor(d.language));
        
        // Add value labels
        const labels = svg.selectAll('.language-label')
            .data(data, d => d.language);
        
        // Remove old labels
        labels.exit().remove();
        
        // Add new labels
        labels.enter()
            .append('text')
            .attr('class', 'language-label')
            .attr('x', d => x(d.count) + 5)
            .attr('y', d => y(d.language) + y.bandwidth() / 2)
            .attr('dy', '0.35em')
            .style('font-size', '12px')
            .text(d => d.count)
            .style('opacity', 0)
            .transition()
            .duration(500)
            .style('opacity', 1);
        
        // Update existing labels
        labels.transition()
            .duration(500)
            .attr('x', d => x(d.count) + 5)
            .attr('y', d => y(d.language) + y.bandwidth() / 2)
            .text(d => d.count);
    }
    
    // Create stars distribution chart
    function createStarsChart() {
        const chartContainer = document.getElementById('stars-chart');
        const width = chartContainer.clientWidth;
        const height = chartContainer.clientHeight;
        const margin = { top: 30, right: 30, bottom: 40, left: 60 };
        
        // Create SVG
        const svg = d3.select(chartContainer)
            .append('svg')
            .attr('width', width)
            .attr('height', height)
            .append('g')
            .attr('transform', `translate(${margin.left},${margin.top})`);
        
        // Create scales (will be updated in updateStarsChart)
        const x = d3.scaleLinear()
            .range([0, width - margin.left - margin.right]);
        
        const y = d3.scaleLinear()
            .range([height - margin.top - margin.bottom, 0]);
        
        // Create axes
        svg.append('g')
            .attr('class', 'x-axis')
            .attr('transform', `translate(0,${height - margin.top - margin.bottom})`);
        
        svg.append('g')
            .attr('class', 'y-axis');
        
        // Initial update
        updateStarsChart();
    }
    
    // Update stars distribution chart
    function updateStarsChart() {
        const chartContainer = document.getElementById('stars-chart');
        const width = chartContainer.clientWidth;
        const height = chartContainer.clientHeight;
        const margin = { top: 30, right: 30, bottom: 40, left: 60 };
        
        // Get all star values
        const stars = filteredData.map(tool => tool.stars || 0);
        
        // Create histogram bins
        const maxStars = d3.max(stars) || 100;
        const binCount = 10;
        const binWidth = Math.ceil(maxStars / binCount);
        
        // Define histogram function
        const histogram = d3.histogram()
            .value(d => d)
            .domain([0, maxStars])
            .thresholds(d3.range(0, maxStars, binWidth));
        
        // Generate bins
        const bins = histogram(stars);
        
        // Update scales
        const x = d3.scaleLinear()
            .domain([0, maxStars])
            .range([0, width - margin.left - margin.right]);
        
        const y = d3.scaleLinear()
            .domain([0, d3.max(bins, d => d.length)])
            .nice()
            .range([height - margin.top - margin.bottom, 0]);
        
        // Update axes
        const svg = d3.select(chartContainer).select('svg').select('g');
        
        svg.select('.x-axis')
            .transition()
            .duration(500)
            .call(d3.axisBottom(x)
                .ticks(binCount)
                .tickFormat(d => d)
                .tickSizeOuter(0));
        
        svg.select('.y-axis')
            .transition()
            .duration(500)
            .call(d3.axisLeft(y)
                .ticks(5)
                .tickSizeOuter(0));
        
        // Update bars
        const bars = svg.selectAll('.stars-bar')
            .data(bins);
        
        // Remove old bars
        bars.exit()
            .transition()
            .duration(300)
            .attr('y', height - margin.top - margin.bottom)
            .attr('height', 0)
            .remove();
        
        // Add new bars
        bars.enter()
            .append('rect')
            .attr('class', 'stars-bar bar')
            .attr('x', d => x(d.x0))
            .attr('width', d => Math.max(0, x(d.x1) - x(d.x0) - 1))
            .attr('y', height - margin.top - margin.bottom)
            .attr('height', 0)
            .attr('fill', '#6610f2')
            .on('mouseover', function(event, d) {
                showTooltip(event, `<strong>${d.x0} - ${d.x1} stars</strong><br>${d.length} tools`);
            })
            .on('mouseout', hideTooltip)
            .transition()
            .duration(500)
            .attr('y', d => y(d.length))
            .attr('height', d => height - margin.top - margin.bottom - y(d.length));
        
        // Update existing bars
        bars.transition()
            .duration(500)
            .attr('x', d => x(d.x0))
            .attr('width', d => Math.max(0, x(d.x1) - x(d.x0) - 1))
            .attr('y', d => y(d.length))
            .attr('height', d => height - margin.top - margin.bottom - y(d.length));
    }
    
    // Create timeline chart
    function createTimelineChart() {
        const chartContainer = document.getElementById('timeline-chart');
        const width = chartContainer.clientWidth;
        const height = chartContainer.clientHeight;
        const margin = { top: 30, right: 30, bottom: 40, left: 60 };
        
        // Create SVG
        const svg = d3.select(chartContainer)
            .append('svg')
            .attr('width', width)
            .attr('height', height)
            .append('g')
            .attr('transform', `translate(${margin.left},${margin.top})`);
        
        // Create scales (will be updated in updateTimelineChart)
        const x = d3.scaleLinear()
            .range([0, width - margin.left - margin.right]);
        
        const y = d3.scaleLinear()
            .range([height - margin.top - margin.bottom, 0]);
        
        // Create axes
        svg.append('g')
            .attr('class', 'x-axis')
            .attr('transform', `translate(0,${height - margin.top - margin.bottom})`);
        
        svg.append('g')
            .attr('class', 'y-axis');
        
        // Initial update
        updateTimelineChart();
    }
    
    // Update timeline chart
    function updateTimelineChart() {
        const chartContainer = document.getElementById('timeline-chart');
        const width = chartContainer.clientWidth;
        const height = chartContainer.clientHeight;
        const margin = { top: 30, right: 30, bottom: 40, left: 60 };
        
        // Filter tools with valid update times
        const data = filteredData.filter(tool => 
            tool.updateTime !== null && tool.updateTime !== undefined
        );
        
        // Group tools by update time
        const updateGroups = d3.group(data, d => Math.floor(d.updateTime));
        
        // Create data points
        const timelineData = Array.from(updateGroups, ([month, tools]) => ({
            month,
            count: tools.length,
            tools
        })).sort((a, b) => a.month - b.month);
        
        // Update scales
        const x = d3.scaleLinear()
            .domain([0, d3.max(timelineData, d => d.month) || 24])
            .range([0, width - margin.left - margin.right]);
        
        const y = d3.scaleLinear()
            .domain([0, d3.max(timelineData, d => d.count) || 1])
            .nice()
            .range([height - margin.top - margin.bottom, 0]);
        
        // Update axes
        const svg = d3.select(chartContainer).select('svg').select('g');
        
        svg.select('.x-axis')
            .transition()
            .duration(500)
            .call(d3.axisBottom(x)
                .ticks(10)
                .tickFormat(d => `${d} mo ago`)
                .tickSizeOuter(0));
        
        svg.select('.y-axis')
            .transition()
            .duration(500)
            .call(d3.axisLeft(y)
                .ticks(5)
                .tickSizeOuter(0));
        
        // Create line generator
        const line = d3.line()
            .x(d => x(d.month))
            .y(d => y(d.count))
            .curve(d3.curveMonotoneX);
        
        // Update line
        const path = svg.selectAll('.timeline-line')
            .data([timelineData]);
        
        if (path.empty()) {
            svg.append('path')
                .attr('class', 'timeline-line')
                .attr('fill', 'none')
                .attr('stroke', '#0d6efd')
                .attr('stroke-width', 2)
                .attr('d', line);
        } else {
            path.transition()
                .duration(500)
                .attr('d', line);
        }
        
        // Update points
        const points = svg.selectAll('.timeline-point')
            .data(timelineData);
        
        // Remove old points
        points.exit()
            .transition()
            .duration(300)
            .attr('r', 0)
            .remove();
        
        // Add new points
        points.enter()
            .append('circle')
            .attr('class', 'timeline-point')
            .attr('cx', d => x(d.month))
            .attr('cy', d => y(d.count))
            .attr('r', 0)
            .attr('fill', d => {
                if (d.month <= 6) return statusColors.recent;
                if (d.month <= 12) return statusColors.moderate;
                return statusColors.stale;
            })
            .on('mouseover', function(event, d) {
                const tooltipContent = `
                    <strong>${d.month} months ago</strong><br>
                    ${d.count} tools updated<br>
                    <small>Click to see tools</small>
                `;
                showTooltip(event, tooltipContent);
            })
            .on('mouseout', hideTooltip)
            .on('click', function(event, d) {
                // Show tools updated in this month
                showToolsForMonth(d.month, d.tools);
            })
            .transition()
            .duration(500)
            .attr('r', 5);
        
        // Update existing points
        points.transition()
            .duration(500)
            .attr('cx', d => x(d.month))
            .attr('cy', d => y(d.count))
            .attr('fill', d => {
                if (d.month <= 6) return statusColors.recent;
                if (d.month <= 12) return statusColors.moderate;
                return statusColors.stale;
            });
        
        // Add axis labels
        if (!svg.select('.x-label').size()) {
            svg.append('text')
                .attr('class', 'x-label')
                .attr('text-anchor', 'middle')
                .attr('x', (width - margin.left - margin.right) / 2)
                .attr('y', height - margin.top)
                .style('font-size', '12px')
                .text('Months since last update');
        }
        
        if (!svg.select('.y-label').size()) {
            svg.append('text')
                .attr('class', 'y-label')
                .attr('text-anchor', 'middle')
                .attr('transform', `translate(${-margin.left / 2},${(height - margin.top - margin.bottom) / 2}) rotate(-90)`)
                .style('font-size', '12px')
                .text('Number of tools');
        }
    }
    
    // Show tools updated in a specific month
    function showToolsForMonth(month, tools) {
        // Create modal to show the tools
        const modal = document.createElement('div');
        modal.className = 'modal fade';
        modal.id = 'timelineModal';
        modal.tabIndex = '-1';
        modal.setAttribute('aria-labelledby', 'timelineModalLabel');
        modal.setAttribute('aria-hidden', 'true');
        
        modal.innerHTML = `
            <div class="modal-dialog modal-lg">
                <div class="modal-content">
                    <div class="modal-header">
                        <h5 class="modal-title" id="timelineModalLabel">Tools Updated ${month} Months Ago</h5>
                        <button type="button" class="btn-close" data-bs-dismiss="modal" aria-label="Close"></button>
                    </div>
                    <div class="modal-body">
                        <div class="list-group">
                            ${tools.map(tool => `
                                <a href="${tool.url}" target="_blank" class="list-group-item list-group-item-action">
                                    <div class="d-flex justify-content-between align-items-center">
                                        <h5 class="mb-1">${tool.name}</h5>
                                        <span class="badge bg-primary rounded-pill">${tool.stars || 0} ★</span>
                                    </div>
                                    <p class="mb-1">${tool.description || ''}</p>
                                    <small>${tool.category}${tool.subcategory ? ' › ' + tool.subcategory : ''} · ${tool.language || 'Unknown language'}</small>
                                </a>
                            `).join('')}
                        </div>
                    </div>
                    <div class="modal-footer">
                        <button type="button" class="btn btn-secondary" data-bs-dismiss="modal">Close</button>
                    </div>
                </div>
            </div>
        `;
        
        // Remove existing modal if any
        const existingModal = document.getElementById('timelineModal');
        if (existingModal) {
            existingModal.remove();
        }
        
        // Add modal to the document
        document.body.appendChild(modal);
        
        // Initialize and show the modal
        const bsModal = new bootstrap.Modal(modal);
        bsModal.show();
    }
    
    // Update top tools section
    function updateTopTools() {
        // Get top 10 tools by stars
        const topTools = [...filteredData]
            .sort((a, b) => (b.stars || 0) - (a.stars || 0))
            .slice(0, 10);
        
        // Update the list
        const topToolsList = document.getElementById('top-tools');
        topToolsList.innerHTML = '';
        
        topTools.forEach(tool => {
            const listItem = document.createElement('a');
            listItem.href = '#';
            listItem.className = 'list-group-item list-group-item-action';
            listItem.dataset.id = tool.id;
            
            // Get update status class
            let statusClass = 'status-unknown';
            if (tool.updateTime !== null && tool.updateTime !== undefined) {
                if (tool.updateTime <= 6) {
                    statusClass = 'status-recent';
                } else if (tool.updateTime <= 12) {
                    statusClass = 'status-moderate';
                } else {
                    statusClass = 'status-stale';
                }
            }
            
            listItem.innerHTML = `
                <div class="d-flex justify-content-between align-items-center">
                    <div>
                        <span class="me-2 ${statusClass}">●</span>
                        ${tool.name}
                    </div>
                    <span class="tool-stars">${tool.stars || 0} ★</span>
                </div>
                <small>${tool.category}${tool.subcategory ? ' › ' + tool.subcategory : ''} · ${tool.language || 'Unknown language'}</small>
            `;
            
            // Add click event to show details
            listItem.addEventListener('click', function(event) {
                event.preventDefault();
                showToolDetails(tool);
                
                // Highlight selected tool
                document.querySelectorAll('#top-tools .list-group-item').forEach(item => {
                    item.classList.remove('active');
                });
                this.classList.add('active');
            });
            
            topToolsList.appendChild(listItem);
        });
    }
    
    // Update tools table
    function updateToolsTable() {
        const tableBody = document.querySelector('#tools-table tbody');
        tableBody.innerHTML = '';
        
        // Paginate if there are many tools
        const pageSize = 50;
        const pageCount = Math.ceil(filteredData.length / pageSize);
        const currentPage = 1;
        
        // Get tools for current page
        const toolsToShow = filteredData.slice(0, pageSize);
        
        // Add rows to the table
        toolsToShow.forEach(tool => {
            const row = document.createElement('tr');
            
            // Get update status class
            let statusClass = 'status-unknown';
            let updateText = 'Unknown';
            
            if (tool.updateTime !== null && tool.updateTime !== undefined) {
                updateText = `${tool.updateTime} months ago`;
                
                if (tool.updateTime <= 6) {
                    statusClass = 'status-recent';
                } else if (tool.updateTime <= 12) {
                    statusClass = 'status-moderate';
                } else {
                    statusClass = 'status-stale';
                }
            }
            
            row.innerHTML = `
                <td><a href="${tool.url}" target="_blank">${tool.name}</a></td>
                <td>${tool.category}${tool.subcategory ? ' › ' + tool.subcategory : ''}</td>
                <td>${tool.language || 'Unknown'}</td>
                <td>${tool.stars || 0}</td>
                <td class="${statusClass}">${updateText}</td>
                <td>
                    <button class="btn btn-sm btn-outline-primary view-details" data-id="${tool.id}">Details</button>
                </td>
            `;
            
            tableBody.appendChild(row);
        });
        
        // Add event listeners to view details buttons
        document.querySelectorAll('#tools-table .view-details').forEach(button => {
            button.addEventListener('click', function() {
                const toolId = this.dataset.id;
                const tool = filteredData.find(t => t.id === toolId);
                if (tool) {
                    showToolDetails(tool);
                }
            });
        });
    }
    
    // Show tool details
    function showToolDetails(tool) {
        selectedTool = tool;
        
        // Get update status
        let statusClass = 'status-unknown';
        let updateText = 'Unknown';
        
        if (tool.updateTime !== null && tool.updateTime !== undefined) {
            updateText = `${tool.updateTime} months ago`;
            
            if (tool.updateTime <= 6) {
                statusClass = 'status-recent';
            } else if (tool.updateTime <= 12) {
                statusClass = 'status-moderate';
            } else {
                statusClass = 'status-stale';
            }
        }
        
        // Update details pane
        const detailsPane = document.getElementById('tool-details');
        
        detailsPane.innerHTML = `
            <h4>${tool.name}</h4>
            <p>${tool.description || 'No description available.'}</p>
            
            <div class="detail-row">
                <div class="detail-label">Category:</div>
                <div class="detail-value">${tool.category}${tool.subcategory ? ' › ' + tool.subcategory : ''}</div>
            </div>
            
            <div class="detail-row">
                <div class="detail-label">Language:</div>
                <div class="detail-value">${tool.language || 'Unknown'}</div>
            </div>
            
            <div class="detail-row">
                <div class="detail-label">Stars:</div>
                <div class="detail-value">${tool.stars || 0} ★</div>
            </div>
            
            <div class="detail-row">
                <div class="detail-label">Last Updated:</div>
                <div class="detail-value ${statusClass}">${updateText}</div>
            </div>
            
            <div class="detail-row">
                <div class="detail-label">Repository:</div>
                <div class="detail-value">
                    <a href="${tool.url}" target="_blank">${tool.url}</a>
                </div>
            </div>
            
            <div class="mt-4">
                <a href="${tool.url}" target="_blank" class="btn btn-primary">
                    <i class="bi bi-github me-2"></i>View on GitHub
                </a>
            </div>
        `;
    }
    
    // Show tooltip
    function showTooltip(event, content) {
        // Create tooltip if it doesn't exist
        let tooltip = document.querySelector('.d3-tooltip');
        
        if (!tooltip) {
            tooltip = document.createElement('div');
            tooltip.className = 'd3-tooltip tooltip';
            document.body.appendChild(tooltip);
        }
        
        // Update content and position
        tooltip.innerHTML = content;
        tooltip.style.left = `${event.pageX + 10}px`;
        tooltip.style.top = `${event.pageY - 10}px`;
        tooltip.style.opacity = 1;
    }
    
    // Hide tooltip
    function hideTooltip() {
        const tooltip = document.querySelector('.d3-tooltip');
        if (tooltip) {
            tooltip.style.opacity = 0;
        }
    }
    
    // Show dashboard content
    function showDashboard() {
        document.getElementById('loading-indicator').style.display = 'none';
        document.getElementById('dashboard-content').style.display = 'block';
    }
    
    // Show error message
    function showError() {
        document.getElementById('loading-indicator').style.display = 'none';
        document.getElementById('error-message').style.display = 'block';
    }
    
    // Helper function: Debounce
    function debounce(func, wait) {
        let timeout;
        return function() {
            const context = this;
            const args = arguments;
            clearTimeout(timeout);
            timeout = setTimeout(() => {
                func.apply(context, args);
            }, wait);
        };
    }
    
    // Handle window resize
    window.addEventListener('resize', debounce(() => {
        updateCharts();
    }, 200));
});