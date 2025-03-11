// Main visualization logic for Awesome-Virome tool network
document.addEventListener('DOMContentLoaded', function() {
    // Initialize visualization parameters
    const width = document.getElementById('visualization').clientWidth;
    const height = document.getElementById('visualization').clientHeight;
    
    // Create SVG container
    const svg = d3.select('#visualization')
        .append('svg')
        .attr('width', width)
        .attr('height', height);
    
    // Add zoom functionality
    const zoom = d3.zoom()
        .scaleExtent([0.5, 5])
        .on('zoom', (event) => {
            g.attr('transform', event.transform);
        });
    
    svg.call(zoom);
    
    // Create a group element for the entire visualization
    const g = svg.append('g');
    
    // Initialize the tooltip
    const tooltip = d3.select('body')
        .append('div')
        .attr('class', 'tooltip')
        .style('opacity', 0);
    
    // Load data
    loadData();
    
    function loadData() {
        // Add a loading spinner
        d3.select('#visualization')
            .append('div')
            .attr('class', 'spinner-container')
            .html(`
                <div class="spinner-border text-primary" role="status">
                    <span class="visually-hidden">Loading...</span>
                </div>
            `);
        
        // Load data from data.json file
        d3.json('data.json')
            .then(data => {
                // Remove loading spinner
                d3.select('.spinner-container').remove();
                
                // Process and visualize the data
                processData(data);
            })
            .catch(error => {
                console.error('Error loading data:', error);
                // If we can't load the data.json file, fall back to sample data
                console.log('Falling back to sample data...');
                createSampleData()
                    .then(sampleData => {
                        d3.select('.spinner-container').remove();
                        processData(sampleData);
                    })
                    .catch(sampleError => {
                        d3.select('.spinner-container').html(`
                            <div class="alert alert-danger">
                                Error loading data. Please try again later.
                            </div>
                        `);
                    });
            });
    }
    
    function createSampleData() {
        // This function creates sample data for demonstration
        // In a real implementation, you would fetch this from repo_updates.json
        return new Promise((resolve) => {
            // Sample categories
            const categories = [
                'Virus and Phage Identification', 
                'Host Prediction', 
                'Genome Analysis', 
                'Taxonomy', 
                'Functional Analysis',
                'CRISPR Analysis',
                'Sequence Analysis',
                'Visualization and Infrastructure',
                'Structural Analysis Tools',
                'Antimicrobial Resistance Analysis'
            ];
            
            // Sample programming languages
            const languages = ['Python', 'R', 'Java', 'C++', 'Perl', 'Nextflow', 'JavaScript'];
            
            // Sample update times (in months ago)
            const updateTimes = [1, 2, 3, 5, 7, 10, 13, 18, 24];
            
            // Generate sample nodes
            const nodes = [];
            
            // Add category nodes
            categories.forEach(category => {
                nodes.push({
                    id: `category-${category}`,
                    name: category,
                    type: 'category',
                    size: 10,
                    color: '#343a40'
                });
            });
            
            // Add tool nodes
            const tools = [];
            const toolCount = 100; // Number of sample tools
            
            for (let i = 0; i < toolCount; i++) {
                const category = categories[Math.floor(Math.random() * categories.length)];
                const language = languages[Math.floor(Math.random() * languages.length)];
                const stars = Math.floor(Math.random() * 200);
                const updateTime = updateTimes[Math.floor(Math.random() * updateTimes.length)];
                
                // Determine color based on update time
                let color;
                if (updateTime <= 6) {
                    color = '#198754'; // green
                } else if (updateTime <= 12) {
                    color = '#ffc107'; // yellow
                } else {
                    color = '#dc3545'; // red
                }
                
                // Calculate node size based on stars
                const size = Math.max(5, Math.min(15, 5 + stars / 20));
                
                const tool = {
                    id: `tool-${i}`,
                    name: `Tool ${i}`,
                    type: 'tool',
                    category: category,
                    language: language,
                    description: `This is a sample description for Tool ${i}, which is a ${language} tool for ${category}.`,
                    stars: stars,
                    updateTime: updateTime,
                    size: size,
                    color: color,
                    url: `https://github.com/example/tool${i}`
                };
                
                tools.push(tool);
                nodes.push(tool);
            }
            
            // Generate links
            const links = [];
            
            // Connect tools to their categories
            tools.forEach(tool => {
                links.push({
                    source: tool.id,
                    target: `category-${tool.category}`,
                    value: 1
                });
            });
            
            // Connect some tools together to show workflows
            // Each tool has a small chance to be connected to another random tool
            tools.forEach(tool => {
                if (Math.random() > 0.7) {
                    const randomTool = tools[Math.floor(Math.random() * tools.length)];
                    if (tool.id !== randomTool.id) {
                        links.push({
                            source: tool.id,
                            target: randomTool.id,
                            value: 0.5,
                            isWorkflow: true
                        });
                    }
                }
            });
            
            // Return the sample data
            resolve({
                nodes: nodes,
                links: links,
                categories: categories,
                languages: languages
            });
        });
    }
    
    function processData(data) {
        // Populate filter dropdowns
        populateFilterDropdowns(data);
        
        // Create the visualization
        createVisualization(data);
        
        // Setup event listeners for filters
        setupFilterEventListeners(data);
    }
    
    function populateFilterDropdowns(data) {
        // Populate category filter
        const categorySelect = document.getElementById('category-filter');
        data.categories.forEach(category => {
            const option = document.createElement('option');
            option.value = category;
            option.textContent = category;
            categorySelect.appendChild(option);
        });
        
        // Populate language filter
        const languageSelect = document.getElementById('language-filter');
        data.languages.forEach(language => {
            const option = document.createElement('option');
            option.value = language;
            option.textContent = language;
            languageSelect.appendChild(option);
        });
    }
    
    let simulation; // Declare simulation variable in a wider scope
    
    function createVisualization(data) {
        // Clear previous visualization
        g.selectAll('*').remove();
        
        // Create a force simulation
        simulation = d3.forceSimulation(data.nodes)
            .force('link', d3.forceLink(data.links).id(d => d.id).distance(d => d.isWorkflow ? 100 : 150))
            .force('charge', d3.forceManyBody().strength(d => d.type === 'category' ? -300 : -100))
            .force('center', d3.forceCenter(width / 2, height / 2))
            .force('collide', d3.forceCollide().radius(d => d.size * 2));
        
        // Create links
        const link = g.append('g')
            .attr('class', 'links')
            .selectAll('line')
            .data(data.links)
            .enter()
            .append('line')
            .attr('class', 'link')
            .attr('stroke-width', d => d.isWorkflow ? 2 : 1)
            .attr('stroke-dasharray', d => d.isWorkflow ? '5,5' : '0');
        
        // Create nodes
        const node = g.append('g')
            .attr('class', 'nodes')
            .selectAll('circle')
            .data(data.nodes)
            .enter()
            .append('circle')
            .attr('class', d => d.type === 'category' ? 'node-category' : 'node')
            .attr('r', d => d.size)
            .attr('fill', d => d.color)
            .on('mouseover', handleMouseOver)
            .on('mouseout', handleMouseOut)
            .on('click', handleClick)
            .call(d3.drag()
                .on('start', dragstarted)
                .on('drag', dragged)
                .on('end', dragended));
        
        // Add node labels
        const labels = g.append('g')
            .attr('class', 'labels')
            .selectAll('text')
            .data(data.nodes)
            .enter()
            .append('text')
            .text(d => d.name)
            .attr('font-size', d => d.type === 'category' ? '12px' : '10px')
            .attr('text-anchor', 'middle')
            .attr('dy', d => d.type === 'category' ? -15 : -10);
        
        // Update positions on each tick of the simulation
        simulation.on('tick', () => {
            // Update link positions
            link
                .attr('x1', d => d.source.x)
                .attr('y1', d => d.source.y)
                .attr('x2', d => d.target.x)
                .attr('y2', d => d.target.y);
            
            // Update node positions
            node
                .attr('cx', d => d.x = Math.max(d.size, Math.min(width - d.size, d.x)))
                .attr('cy', d => d.y = Math.max(d.size, Math.min(height - d.size, d.y)));
            
            // Update label positions
            labels
                .attr('x', d => d.x)
                .attr('y', d => d.y);
        });
        
        // Zoom to fit content
        zoomToFit();
    }
    
    function handleMouseOver(event, d) {
        // Show tooltip with basic information
        tooltip.transition()
            .duration(200)
            .style('opacity', 0.9);
        
        const content = d.type === 'category' 
            ? `<h5>${d.name}</h5><p>Category</p>` 
            : `<h5>${d.name}</h5>
               <p><strong>Category:</strong> ${d.category}</p>
               <p><strong>Language:</strong> ${d.language}</p>
               <p><strong>Stars:</strong> ${d.stars}</p>
               <p><strong>Last Updated:</strong> ${d.updateTime} months ago</p>
               <p>Click for more details</p>`;
        
        tooltip.html(content)
            .style('left', (event.pageX + 10) + 'px')
            .style('top', (event.pageY - 28) + 'px');
    }
    
    function handleMouseOut() {
        // Hide tooltip
        tooltip.transition()
            .duration(500)
            .style('opacity', 0);
    }
    
    function handleClick(event, d) {
        // Only show details for tools, not categories
        if (d.type !== 'tool') return;
        
        // Display tool details in the details section
        const toolDetails = document.getElementById('tool-details');
        const toolDetailsCard = document.getElementById('tool-details-card');
        
        // Create detail content
        const detailContent = `
            <table class="table table-striped">
                <tbody>
                    <tr><td>Name</td><td>${d.name}</td></tr>
                    <tr><td>Category</td><td>${d.category}</td></tr>
                    <tr><td>Language</td><td>${d.language}</td></tr>
                    <tr><td>Description</td><td>${d.description}</td></tr>
                    <tr><td>GitHub Stars</td><td>${d.stars}</td></tr>
                    <tr><td>Last Updated</td><td>${d.updateTime} months ago</td></tr>
                    <tr><td>Repository</td><td><a href="${d.url}" target="_blank">${d.url}</a></td></tr>
                </tbody>
            </table>
        `;
        
        toolDetails.innerHTML = detailContent;
        toolDetailsCard.style.display = 'block';
        
        // Scroll to details
        toolDetailsCard.scrollIntoView({ behavior: 'smooth' });
    }
    
    function setupFilterEventListeners(data) {
        // Category filter
        document.getElementById('category-filter').addEventListener('change', function() {
            applyFilters(data);
        });
        
        // Language filter
        document.getElementById('language-filter').addEventListener('change', function() {
            applyFilters(data);
        });
        
        // Update filter
        document.getElementById('update-filter').addEventListener('change', function() {
            applyFilters(data);
        });
        
        // Stars filter
        document.getElementById('stars-filter').addEventListener('input', function() {
            document.getElementById('stars-value').textContent = `${this.value}+ stars`;
            applyFilters(data);
        });
        
        // Search input
        document.getElementById('search-input').addEventListener('input', function() {
            applyFilters(data);
        });
        
        // Reset filters button
        document.getElementById('reset-filters').addEventListener('click', function() {
            // Reset all filter inputs
            document.getElementById('category-filter').value = 'all';
            document.getElementById('language-filter').value = 'all';
            document.getElementById('update-filter').value = 'all';
            document.getElementById('stars-filter').value = 0;
            document.getElementById('stars-value').textContent = '0+ stars';
            document.getElementById('search-input').value = '';
            
            // Apply the reset filters
            applyFilters(data);
        });
    }
    
    function applyFilters(data) {
        // Get filter values
        const categoryFilter = document.getElementById('category-filter').value;
        const languageFilter = document.getElementById('language-filter').value;
        const updateFilter = document.getElementById('update-filter').value;
        const starsFilter = parseInt(document.getElementById('stars-filter').value);
        const searchTerm = document.getElementById('search-input').value.toLowerCase();
        
        // Filter nodes
        const filteredNodes = data.nodes.filter(node => {
            // Always include category nodes
            if (node.type === 'category') return true;
            
            // Apply category filter
            if (categoryFilter !== 'all' && node.category !== categoryFilter) return false;
            
            // Apply language filter
            if (languageFilter !== 'all' && node.language !== languageFilter) return false;
            
            // Apply update filter
            if (updateFilter !== 'all' && node.updateTime > parseInt(updateFilter)) return false;
            
            // Apply stars filter
            if (node.stars < starsFilter) return false;
            
            // Apply search filter
            if (searchTerm && !node.name.toLowerCase().includes(searchTerm) && 
                !node.description.toLowerCase().includes(searchTerm)) return false;
            
            return true;
        });
        
        // Get IDs of filtered nodes
        const filteredNodeIds = filteredNodes.map(node => node.id);
        
        // Filter links to only include connections between filtered nodes
        const filteredLinks = data.links.filter(link => {
            const sourceId = typeof link.source === 'object' ? link.source.id : link.source;
            const targetId = typeof link.target === 'object' ? link.target.id : link.target;
            
            return filteredNodeIds.includes(sourceId) && filteredNodeIds.includes(targetId);
        });
        
        // Update nodes and links in the visualization
        updateVisualization({
            nodes: filteredNodes,
            links: filteredLinks
        });
    }
    
    function updateVisualization(filteredData) {
        // Create a new simulation with filtered data
        simulation.stop();
        
        simulation = d3.forceSimulation(filteredData.nodes)
            .force('link', d3.forceLink(filteredData.links).id(d => d.id).distance(d => d.isWorkflow ? 100 : 150))
            .force('charge', d3.forceManyBody().strength(d => d.type === 'category' ? -300 : -100))
            .force('center', d3.forceCenter(width / 2, height / 2))
            .force('collide', d3.forceCollide().radius(d => d.size * 2));
        
        // Update links
        const link = g.select('.links')
            .selectAll('line')
            .data(filteredData.links, d => `${d.source.id || d.source}-${d.target.id || d.target}`);
        
        link.exit().remove();
        
        const linkEnter = link.enter()
            .append('line')
            .attr('class', 'link')
            .attr('stroke-width', d => d.isWorkflow ? 2 : 1)
            .attr('stroke-dasharray', d => d.isWorkflow ? '5,5' : '0');
        
        // Update nodes
        const node = g.select('.nodes')
            .selectAll('circle')
            .data(filteredData.nodes, d => d.id);
        
        node.exit().remove();
        
        const nodeEnter = node.enter()
            .append('circle')
            .attr('class', d => d.type === 'category' ? 'node-category' : 'node')
            .attr('r', d => d.size)
            .attr('fill', d => d.color)
            .on('mouseover', handleMouseOver)
            .on('mouseout', handleMouseOut)
            .on('click', handleClick)
            .call(d3.drag()
                .on('start', dragstarted)
                .on('drag', dragged)
                .on('end', dragended));
        
        // Update labels
        const labels = g.select('.labels')
            .selectAll('text')
            .data(filteredData.nodes, d => d.id);
        
        labels.exit().remove();
        
        const labelEnter = labels.enter()
            .append('text')
            .text(d => d.name)
            .attr('font-size', d => d.type === 'category' ? '12px' : '10px')
            .attr('text-anchor', 'middle')
            .attr('dy', d => d.type === 'category' ? -15 : -10);
        
        // Apply search highlighting
        const searchTerm = document.getElementById('search-input').value.toLowerCase();
        g.selectAll('.node')
            .classed('search-highlight', d => {
                if (!searchTerm) return false;
                return d.name.toLowerCase().includes(searchTerm) ||
                       d.description.toLowerCase().includes(searchTerm);
            });
        
        // Update positions on each tick of the simulation
        simulation.on('tick', () => {
            // Link positions
            g.selectAll('.links line')
                .attr('x1', d => d.source.x)
                .attr('y1', d => d.source.y)
                .attr('x2', d => d.target.x)
                .attr('y2', d => d.target.y);
            
            // Node positions
            g.selectAll('.nodes circle')
                .attr('cx', d => d.x = Math.max(d.size, Math.min(width - d.size, d.x)))
                .attr('cy', d => d.y = Math.max(d.size, Math.min(height - d.size, d.y)));
            
            // Label positions
            g.selectAll('.labels text')
                .attr('x', d => d.x)
                .attr('y', d => d.y);
        });
        
        // Zoom to fit content
        zoomToFit();
    }
    
    function dragstarted(event, d) {
        if (!event.active) simulation.alphaTarget(0.3).restart();
        d.fx = d.x;
        d.fy = d.y;
    }
    
    function dragged(event, d) {
        d.fx = event.x;
        d.fy = event.y;
    }
    
    function dragended(event, d) {
        if (!event.active) simulation.alphaTarget(0);
        d.fx = null;
        d.fy = null;
    }
    
    function zoomToFit() {
        const bounds = g.node().getBBox();
        const dx = bounds.width;
        const dy = bounds.height;
        const x = bounds.x + (bounds.width / 2);
        const y = bounds.y + (bounds.height / 2);
        
        const scale = 0.8 / Math.max(dx / width, dy / height);
        const translate = [width / 2 - scale * x, height / 2 - scale * y];
        
        svg.transition()
            .duration(500)
            .call(zoom.transform, d3.zoomIdentity
                .translate(translate[0], translate[1])
                .scale(scale));
    }
});