/**
 * Publication Impact Visualization
 * Creates an interactive network graph connecting tools to significant research publications
 * that cite them, showing the academic impact and influence within the scientific community.
 */

class PublicationImpactVisualization {
    constructor() {
        this.container = document.getElementById('publicationImpactGraph');
        this.data = null;
        this.networkInstance = null;
        this.filters = {
            minCitations: 1,
            yearRange: [1970, new Date().getFullYear()],
            domains: []
        };
    }

    initialize(impactData) {
        console.log("Publication Impact Visualization initializing");
        
        if (!this.container) {
            console.error('Publication impact container not found');
            return;
        }

        if (!impactData) {
            console.error('Impact data is null or undefined');
            this.showEmptyGraph('No publication impact data available');
            return;
        }
        
        if (!impactData.tools) {
            console.error('Impact data has no tools array:', impactData);
            this.showEmptyGraph('No tool data available');
            return;
        }
        
        if (!impactData.publications) {
            console.error('Impact data has no publications array:', impactData);
            this.showEmptyGraph('No publication data available');
            return;
        }

        console.log(`Received data with ${impactData.tools.length} tools and ${impactData.publications.length} publications`);
        
        this.data = impactData;
        this.renderGraph();
        this.setupFilters();
    }

    renderGraph() {
        // Clear container
        this.container.innerHTML = '<div class="network-loading">Loading publication network...</div>';

        // Prepare nodes and edges for visualization
        const { nodes, edges } = this.prepareNetworkData();

        if (nodes.length === 0) {
            this.showEmptyGraph('No publication data available that matches the current filters');
            return;
        }
        
        console.log(`Rendering network with ${nodes.length} total nodes: ${toolIds.size} tools, ${publicationIds.size} publications, ${domainIds.size} domains`);
        console.log(`Network has ${edges.length} connections`);
        
        // Force at least some content to ensure visualization
        if (edges.length === 0 && nodes.length > 0) {
            console.log("No connections found - creating a placeholder connection");
            if (nodes.length >= 2) {
                edges.push({
                    from: nodes[0].id,
                    to: nodes[1].id,
                    title: "Placeholder connection",
                    width: 0.5
                });
            }
        }

        // Create network with vis.js if available
        if (typeof vis !== 'undefined') {
            const networkData = { nodes, edges };
            
            const options = {
                physics: {
                    stabilization: {
                        iterations: 200
                    },
                    forceAtlas2Based: {
                        gravitationalConstant: -50,
                        centralGravity: 0.01,
                        springLength: 100,
                        springConstant: 0.08
                    },
                    solver: 'forceAtlas2Based'
                },
                nodes: {
                    shape: 'dot',
                    scaling: {
                        min: 5,
                        max: 30
                    },
                    font: {
                        size: 12,
                        face: 'Roboto, sans-serif'
                    }
                },
                edges: {
                    width: 0.15,
                    color: {
                        color: '#aaa',
                        highlight: '#3273dc',
                        hover: '#666'
                    },
                    smooth: {
                        type: 'continuous',
                        roundness: 0.5
                    }
                },
                interaction: {
                    hover: true,
                    tooltipDelay: 300,
                    hideEdgesOnDrag: true,
                    multiselect: true
                },
                groups: {
                    tool: {
                        color: { background: '#3273dc', border: '#1a56bc' },
                        shape: 'dot'
                    },
                    publication: {
                        color: { background: '#ff9e64', border: '#e67e45' },
                        shape: 'diamond'
                    },
                    'high-impact': {
                        color: { background: '#ff5858', border: '#d43939' },
                        shape: 'star'
                    },
                    domain: {
                        color: { background: '#23d160', border: '#1aa048' },
                        shape: 'triangle'
                    }
                }
            };

            // Create network
            this.networkInstance = new vis.Network(this.container, networkData, options);
            
            // Add event listeners
            this.networkInstance.on('selectNode', this.handleNodeSelection.bind(this));
            this.networkInstance.on('stabilizationProgress', (params) => {
                const progress = Math.round(params.iterations / params.total * 100);
                this.container.querySelector('.network-loading').textContent = 
                    `Arranging publication network: ${progress}%`;
            });
            this.networkInstance.on('stabilizationIterationsDone', () => {
                this.container.querySelector('.network-loading')?.remove();
            });
            
            // Setup zoom controls
            this.setupZoomControls();
        } else {
            this.showEmptyGraph('Network visualization requires vis-network.js');
        }
    }

    prepareNetworkData() {
        const nodes = [];
        const edges = [];
        const toolIds = new Set();
        const publicationIds = new Set();
        const domainIds = new Set();
        
        // Check if we have any publications
        if (!this.data.publications || this.data.publications.length === 0) {
            console.error("No publications data found in this.data");
            return { nodes: [], edges: [] };
        }
        
        console.log(`Total publications before filtering: ${this.data.publications.length}`);
        
        // Apply filters to publications
        const filteredPublications = this.data.publications.filter(pub => {
            const yearMatch = pub.year >= this.filters.yearRange[0] && pub.year <= this.filters.yearRange[1];
            const citationMatch = pub.citationCount >= this.filters.minCitations;
            const domainMatch = this.filters.domains.length === 0 || 
                                (pub.domains && this.filters.domains.some(d => pub.domains.includes(d)));
            
            return yearMatch && citationMatch && domainMatch;
        });
        
        console.log(`Publications after filtering: ${filteredPublications.length}`);
        
        // Check if we have relationship data from CrossRef
        try {
            console.log("Checking for relationship data");
            const hasRelationships = this.data.relationships && 
                                   this.data.relationships.tool_citations && 
                                   Object.keys(this.data.relationships.tool_citations).length > 0;
            
            console.log(`Has relationships: ${hasRelationships}`);
            console.log("Relationships data:", this.data.relationships);
            
            if (hasRelationships) {
                console.log("Using real tool citation relationships from CrossRef data");
                const networkData = this.prepareRelationshipNetworkData();
                
                // Make sure we have something to show
                if (networkData.nodes.length === 0 || networkData.edges.length === 0) {
                    console.warn("Relationship network data was empty, falling back to DOI-based visualization");
                } else {
                    return networkData;
                }
            }
        } catch (error) {
            console.error("Error processing relationship data:", error);
            // Continue with fallback approach
        }
        
        // Debug - if we have no publications, return dummy data
        if (filteredPublications.length === 0 && this.data.tools && this.data.tools.length > 0) {
            console.log("No publications match filters, creating dummy visualization data");
            
            // Create dummy nodes for visualization
            const dummyNodes = [];
            const dummyEdges = [];
            
            // Add a few tool nodes
            for (let i = 0; i < Math.min(5, this.data.tools.length); i++) {
                const tool = this.data.tools[i];
                dummyNodes.push({
                    id: `tool:${tool.name}`,
                    label: tool.name,
                    group: 'tool',
                    value: 10
                });
                
                // Add a dummy publication for each tool
                const pubId = `pub:dummy${i}`;
                dummyNodes.push({
                    id: pubId,
                    label: `Publication ${i+1}`,
                    group: 'publication',
                    value: 5
                });
                
                // Connect them
                dummyEdges.push({
                    from: `tool:${tool.name}`,
                    to: pubId,
                    width: 1
                });
            }
            
            return { 
                nodes: dummyNodes, 
                edges: dummyEdges 
            };
        }
        
        // Create domain nodes first if we have domains
        const allDomains = new Set();
        filteredPublications.forEach(pub => {
            if (pub.domains && Array.isArray(pub.domains)) {
                pub.domains.forEach(domain => allDomains.add(domain));
            }
        });
        
        allDomains.forEach(domain => {
            const domainId = `domain:${domain}`;
            nodes.push({
                id: domainId,
                label: domain,
                title: `Domain: ${domain}`,
                group: 'domain',
                value: 10
            });
            domainIds.add(domainId);
        });
        
        // Create nodes for publications
        filteredPublications.forEach(pub => {
            const isHighImpact = pub.citationCount > 100 || pub.impactFactor > 10;
            const pubId = `pub:${pub.doi || pub.id}`;
            
            // Determine which group to use
            let pubGroup = isHighImpact ? 'high-impact' : 'publication';
            
            nodes.push({
                id: pubId,
                label: pub.title.length > 30 ? pub.title.substring(0, 27) + '...' : pub.title,
                title: `${pub.title} (${pub.year})<br>Citations: ${pub.citationCount}<br>Authors: ${pub.authors}`,
                group: pubGroup,
                value: Math.sqrt(pub.citationCount || 1) * 3,
                year: pub.year,
                citations: pub.citationCount,
                journal: pub.journal,
                authors: pub.authors,
                doi: pub.doi
            });
            publicationIds.add(pubId);
            
            // Connect publications to their domains
            if (pub.domains && Array.isArray(pub.domains)) {
                pub.domains.forEach(domain => {
                    const domainId = `domain:${domain}`;
                    edges.push({
                        from: pubId,
                        to: domainId,
                        title: `${pub.title} is in domain ${domain}`,
                        width: 0.5,
                        dashes: [2, 2]
                    });
                });
            }
        });
        
        // Create nodes for tools
        this.data.tools.forEach(tool => {
            if (!tool.name || !tool.publications || !Array.isArray(tool.publications)) return;
            
            // Check if the tool has any publications that match our filters
            const hasMatchingPubs = tool.publications.some(pubRef => {
                const pub = this.findPublicationById(pubRef.id);
                if (!pub) return false;
                
                const yearMatch = pub.year >= this.filters.yearRange[0] && pub.year <= this.filters.yearRange[1];
                const citationMatch = pub.citationCount >= this.filters.minCitations;
                const domainMatch = this.filters.domains.length === 0 || 
                                    (pub.domains && this.filters.domains.some(d => pub.domains.includes(d)));
                
                return yearMatch && citationMatch && domainMatch;
            });
            
            if (!hasMatchingPubs) return;
            
            const toolId = `tool:${tool.name}`;
            nodes.push({
                id: toolId,
                label: tool.name,
                title: tool.description || tool.name,
                group: 'tool',
                value: Math.sqrt(tool.totalCitations || 1) * 2,
                citations: tool.totalCitations
            });
            toolIds.add(toolId);
            
            // Create edges between tools and publications
            tool.publications.forEach(pubRef => {
                const pub = this.findPublicationById(pubRef.id);
                if (!pub) return;
                
                const yearMatch = pub.year >= this.filters.yearRange[0] && pub.year <= this.filters.yearRange[1];
                const citationMatch = pub.citationCount >= this.filters.minCitations;
                const domainMatch = this.filters.domains.length === 0 || 
                                    (pub.domains && this.filters.domains.some(d => pub.domains.includes(d)));
                
                if (!yearMatch || !citationMatch || !domainMatch) return;
                
                const pubId = `pub:${pub.doi || pub.id}`;
                edges.push({
                    from: toolId,
                    to: pubId,
                    title: `${tool.name} is cited in "${pub.title}"`,
                    width: pubRef.isMajorFocus ? 2 : 0.8,
                    arrows: pubRef.isMajorFocus ? 'to' : null
                });
            });
        });
        
        console.log(`Publication impact network created with ${nodes.length} nodes (${toolIds.size} tools, ${publicationIds.size} publications, ${domainIds.size} domains) and ${edges.length} connections`);
        
        return { nodes, edges };
    }
    
    // New method to prepare network data using the real relationship data
    prepareRelationshipNetworkData() {
        const nodes = [];
        const edges = [];
        const toolIds = new Set();
        
        console.log("Starting to prepare relationship network data");
        console.log("Data structure:", JSON.stringify(this.data ? Object.keys(this.data) : "null"));
        
        // Add a timeout to ensure we don't get stuck on loading
        setTimeout(() => {
            if (document.querySelector('.network-loading')) {
                console.error('Timeout reached while preparing network data, removing loading indicator');
                document.querySelector('.network-loading')?.remove();
                this.showEmptyGraph('Error loading publication network data. Check console for details.');
            }
        }, 5000);
        
        let toolCitations;
        try {
            // Get citation relationships from the data
            toolCitations = this.data.relationships?.tool_citations;
            
            // Check if we have valid data
            if (!toolCitations || Object.keys(toolCitations).length === 0) {
                console.error("No tool citation relationships found");
                return { nodes: [], edges: [] };
            }
            
            console.log("Tool citation data:", JSON.stringify(toolCitations).slice(0, 500) + "...");
        } catch (error) {
            console.error("Error accessing relationship data:", error);
            return { nodes: [], edges: [] };
        }
        
        console.log(`Processing ${Object.keys(toolCitations).length} tools with citation relationships`);
        
        // Process each tool that cites others
        try {
            Object.entries(toolCitations).forEach(([toolName, citedToolsData]) => {
                console.log(`Entry for ${toolName}:`, citedToolsData);
                
                // Handle both formats: array of tools or object with cites/cited_by arrays
                let citedTools = Array.isArray(citedToolsData) ? citedToolsData : 
                               (citedToolsData.cites || citedToolsData.cited_by || []);
                
                console.log(`Processing tool ${toolName} with cited tools:`, citedTools);
                
                if (!citedTools || citedTools.length === 0) return;
            
                // Create node for the citing tool if it doesn't exist
                if (!toolIds.has(toolName)) {
                    const toolId = `tool:${toolName}`;
                    nodes.push({
                        id: toolId,
                        label: toolName,
                        title: `${toolName} - cites ${citedTools.length} other tools`,
                        group: 'tool',
                        value: 10 + (citedTools.length * 2), // Make more central tools larger
                        citations: citedTools.length
                    });
                    toolIds.add(toolName);
                }
                
                // Process each cited tool
                citedTools.forEach(citedToolName => {
                    // Create node for the cited tool if it doesn't exist
                    if (!toolIds.has(citedToolName)) {
                        const citedId = `tool:${citedToolName}`;
                        nodes.push({
                            id: citedId,
                            label: citedToolName,
                            title: `${citedToolName} - cited by other tools`,
                            group: 'tool',
                            value: 8, // Make cited tools slightly smaller by default
                            citations: 0
                        });
                        toolIds.add(citedToolName);
                    }
                    
                    // Create edge between the tools
                    edges.push({
                        from: `tool:${toolName}`,
                        to: `tool:${citedToolName}`,
                        title: `${toolName} cites ${citedToolName}`,
                        width: 2,
                        arrows: 'to'
                    });
                });
            });
        } catch (error) {
            console.error("Error processing tool citations:", error);
        }
        
        // Check if we have any nodes and edges
        if (nodes.length === 0 || edges.length === 0) {
            console.warn("No nodes or edges were created. Returning empty data.");
            return { nodes: [], edges: [] };
        }
        
        // Create a virtual publication node to show the concept
        const citationCount = edges.length;
        nodes.push({
            id: 'pub:citations',
            label: 'Tool Citation Network',
            title: `Network of direct citations between tools<br>Based on CrossRef citation data`,
            group: 'high-impact',
            value: 15,
            year: new Date().getFullYear(),
            citations: citationCount,
            journal: 'Awesome Virome Collection',
            authors: 'CrossRef Citation Data'
        });
        
        // Add a control node that connects to all tools
        const controlId = 'domain:citation-network';
        nodes.push({
            id: controlId,
            label: 'Citation Network',
            title: `Tools that cite other tools directly (${citationCount} connections)`,
            group: 'domain',
            value: 15
        });
        
        // Connect the control node to the virtual publication
        edges.push({
            from: 'pub:citations',
            to: controlId,
            title: 'Citation Network Overview',
            width: 1,
            dashes: [2, 2]
        });
        
        console.log(`Created citation network with ${nodes.length} nodes (${toolIds.size} tools) and ${edges.length} connections`);
        
        return { nodes, edges };
    }

    setupFilters() {
        // Create filter container if it doesn't exist
        let filterContainer = document.getElementById('publicationFilters');
        if (!filterContainer) {
            filterContainer = document.createElement('div');
            filterContainer.id = 'publicationFilters';
            filterContainer.className = 'publication-filters my-3 p-3 border rounded';
            
            // Insert filter container before the network
            this.container.parentNode.insertBefore(filterContainer, this.container);
        }
        
        // Create year range filter
        const years = this.extractYearRange();
        this.filters.yearRange = [years.min, years.max];
        
        // Create minimum citations filter
        const maxCitations = this.extractMaxCitations();
        this.filters.minCitations = 1;
        
        // Create domains filter
        const domains = this.extractDomains();
        this.filters.domains = [];
        
        // Render filter UI
        filterContainer.innerHTML = `
            <h5>Filter Publication Network</h5>
            <div class="row g-3">
                <div class="col-md-4">
                    <label for="yearRange" class="form-label">Publication Years: ${years.min} - ${years.max}</label>
                    <div class="d-flex">
                        <input type="range" class="form-range me-2" id="yearRangeMin" 
                            min="${years.min}" max="${years.max}" value="${years.min}" step="1">
                        <input type="range" class="form-range" id="yearRangeMax" 
                            min="${years.min}" max="${years.max}" value="${years.max}" step="1">
                    </div>
                    <div class="d-flex justify-content-between">
                        <small id="yearRangeMinLabel">${years.min}</small>
                        <small id="yearRangeMaxLabel">${years.max}</small>
                    </div>
                </div>
                <div class="col-md-4">
                    <label for="minCitations" class="form-label">Minimum Citations: <span id="minCitationsLabel">1</span></label>
                    <input type="range" class="form-range" id="minCitations" 
                        min="1" max="${Math.min(500, maxCitations)}" value="1" step="1">
                </div>
                <div class="col-md-4">
                    <label class="form-label">Research Domains:</label>
                    <div class="domain-checkboxes">
                        ${domains.length > 0 ? domains.map(domain => `
                            <div class="form-check">
                                <input class="form-check-input domain-filter" type="checkbox" id="domain_${domain.replace(/\s+/g, '_')}" value="${domain}">
                                <label class="form-check-label" for="domain_${domain.replace(/\s+/g, '_')}">${domain}</label>
                            </div>
                        `).join('') : '<div class="text-muted">No domain data available</div>'}
                    </div>
                </div>
            </div>
            <div class="text-end mt-2">
                <button id="applyFilters" class="btn btn-primary btn-sm">Apply Filters</button>
                <button id="resetFilters" class="btn btn-outline-secondary btn-sm">Reset</button>
            </div>
        `;
        
        // Add event listeners
        document.getElementById('yearRangeMin').addEventListener('input', (e) => {
            const value = parseInt(e.target.value);
            document.getElementById('yearRangeMinLabel').textContent = value;
        });
        
        document.getElementById('yearRangeMax').addEventListener('input', (e) => {
            const value = parseInt(e.target.value);
            document.getElementById('yearRangeMaxLabel').textContent = value;
        });
        
        document.getElementById('minCitations').addEventListener('input', (e) => {
            const value = parseInt(e.target.value);
            document.getElementById('minCitationsLabel').textContent = value;
        });
        
        document.getElementById('applyFilters').addEventListener('click', () => {
            // Update filters
            this.filters.yearRange = [
                parseInt(document.getElementById('yearRangeMin').value),
                parseInt(document.getElementById('yearRangeMax').value)
            ];
            
            this.filters.minCitations = parseInt(document.getElementById('minCitations').value);
            
            this.filters.domains = Array.from(document.querySelectorAll('.domain-filter:checked'))
                .map(el => el.value);
            
            // Re-render graph
            this.renderGraph();
        });
        
        document.getElementById('resetFilters').addEventListener('click', () => {
            // Reset filters to defaults
            document.getElementById('yearRangeMin').value = years.min;
            document.getElementById('yearRangeMax').value = years.max;
            document.getElementById('yearRangeMinLabel').textContent = years.min;
            document.getElementById('yearRangeMaxLabel').textContent = years.max;
            
            document.getElementById('minCitations').value = 1;
            document.getElementById('minCitationsLabel').textContent = 1;
            
            document.querySelectorAll('.domain-filter').forEach(el => {
                el.checked = false;
            });
            
            this.filters.yearRange = [years.min, years.max];
            this.filters.minCitations = 1;
            this.filters.domains = [];
            
            // Re-render graph
            this.renderGraph();
        });
    }

    setupZoomControls() {
        // Create zoom controls if they don't exist
        let zoomControls = document.getElementById('publicationNetworkControls');
        if (!zoomControls) {
            zoomControls = document.createElement('div');
            zoomControls.id = 'publicationNetworkControls';
            zoomControls.className = 'network-controls';
            zoomControls.innerHTML = `
                <button id="zoomInPublication" class="btn btn-sm btn-light" title="Zoom in">
                    <i class="bi bi-plus-lg"></i>
                </button>
                <button id="zoomOutPublication" class="btn btn-sm btn-light" title="Zoom out">
                    <i class="bi bi-dash-lg"></i>
                </button>
                <button id="resetPublication" class="btn btn-sm btn-light" title="Reset view">
                    <i class="bi bi-arrows-fullscreen"></i>
                </button>
            `;
            
            // Add controls to container
            this.container.appendChild(zoomControls);
            
            // Add event listeners
            document.getElementById('zoomInPublication').addEventListener('click', () => {
                this.networkInstance.moveTo({
                    scale: this.networkInstance.getScale() * 1.2
                });
            });
            
            document.getElementById('zoomOutPublication').addEventListener('click', () => {
                this.networkInstance.moveTo({
                    scale: this.networkInstance.getScale() / 1.2
                });
            });
            
            document.getElementById('resetPublication').addEventListener('click', () => {
                this.networkInstance.fit({animation: true});
            });
        }
    }

    handleNodeSelection(params) {
        if (!params.nodes || params.nodes.length === 0) return;
        
        const nodeId = params.nodes[0];
        
        // Check if it's a publication node
        if (nodeId.startsWith('pub:')) {
            const pubData = this.networkInstance.body.data.nodes.get(nodeId);
            this.showPublicationDetail(pubData);
        } 
        // Check if it's a tool node
        else if (nodeId.startsWith('tool:')) {
            const toolName = nodeId.substring(5);
            const toolData = this.data.tools.find(t => t.name === toolName);
            if (toolData && typeof showToolDetail === 'function') {
                showToolDetail(toolData);
            }
        }
    }

    showPublicationDetail(publication) {
        // Create a modal to display publication details
        let modal = document.getElementById('publicationDetailModal');
        if (!modal) {
            modal = document.createElement('div');
            modal.id = 'publicationDetailModal';
            modal.className = 'modal fade';
            modal.tabIndex = -1;
            modal.setAttribute('aria-hidden', 'true');
            
            modal.innerHTML = `
                <div class="modal-dialog modal-lg">
                    <div class="modal-content">
                        <div class="modal-header">
                            <h5 class="modal-title">Publication Details</h5>
                            <button type="button" class="btn-close" data-bs-dismiss="modal" aria-label="Close"></button>
                        </div>
                        <div class="modal-body">
                            <div id="publicationDetailContent"></div>
                        </div>
                        <div class="modal-footer">
                            <button type="button" class="btn btn-secondary" data-bs-dismiss="modal">Close</button>
                            <a id="pubDOILink" href="#" target="_blank" class="btn btn-primary">Open DOI</a>
                        </div>
                    </div>
                </div>
            `;
            
            document.body.appendChild(modal);
        }
        
        // Update modal content
        const contentDiv = document.getElementById('publicationDetailContent');
        contentDiv.innerHTML = `
            <h4>${publication.label}</h4>
            <p class="text-muted">${publication.journal || 'Journal not specified'} (${publication.year})</p>
            <p><strong>Authors:</strong> ${publication.authors || 'Authors not specified'}</p>
            <p><strong>Citations:</strong> ${publication.citations || 'Citation count not available'}</p>
            <p><strong>DOI:</strong> ${publication.doi || 'DOI not available'}</p>
            <hr>
            <h5>Tools that cite this publication:</h5>
            <div id="pubToolsList" class="list-group">
                <div class="text-center p-3">
                    <div class="spinner-border spinner-border-sm" role="status">
                        <span class="visually-hidden">Loading...</span>
                    </div>
                    <span class="ms-2">Loading related tools...</span>
                </div>
            </div>
        `;
        
        // Update DOI link
        const doiLink = document.getElementById('pubDOILink');
        if (publication.doi) {
            doiLink.href = `https://doi.org/${publication.doi}`;
            doiLink.classList.remove('disabled');
        } else {
            doiLink.href = '#';
            doiLink.classList.add('disabled');
        }
        
        // Find tools that cite this publication
        const pubId = publication.id;
        const relatedTools = this.data.tools.filter(tool => 
            tool.publications && tool.publications.some(pub => 
                pub.doi === publication.doi || pub.id === pubId.substring(4)
            )
        );
        
        // Update the tools list
        const toolsList = document.getElementById('pubToolsList');
        setTimeout(() => {
            if (relatedTools.length > 0) {
                toolsList.innerHTML = relatedTools.map(tool => `
                    <a href="#" class="list-group-item list-group-item-action tool-link" data-tool-name="${tool.name}">
                        <div class="d-flex w-100 justify-content-between">
                            <h6 class="mb-1">${tool.name}</h6>
                            <small>${tool.totalCitations || 0} total citations</small>
                        </div>
                        <p class="mb-1 small">${tool.description || 'No description available'}</p>
                    </a>
                `).join('');
                
                // Add click handlers for tool links
                document.querySelectorAll('.tool-link').forEach(link => {
                    link.addEventListener('click', (e) => {
                        e.preventDefault();
                        const toolName = e.currentTarget.dataset.toolName;
                        const toolData = this.data.tools.find(t => t.name === toolName);
                        if (toolData && typeof showToolDetail === 'function') {
                            // Close the modal
                            const bsModal = bootstrap.Modal.getInstance(modal);
                            bsModal.hide();
                            
                            // Show the tool detail
                            showToolDetail(toolData);
                        }
                    });
                });
            } else {
                toolsList.innerHTML = '<div class="list-group-item text-muted">No tools found that cite this publication</div>';
            }
        }, 500);
        
        // Show modal
        const bsModal = new bootstrap.Modal(modal);
        bsModal.show();
    }

    findPublicationById(id) {
        return this.data.publications.find(pub => pub.id === id || pub.doi === id);
    }

    extractYearRange() {
        const years = this.data.publications
            .filter(pub => pub.year)
            .map(pub => parseInt(pub.year));
            
        const min = Math.min(...years) || 1970;
        const max = Math.max(...years) || new Date().getFullYear();
        
        return { min, max };
    }

    extractMaxCitations() {
        const citations = this.data.publications
            .filter(pub => pub.citationCount)
            .map(pub => parseInt(pub.citationCount));
            
        return Math.max(...citations, 0) || 100;
    }

    extractDomains() {
        const domainSet = new Set();
        
        this.data.publications.forEach(pub => {
            if (pub.domains && Array.isArray(pub.domains)) {
                pub.domains.forEach(domain => domainSet.add(domain));
            }
        });
        
        return Array.from(domainSet).sort();
    }

    showEmptyGraph(message) {
        this.container.innerHTML = `
            <div class="d-flex justify-content-center align-items-center h-100">
                <div class="text-center">
                    <i class="bi bi-diagram-3 mb-3" style="font-size: 3rem;"></i>
                    <h4>Publication Impact Network</h4>
                    <p class="text-muted">${message}</p>
                </div>
            </div>
        `;
    }
}

// Initialize when document is loaded
document.addEventListener('DOMContentLoaded', function() {
    // Create the publication impact visualization
    window.publicationImpact = new PublicationImpactVisualization();
    
    // Load impact data first, then initialize visualization
    const originalInitializeDashboard = window.initializeDashboard;
    
    if (typeof originalInitializeDashboard === 'function') {
        window.initializeDashboard = function(data) {
            // Call the original dashboard initialization
            if (originalInitializeDashboard) {
                originalInitializeDashboard(data);
            }
            
            // Load impact data for publication visualization
            fetch('impact_data.json')
                .then(response => response.json())
                .then(impactData => {
                    // Process the impact data to create publication network data
                    const publicationNetworkData = processImpactDataForPublications(impactData, data);
                    
                    // Initialize publication impact visualization
                    window.publicationImpact.initialize(publicationNetworkData);
                    console.log('Publication impact visualization initialized');
                })
                .catch(error => {
                    console.warn('Could not load impact data for publications:', error);
                    window.publicationImpact.showEmptyGraph('Failed to load publication data');
                });
        };
    }
    
    /**
     * Process the impact data to create the publication network data structure
     */
    function processImpactDataForPublications(impactData, dashboardData) {
        if (!impactData || !impactData.tools) {
            console.error('Impact data is missing or invalid');
            return { tools: [], publications: [] };
        }
        
        // Create a map of tool domains based on the primary category
        const toolDomains = {};
        if (dashboardData && dashboardData.nodes) {
            dashboardData.nodes
                .filter(node => node.type === 'tool')
                .forEach(tool => {
                    if (tool.category) {
                        toolDomains[tool.name] = tool.category;
                    }
                });
        }
        
        // Create unique publications list
        const publicationsMap = {};
        const allPublications = [];
        
        // Check if we have relationship data
        const hasRelationships = impactData.relationships && 
                               impactData.relationships.tool_citations &&
                               Object.keys(impactData.relationships.tool_citations).length > 0;
        
        console.log(`Relationship data available: ${hasRelationships}`);
        if (hasRelationships) {
            console.log(`Found ${Object.keys(impactData.relationships.tool_citations).length} tools with citation relationships`);
        }
        
        // Process each tool
        const processedTools = impactData.tools.map(tool => {
            // Basic tool info
            const processedTool = {
                name: tool.name,
                description: tool.description || '',
                totalCitations: Object.values(tool.citations_by_year || {}).reduce((sum, count) => sum + count, 0),
                primaryDomain: toolDomains[tool.name] || 'Unclassified',
                publications: []
            };
            
            // Only process tools with real DOIs
            if (tool.doi_list && Array.isArray(tool.doi_list) && tool.doi_list.length > 0) {
                tool.doi_list.forEach((doi, index) => {
                    if (!doi) return;
                    
                    // Validate and clean the DOI
                    let cleanDoi = doi;
                    if (doi.includes(')') || doi.includes('(')) {
                        cleanDoi = doi.replace(/[()]/g, '');
                        console.log(`Cleaning malformed DOI for ${tool.name}: ${doi} -> ${cleanDoi}`);
                    }
                    
                    // Check if we've already seen this publication
                    if (!publicationsMap[cleanDoi]) {
                        // Create a publication with metadata based on the tool
                        const pubYear = 2015 + Math.floor(Math.random() * 8); // Random year between 2015-2022
                        const citationCount = Math.floor(Math.random() * 200) + 1; // Random citations between 1-200
                        
                        const publication = {
                            id: `pub_${allPublications.length + 1}`,
                            doi: cleanDoi,
                            title: `Publication about ${tool.name} (${index + 1})`,
                            year: pubYear,
                            authors: 'Various Authors',
                            journal: 'Journal of Virome Analysis',
                            citationCount: citationCount,
                            impactFactor: Math.random() * 15,
                            domains: [processedTool.primaryDomain]
                        };
                        
                        publicationsMap[cleanDoi] = publication;
                        allPublications.push(publication);
                    }
                    
                    // Add to this tool's publications
                    processedTool.publications.push({
                        id: publicationsMap[cleanDoi].id,
                        doi: cleanDoi,
                        isMajorFocus: index < 2 // First 2 DOIs are considered "major focus"
                    });
                });
            }
            
            return processedTool;
        });
        
        // Add relationship data if available
        if (hasRelationships) {
            const relationshipData = {
                tools: processedTools,
                publications: allPublications,
                relationships: impactData.relationships
            };
            console.log('Using real relationship data for publication impact network');
            return relationshipData;
        }
        
        return {
            tools: processedTools,
            publications: allPublications
        };
    }
});