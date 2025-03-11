// This script parses the README.md file and generates a JSON data structure for the visualization
// It can be run with Node.js: node generate_data.js > data.json

const fs = require('fs');
const path = require('path');

// Read the README.md file
const readmePath = path.join(__dirname, 'README.md');
const repoUpdatesPath = path.join(__dirname, 'repo_updates.json'); 

// Check if both files exist
if (!fs.existsSync(readmePath)) {
    console.error('README.md file not found');
    process.exit(1);
}

let repoUpdates = {};
if (fs.existsSync(repoUpdatesPath)) {
    try {
        repoUpdates = JSON.parse(fs.readFileSync(repoUpdatesPath, 'utf8'));
    } catch (error) {
        console.error('Error parsing repo_updates.json:', error);
        // Continue without repo updates
    }
}

const readmeContent = fs.readFileSync(readmePath, 'utf8');

// Parse README.md to extract categories and tools
function parseReadme(content) {
    const lines = content.split('\n');
    
    const nodes = [];
    const links = [];
    const categories = [];
    const languages = new Set();
    
    let currentCategory = null;
    let currentSubcategory = null;
    
    lines.forEach(line => {
        // Match category headers (## Title)
        const categoryMatch = line.match(/^## (.+)$/);
        if (categoryMatch && !categoryMatch[1].includes('License') && 
            !categoryMatch[1].includes('Contributing') && 
            !categoryMatch[1].includes('Acknowledgments') && 
            !categoryMatch[1].includes('Maintenance') &&
            !categoryMatch[1].includes('Last Updated') &&
            !categoryMatch[1].includes('Table of Contents') &&
            !categoryMatch[1].includes('Getting Started') &&
            !categoryMatch[1].includes('Typical Workflows') &&
            !categoryMatch[1].includes('Introduction to Virome Analysis') &&
            !categoryMatch[1].includes('Popular Packages')) {
            
            currentCategory = categoryMatch[1].trim();
            currentSubcategory = null;
            
            // Add category node
            const categoryId = `category-${currentCategory}`;
            nodes.push({
                id: categoryId,
                name: currentCategory,
                type: 'category',
                size: 15,
                color: '#343a40'
            });
            
            categories.push(currentCategory);
        }
        
        // Match subcategory headers (### Title)
        const subcategoryMatch = line.match(/^### (.+)$/);
        if (subcategoryMatch && currentCategory) {
            currentSubcategory = subcategoryMatch[1].trim();
            
            // Add subcategory node
            const subcategoryId = `subcategory-${currentSubcategory}`;
            nodes.push({
                id: subcategoryId,
                name: currentSubcategory,
                type: 'subcategory',
                parent: currentCategory,
                size: 10,
                color: '#495057'
            });
            
            // Link subcategory to its parent category
            links.push({
                source: subcategoryId,
                target: `category-${currentCategory}`,
                value: 1
            });
        }
        
        // Match tool entries (- [ToolName](URL) - Description. [tag1] [tag2])
        const toolMatch = line.match(/^- \[(.+?)\]\((.+?)\)(?:\s+\[Updated: .+?\])*\s+-\s+(.+?)(?:\s+\[(.+?)\])?(?:\s+\[(.+?)\])?/);
        if (toolMatch && (currentCategory || currentSubcategory)) {
            const toolName = toolMatch[1].trim();
            const toolUrl = toolMatch[2].trim();
            const description = toolMatch[3].trim();
            
            // Extract package manager and language information if available
            const tags = [];
            if (toolMatch[4]) tags.push(toolMatch[4].trim());
            if (toolMatch[5]) tags.push(toolMatch[5].trim());
            
            // Determine language
            let language = 'Unknown';
            tags.forEach(tag => {
                // Check if tag matches a common programming language
                const langMatches = tag.match(/^(Python|Java|C\+\+|R|Perl|JavaScript|Ruby|Nextflow|Snakemake|Go|Julia|C#|PHP)$/i);
                if (langMatches) {
                    language = langMatches[1];
                }
            });
            
            languages.add(language);
            
            // Get update info from repo_updates.json if available
            let updateTime = null;
            let stars = 0;
            
            if (repoUpdates[toolUrl]) {
                if (repoUpdates[toolUrl].lastUpdate) {
                    // Convert date string to months ago
                    const lastUpdate = new Date(repoUpdates[toolUrl].lastUpdate);
                    const now = new Date();
                    const monthsAgo = (now.getFullYear() - lastUpdate.getFullYear()) * 12 + 
                                      (now.getMonth() - lastUpdate.getMonth());
                    updateTime = monthsAgo;
                }
                
                if (repoUpdates[toolUrl].stars) {
                    stars = repoUpdates[toolUrl].stars;
                }
            }
            
            // Determine color based on update time
            let color;
            if (updateTime === null) {
                color = '#6c757d'; // gray for unknown
            } else if (updateTime <= 6) {
                color = '#198754'; // green for recent
            } else if (updateTime <= 12) {
                color = '#ffc107'; // yellow for moderately recent
            } else {
                color = '#dc3545'; // red for old
            }
            
            // Calculate node size based on stars
            const size = Math.max(5, Math.min(15, 5 + stars / 20));
            
            // Add tool node
            const toolId = `tool-${toolName.replace(/[^a-zA-Z0-9]/g, '')}`;
            nodes.push({
                id: toolId,
                name: toolName,
                type: 'tool',
                category: currentCategory,
                subcategory: currentSubcategory,
                description: description,
                language: language,
                stars: stars,
                updateTime: updateTime,
                size: size,
                color: color,
                url: toolUrl
            });
            
            // Link tool to its category/subcategory
            if (currentSubcategory) {
                links.push({
                    source: toolId,
                    target: `subcategory-${currentSubcategory}`,
                    value: 1
                });
            } else {
                links.push({
                    source: toolId,
                    target: `category-${currentCategory}`,
                    value: 1
                });
            }
        }
    });
    
    return {
        nodes,
        links,
        categories,
        languages: Array.from(languages)
    };
}

// Add workflow connections between tools that are commonly used together
function addWorkflowConnections(data) {
    // Define common workflows based on mentioned combinations in the README
    // Each workflow is an array of tool names that are often used together
    const workflows = [
        ['CheckV', 'VirSorter2', 'VIBRANT'],
        ['iPHoP', 'CHERRY', 'vConTACT2'],
        ['Pharokka', 'DRAMv'],
        ['metaviralSPAdes', 'coronaSPAdes'],
        ['MARVEL', 'PPR-Meta', 'DeepVirFinder'],
        ['PhiSpy', 'PHASTER', 'Prophinder'],
        ['VirFinder', 'viralVerify', 'CheckV']
    ];
    
    // Create lookup table for tool IDs by name
    const toolLookup = {};
    data.nodes.forEach(node => {
        if (node.type === 'tool') {
            toolLookup[node.name] = node.id;
        }
    });
    
    // Add workflow connections
    workflows.forEach(workflow => {
        for (let i = 0; i < workflow.length; i++) {
            for (let j = i + 1; j < workflow.length; j++) {
                const source = toolLookup[workflow[i]];
                const target = toolLookup[workflow[j]];
                
                if (source && target) {
                    data.links.push({
                        source: source,
                        target: target,
                        value: 0.5,
                        isWorkflow: true
                    });
                }
            }
        }
    });
    
    return data;
}

// Parse README and generate the data structure
const data = parseReadme(readmeContent);
const dataWithWorkflows = addWorkflowConnections(data);

// Output the data as JSON
console.log(JSON.stringify(dataWithWorkflows, null, 2));