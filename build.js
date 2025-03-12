// Build script for GitHub Pages deployment
// Generates a self-contained index.html with embedded data
// Usage: node build.js

const fs = require('fs');
const path = require('path');
const { exec } = require('child_process');

// Paths
const DATA_PATH = './data.json';
const INDEX_PATH = './index.html';
const OUTPUT_PATH = './dist';
const OUTPUT_INDEX = path.join(OUTPUT_PATH, 'index.html');

console.log('🔨 Building Awesome-Virome Dashboard for GitHub Pages...');

// Create dist directory if it doesn't exist
if (!fs.existsSync(OUTPUT_PATH)) {
    fs.mkdirSync(OUTPUT_PATH, { recursive: true });
    console.log(`✅ Created output directory: ${OUTPUT_PATH}`);
}

// Function to generate data.json from README if it doesn't exist
function ensureDataExists() {
    return new Promise((resolve, reject) => {
        if (fs.existsSync(DATA_PATH)) {
            console.log('✅ data.json already exists');
            resolve();
            return;
        }
        
        console.log('⚠️ data.json not found, attempting to generate it...');
        
        // Check if generate_data.js exists
        if (!fs.existsSync('./generate_data.js')) {
            console.error('❌ generate_data.js not found, cannot generate data');
            reject(new Error('generate_data.js not found'));
            return;
        }
        
        // Run generate_data.js to create data.json
        exec('node generate_data.js > data.json', (error, stdout, stderr) => {
            if (error) {
                console.error(`❌ Error generating data.json: ${error.message}`);
                reject(error);
                return;
            }
            if (stderr) {
                console.warn(`⚠️ Warning while generating data.json: ${stderr}`);
            }
            console.log('✅ Generated data.json from README');
            resolve();
        });
    });
}

// Main build function
async function build() {
    try {
        // Ensure data.json exists
        await ensureDataExists();
        
        // Read the data.json file
        const dataRaw = fs.readFileSync(DATA_PATH, 'utf8');
        
        // Parse and then stringify the data to ensure proper JSON format
        let data;
        try {
            data = JSON.parse(dataRaw);
            console.log(`✅ Successfully parsed data.json (${Object.keys(data).length} keys found)`);
        } catch (error) {
            console.error(`❌ Error parsing data.json: ${error.message}`);
            throw new Error(`Invalid JSON in data.json: ${error.message}`);
        }
        
        // Read the index.html file
        let html = fs.readFileSync(INDEX_PATH, 'utf8');
        
        // Replace the data placeholder with properly JSON stringified data
        // Use a clean JSON string without any HTML entities or special chars
        const jsonString = JSON.stringify(data, null, 0);
        console.log(`✅ Data serialized for embedding (${jsonString.length} characters)`);
        
        // Directly replace the placeholder with the clean JSON
        html = html.replace('<!-- DATA_PLACEHOLDER -->', jsonString);
        
        // Verify that the replacement was successful
        if (html.includes('<!-- DATA_PLACEHOLDER -->')) {
            console.error('❌ Error: Failed to replace data placeholder in HTML');
            throw new Error('Data placeholder was not replaced');
        }
        
        // Double check that the replacement actually worked and the script contains valid JSON
        const jsonStart = html.indexOf('<script id="embedded-data" type="application/json">') + 
                         '<script id="embedded-data" type="application/json">'.length;
        const jsonEnd = html.indexOf('</script>', jsonStart);
        const extractedJson = html.substring(jsonStart, jsonEnd).trim();
        
        // Validate the JSON
        try {
            JSON.parse(extractedJson);
            console.log('✅ Successfully verified that the embedded JSON is valid');
        } catch (error) {
            console.error(`❌ Error: Embedded JSON is not valid: ${error.message}`);
            console.error(`JSON snippet: ${extractedJson.substring(0, 100)}...`);
            throw new Error(`Invalid JSON was embedded: ${error.message}`);
        }
        
        // Write the result to the output directory
        fs.writeFileSync(OUTPUT_INDEX, html);
        console.log(`✅ Created ${OUTPUT_INDEX} with embedded data`);
        
        // Verify the data was properly embedded
        const outputHtml = fs.readFileSync(OUTPUT_INDEX, 'utf8');
        if (!outputHtml.includes('"nodes":') && !outputHtml.includes('"links":') && 
            !outputHtml.includes('"categories":')) {
            console.warn('⚠️ Warning: Output HTML may not contain expected data structure');
        } else {
            console.log('✅ Data appears to be properly embedded in HTML');
        }
        
        // Copy other necessary files to the output directory
        const filesToCopy = ['styles.css', 'visualization.js', 'data.json'];
        
        filesToCopy.forEach(file => {
            if (fs.existsSync(file)) {
                fs.copyFileSync(file, path.join(OUTPUT_PATH, file));
                console.log(`✅ Copied ${file} to ${OUTPUT_PATH}`);
            } else {
                console.warn(`⚠️ Warning: ${file} not found, skipping`);
            }
        });
        
        // Add a fallback data loading option with direct link to data.json
        const fallbackScript = `
        <script>
        // Fallback script to handle data loading failures
        window.addEventListener('DOMContentLoaded', function() {
            // After 5 seconds, check if dashboard is still loading
            setTimeout(function() {
                var loadingIndicator = document.getElementById('loading-indicator');
                var dashboardContent = document.getElementById('dashboard-content');
                
                if (loadingIndicator && loadingIndicator.style.display !== 'none' && 
                    dashboardContent && dashboardContent.style.display === 'none') {
                    console.warn('Dashboard still loading after 5s, adding data load fallback button');
                    
                    // Create fallback button
                    var fallbackButton = document.createElement('button');
                    fallbackButton.textContent = 'Try Loading External Data';
                    fallbackButton.className = 'btn btn-warning mt-3';
                    fallbackButton.onclick = function() {
                        // Directly load data.json using fetch
                        fetch('./data.json')
                            .then(function(response) { return response.json(); })
                            .then(function(data) {
                                console.log('Successfully loaded data via fallback button');
                                // Find visualization.js script and execute loadData function
                                if (window.processData && window.showDashboard) {
                                    window.processData(data);
                                    window.showDashboard();
                                } else {
                                    console.error('Cannot find visualization functions');
                                    alert('Could not process data. Please check console for errors.');
                                }
                            })
                            .catch(function(error) {
                                console.error('Fallback data loading failed:', error);
                                alert('Failed to load data: ' + error.message);
                            });
                    };
                    
                    // Add to loading indicator
                    loadingIndicator.appendChild(document.createElement('div')).appendChild(fallbackButton);
                }
            }, 5000);
        });
        </script>
        `;
        
        // Add fallback script to output HTML
        const outputHtmlWithFallback = fs.readFileSync(OUTPUT_INDEX, 'utf8') + fallbackScript;
        fs.writeFileSync(OUTPUT_INDEX, outputHtmlWithFallback);
        console.log('✅ Added fallback data loading script to output HTML');
        
        console.log('✅ Build completed successfully!');
        console.log('');
        console.log('To deploy to GitHub Pages:');
        console.log('1. Copy all files from the dist/ directory to your gh-pages branch');
        console.log('2. Push the gh-pages branch to GitHub');
        console.log('3. Configure your repository to serve from the gh-pages branch');
        console.log('');
        console.log('You can also use the gh-pages npm package for easier deployment:');
        console.log('npm install -g gh-pages');
        console.log('gh-pages -d dist');
        
    } catch (error) {
        console.error(`❌ Build failed: ${error.message}`);
        process.exit(1);
    }
}

// Run the build
build();