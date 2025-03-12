# Deploying the Awesome-Virome Dashboard to GitHub Pages

This guide explains how to deploy the Awesome-Virome Dashboard to GitHub Pages so it can be publicly accessed online.

## Option 1: Manual Deployment

### Prerequisites
- Your awesome-virome repository cloned locally
- Node.js installed for running the build script

### Steps

1. **Build the dashboard files**

   ```bash
   # From the root of your awesome-virome repository
   node build.js
   ```
   
   This will create a `dist` directory with all necessary files for GitHub Pages.

2. **Create a gh-pages branch (if you don't have one already)**

   ```bash
   git checkout -b gh-pages
   git rm -rf .
   git clean -fdx
   ```

3. **Copy built files to the gh-pages branch**

   ```bash
   cp -r dist/* .
   git add .
   git commit -m "Deploy dashboard to GitHub Pages"
   git push origin gh-pages
   ```

4. **Configure GitHub Pages in your repository settings**

   - Go to your repository on GitHub
   - Click on "Settings" > "Pages" 
   - Select "gh-pages" as the source branch
   - Click "Save"

5. **Access your dashboard**

   After a few minutes, your dashboard will be available at:
   ```
   https://yourusername.github.io/awesome-virome/
   ```

## Option 2: Automated Deployment with gh-pages package

### Prerequisites
- Your awesome-virome repository cloned locally
- Node.js installed

### Steps

1. **Install the gh-pages package**

   ```bash
   # Install globally
   npm install -g gh-pages
   # Or install locally in your project
   npm install --save-dev gh-pages
   ```

2. **Build the dashboard files**

   ```bash
   node build.js
   ```

3. **Deploy to GitHub Pages with one command**

   ```bash
   gh-pages -d dist
   ```

4. **Configure GitHub Pages in your repository settings**

   - Same as in Option 1, step 4

## Option 3: GitHub Actions Automated Deployment

For the most automated solution, you can set up a GitHub Actions workflow:

1. **Create a GitHub Actions workflow file**

   Create a file in your repository at `.github/workflows/deploy-dashboard.yml` with the following content:

   ```yaml
   name: Deploy Dashboard

   on:
     push:
       branches:
         - main
       paths:
         - 'README.md'
         - 'index.html'
         - 'styles.css'
         - 'visualization.js'
         - 'generate_data.js'
         - 'build.js'

   jobs:
     build-and-deploy:
       runs-on: ubuntu-latest
       steps:
         - name: Checkout
           uses: actions/checkout@v3
           
         - name: Setup Node.js
           uses: actions/setup-node@v3
           with:
             node-version: '16'
             
         - name: Build
           run: |
             node generate_data.js > data.json
             node build.js
             
         - name: Deploy
           uses: JamesIves/github-pages-deploy-action@4.1.5
           with:
             branch: gh-pages
             folder: dist
   ```

2. **Commit and push this workflow file**

   ```bash
   git add .github/workflows/deploy-dashboard.yml
   git commit -m "Add GitHub Actions workflow for dashboard deployment"
   git push origin main
   ```

3. **The dashboard will automatically be deployed** whenever you update the README or dashboard files.

## Troubleshooting

### Dashboard shows "Error loading data"

If your deployed dashboard shows an error loading data, it might be because:

1. **The data.json file was not embedded properly:**
   - Make sure you're using the build script to generate the HTML with embedded data
   - Check that the data.json file exists before running the build script

2. **CORS issues:**
   - The dashboard is designed to handle this with multiple fallback methods
   - Click the "Load Sample Data" button to see a demonstration with sample data

3. **Path issues:**
   - If your repository is not at root (e.g., yourusername.github.io/awesome-virome), the dashboard includes code to detect and adjust paths automatically
   - Check the browser console for any path-related errors

### Other common issues

1. **Files not found:**
   - Make sure all files (index.html, styles.css, visualization.js) are in the gh-pages branch
   
2. **GitHub Pages not enabled:**
   - Verify in repository settings that GitHub Pages is enabled for the gh-pages branch

3. **Cache issues:**
   - Try clearing your browser cache or viewing in an incognito/private window

## Updating the Dashboard

When you make changes to your README.md or want to update the dashboard:

1. Run the build process again: `node build.js`
2. Deploy the updated dist directory using your chosen method

For GitHub Actions (Option 3), updates will happen automatically when you push changes to the main branch.