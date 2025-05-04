# Updating Documentation

This guide explains how to update the Awesome-Virome documentation.

## Prerequisites

1. Python 3.7+ installed
2. Git installed
3. Clone of the repository

## Setting Up the Environment

```bash
# Navigate to the mkdocs directory
cd mkdocs

# Create a virtual environment
python -m venv venv

# Activate the virtual environment
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

## Making Documentation Changes

1. Edit Markdown files in the `docs/` directory
2. Add new files as needed (remember to update navigation in `mkdocs.yml`)
3. Test your changes locally:

```bash
mkdocs serve
```

4. View your changes at http://127.0.0.1:8000/

## Creating New Versions

We use [mike](https://github.com/jimporter/mike) for versioning:

```bash
# Deploy a new version (e.g., 1.1.0)
mike deploy 1.1.0

# Update the 'latest' alias
mike alias 1.1.0 latest

# Set the default version
mike set-default 1.1.0
```

## Submitting Changes

1. Commit your changes:

```bash
git add .
git commit -m "Update documentation: [brief description]"
```

2. Push to your fork:

```bash
git push origin your-branch-name
```

3. Create a pull request

## Documentation Structure

- `docs/index.md` - Home page
- `docs/intro/` - Introduction and getting started 
- `docs/tools/` - Tool documentation
- `docs/api/` - API documentation
- `docs/contributing/` - Contribution guidelines

## Style Guidelines

- Use Markdown headers appropriately (# for title, ## for sections)
- Include code examples when relevant
- Use relative links when linking to other documentation pages
- Use tables for presenting comparative information
- Include screenshots when helpful