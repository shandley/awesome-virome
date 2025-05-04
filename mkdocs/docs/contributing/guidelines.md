# Contributing to Awesome-Virome

Thank you for considering contributing to Awesome-Virome! This document provides guidelines and instructions for contributing to this curated list of virome analysis tools.

## Ways to Contribute

There are several ways to contribute to Awesome-Virome:

1. **Add New Tools**: Submit entries for tools not yet included in the list
2. **Update Information**: Help keep tool descriptions and links up to date
3. **Report Issues**: Let us know about broken links or outdated information
4. **Improve Documentation**: Suggest enhancements to the documentation
5. **Fix Bugs**: Help resolve issues with the website or dashboard

## Contribution Process

### Option 1: Submit an Issue

The easiest way to contribute is by submitting an issue:

1. Go to the [Issues page](https://github.com/shandley/awesome-virome/issues) and click "New Issue"
2. Select the appropriate template:
   - **Tool Addition**: For suggesting new tools
   - **Tool Update**: For updating existing tool information
   - **Bug Report**: For reporting issues with the repository
   - **Feature Request**: For suggesting new features
3. Fill out all required fields
4. Submit the issue

Our team will review your submission and incorporate valid contributions.

### Option 2: Submit a Pull Request

For more direct contributions:

1. Fork the repository
2. Create a new branch for your addition (`git checkout -b add-new-tool`)
3. Make your changes
4. Submit a pull request

## Guidelines for Tool Entries

Each tool entry should follow this format:

```markdown
- [Tool Name](link-to-tool) - Brief description of the tool. [package-manager] [language]
```

Where:
- **Tool Name**: The name of the software, tool, or database
- **link-to-tool**: A URL to the tool's website, GitHub repository, or paper
- **Brief description**: Explains what the tool does in 1-2 sentences
- **package-manager**: Indicates how the tool can be installed (conda, pip, etc.)
- **language**: Indicates the programming language used (optional)

Example:

```markdown
- [ViralMSA](https://github.com/niemasd/ViralMSA) - Python script for viral multiple sequence alignment using read mappers. [source] [Python]
```

## Detailed Requirements

### Required Information

Every tool submission must include:

1. **Tool Name**: The official name of the tool
2. **URL**: A link to the tool's repository, website, or publication
3. **Description**: A clear, concise description of what the tool does
4. **Category**: The appropriate category for the tool

### Optional But Recommended

For higher quality entries, include:

1. **Version Information**: The latest release version, if available
2. **Installation Method**: How to install the tool (conda, pip, docker, etc.)
3. **Programming Language**: The primary language the tool is written in
4. **License Information**: The license under which the tool is distributed
5. **Citation Information**: DOI or publication reference
6. **GitHub Stars**: For GitHub repositories (automatically collected)
7. **Maintenance Status**: Whether the tool is actively maintained

## Quality Requirements

Your submission will be evaluated against these key criteria:

1. **Relevance**: The tool must be relevant to virome analysis
2. **Accessibility**: The URL must be accessible and correct
3. **Clarity**: The description must be clear and informative
4. **Categorization**: The tool must be placed in the correct category

## Automated Validation

All submissions undergo automated validation:

1. URLs are checked for accessibility
2. Descriptions are analyzed for clarity and relevance
3. Submissions are checked for duplicates
4. DOIs and citations are validated when provided

A quality score (0-100) is assigned based on these checks. Submissions scoring below 70 require improvements before acceptance.

## Enhanced Metadata Collection

For tools hosted on GitHub, GitLab, or Bitbucket, we automatically collect:

- Repository statistics (stars, forks, open issues)
- License information
- Programming languages
- Repository topics/tags
- Release information
- Creation and update dates

This metadata is updated weekly/monthly and used in our visualizations.

## Code of Conduct

By participating in this project, you agree to abide by our Code of Conduct:

- Be respectful and inclusive
- Focus on the technical merits of submissions
- Provide constructive feedback
- Maintain a positive and helpful environment

## Attribution

We appreciate all contributions to Awesome-Virome. Contributors are recognized in the following ways:

- Listed in the Contributors section of the repository
- Acknowledged in release notes for significant contributions
- Credited in the changelog

## Questions?

If you have questions about contributing, please:

1. Check the [FAQ](https://github.com/shandley/awesome-virome/wiki/FAQ)
2. Ask in [GitHub Discussions](https://github.com/shandley/awesome-virome/discussions)
3. Open an issue labeled "question"

Thank you for your contributions to the virome analysis community!