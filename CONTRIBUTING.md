# Contributing to Awesome-Virome

Thank you for considering contributing to this curated list of software, tools, and databases for virome analysis!

## How to Contribute

There are two ways to contribute a new tool:

### Option 1: Submit an Issue
1. Go to the [Issues page](https://github.com/shandley/awesome-virome/issues) and click "New Issue"
2. Select the "New Tool Submission" template
3. Fill out all required fields with information about your tool
4. Submit the issue

### Option 2: Submit a Pull Request
1. Fork the repository
2. Create a new branch for your addition (`git checkout -b add-new-tool`)
3. Add your entry to the appropriate section in README.md
4. Submit a pull request

## Automated Validation

All tool submissions undergo automated validation to ensure data quality:

1. **Initial Validation**: When you submit an issue or PR, our validation system will:
   - Check that all required fields are provided
   - Verify the tool URL is accessible
   - Validate any DOIs or citations
   - Analyze description quality
   - Check for duplicate entries
   - Assign a quality score (0-100) to your submission

2. **Feedback Process**: 
   - You'll receive immediate feedback as a comment on your issue/PR
   - If validation fails, you'll see specific errors that need to be fixed
   - Warnings and suggestions will help you improve your submission quality
   - Submissions with a quality score below 70 will need improvements

3. **Final Review**:
   - After passing automated validation, a maintainer will manually review your submission
   - Once approved, your tool will be added to the collection

## Guidelines for Entries

Each entry should follow this format:
```
- [Tool Name](link-to-tool) - Brief description of the tool. [package-manager] [language]
```

Where:
- **Tool Name** is the name of the software, tool, or database
- **link-to-tool** is a URL to the tool's website, GitHub repository, or paper
- **Brief description** explains what the tool does in 1-2 sentences
- **package-manager** indicates how the tool can be installed (conda, pip, etc.)
- **language** indicates the programming language used (optional)

## Example Entry

```
- [ViralMSA](https://github.com/niemasd/ViralMSA) - Python script for viral multiple sequence alignment using read mappers. [source] [Python]
```

## Key Quality Requirements

Your submission will be evaluated against these key criteria:

1. **Basic Requirements**:
   - Tool must be relevant to virome analysis
   - URL must be accessible and point to the correct resource
   - Description must be clear and include virome-related terms
   - Tools must be placed in the correct category

2. **Quality Enhancement**:
   - Include DOI or citation information when available
   - Provide installation instructions (conda, pip, docker, etc.)
   - Mention programming language and dependencies
   - Note license information if applicable
   - Indicate maintenance status (especially if no longer maintained)

3. **Repository Information**:
   - For GitHub, GitLab, or Bitbucket repositories:
     - Ensure the repository has a README with usage information
     - Check that the tool is not archived unless historically significant
     - Verify recent activity or note if not actively maintained

## Enhanced Metadata Collection

This repository automatically collects enhanced metadata from GitHub, GitLab, and Bitbucket repositories:

1. **Weekly Updates**: Basic repository information is updated weekly (stars, update time)
2. **Monthly Metadata**: Detailed metadata is collected monthly including:
   - Programming languages
   - License information
   - Repository topics/tags
   - Release information
   - Creation and update dates
   - Dependencies (when available)

If you add a tool that is hosted on GitHub, GitLab, or Bitbucket, this metadata will be automatically collected and incorporated into the data visualization. For other tools, only the information you provide in the README will be used.

## Citation Validation

All DOIs and citations undergo validation to ensure they remain accurate:

1. DOIs are checked for proper formatting and resolveability
2. Citations are validated for compliance with standard formats
3. Consistency is ensured between DOIs and citation text

Providing accurate citation information helps researchers properly credit the tools they use.

Thank you for your contributions!