# Submission Validation Process

To maintain high quality data in the Awesome-Virome collection, all tool submissions undergo a comprehensive validation process. This page explains how submissions are validated and provides guidelines to help your submission pass validation.

## Validation Overview

The validation process includes:

1. **Automated Validation**: Initial checks performed by our validation system
2. **Manual Review**: Human review of submissions that pass automated validation
3. **Feedback Process**: Communication about issues or improvements needed
4. **Final Approval**: Acceptance and integration into the collection

## Automated Validation

### What Gets Validated

Our automated system checks:

- **Required Fields**: Ensuring all required information is provided
- **URL Accessibility**: Verifying that the tool URL is accessible
- **DOI/Citation Validity**: Validating any provided DOIs or citations
- **Description Quality**: Analyzing the clarity and informativeness of descriptions
- **Duplicate Detection**: Checking for existing entries with the same name or URL
- **Category Correctness**: Verifying appropriate categorization

### Validation Scoring

Each submission receives a quality score (0-100) based on:

| Component | Weight | Description |
|-----------|--------|-------------|
| Completeness | 30% | Presence of all required and recommended fields |
| URL Validity | 20% | Accessibility and correctness of the provided URL |
| Description Quality | 25% | Clarity, relevance, and informativeness |
| Metadata Accuracy | 15% | Correctness of provided metadata (version, language, etc.) |
| Citation Validity | 10% | Validity of provided DOIs or citations |

Submissions must achieve a score of at least 70 to proceed to manual review.

## Validation Process

### Step 1: Initial Submission

When you submit a new tool (via issue or pull request), the validation process begins automatically.

### Step 2: Automated Checks

Our validation system performs the following checks:

1. **URL Check**: The tool URL is tested to ensure it's accessible
2. **Required Fields Check**: Verifies all required fields are provided
3. **Duplicate Check**: Searches for existing entries with similar names or URLs
4. **DOI Validation**: If a DOI is provided, it's checked for validity
5. **Description Analysis**: Evaluates the description for clarity and relevance
6. **Repository Analysis**: For GitHub/GitLab/Bitbucket repos, additional metadata is collected

### Step 3: Validation Report

After validation completes, you'll receive a comment with:

- Overall validation score
- Pass/fail status for each validation check
- Specific errors or warnings that need addressing
- Suggestions for improving your submission

Example validation report:

```
## Validation Report for ToolName

Overall Score: 85/100 ✅

### Validation Results:
- Required Fields: ✅ All required fields present
- URL Validity: ✅ URL is accessible
- Duplicate Check: ✅ No duplicates found
- Description Quality: ⚠️ Description could be more detailed
- Repository Analysis: ✅ GitHub metadata successfully collected

### Suggestions:
- Add more specific details about what the tool does
- Consider including installation instructions
- Add package manager information if available

This submission has passed validation and will be reviewed by a maintainer.
```

### Step 4: Addressing Feedback

If your submission receives a failing score or has errors:

1. Review the validation report carefully
2. Make the necessary changes to your submission
3. The validation will run again automatically after your changes

### Step 5: Manual Review

Submissions that pass automated validation undergo manual review by a maintainer, who checks:

- Relevance to virome analysis
- Appropriate categorization
- Description accuracy and clarity
- Overall quality and usefulness

### Step 6: Final Decision

After manual review, your submission will either be:

- **Approved**: Integrated into the collection
- **Requested Changes**: Sent back with specific improvement requests
- **Rejected**: Declined with a clear explanation of why

## Common Validation Issues

### URL Validation Failures

**Issue**: The provided URL is not accessible or returns an error.

**Solution**:
- Verify the URL is correct and not missing components
- Ensure the repository or website is public
- Check if the tool has moved to a new location
- For archived tools, provide an archive URL (e.g., Internet Archive)

### Description Quality Issues

**Issue**: The description is too vague, too short, or lacks relevant information.

**Solution**:
- Ensure the description is at least 100 characters
- Include specific functionality the tool provides
- Mention what type of data it works with
- Include unique features or advantages

### Duplicate Detection

**Issue**: A tool with the same or very similar name already exists.

**Solution**:
- Check if you're submitting an updated version of an existing tool
- Clarify the name if it's different from an existing tool
- Add version information if it's a new version
- Consider updating the existing entry instead

### Category Mismatches

**Issue**: The tool is placed in an inappropriate category.

**Solution**:
- Review the category definitions in the README
- Consider the primary functionality of the tool
- Place tools with multiple functions in the most relevant category
- Suggest a new category if none of the existing ones fit

## Example of a Good Submission

Here's an example of a submission that would pass validation:

```
## Tool Information

- **Name**: ViralMSA
- **URL**: https://github.com/niemasd/ViralMSA
- **Description**: A Python tool for viral multiple sequence alignment using various read mappers. It's optimized for viral genome analysis and supports multiple reference genomes.
- **Category**: Sequence Analysis
- **Subcategory**: Multiple Sequence Alignment
- **Installation**: pip install viralmsa
- **Language**: Python
- **Version**: v1.1.2 (2023-02-15)
- **License**: GNU GPL v3.0
- **Paper DOI**: 10.1093/bioinformatics/btaa743
```

This submission includes all required fields, provides comprehensive information, and correctly categorizes the tool.

## Validation Tools

For contributors who want to validate submissions before submitting:

- Our validation script is available in the repository at `scripts/tool_validator.py`
- You can run it locally to check your submission before creating an issue or PR
- Example usage: `python scripts/tool_validator.py --name "ToolName" --url "https://example.com/tool" --description "Tool description" --category "Category"`

## Questions About Validation

If you have questions about the validation process:

- Check the [validation FAQ](https://github.com/shandley/awesome-virome/wiki/Validation-FAQ)
- Ask in [GitHub Discussions](https://github.com/shandley/awesome-virome/discussions/categories/tool-submission)
- Contact the maintainers through the repository

We appreciate your contributions to Awesome-Virome and are committed to maintaining a high-quality collection of virome analysis tools.