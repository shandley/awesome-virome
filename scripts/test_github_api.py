#!/usr/bin/env python3
"""
Test script for GitHub API integration - Phase 1 Citation System V2

Tests the GitHub metrics collector with a small sample of tools.
This ensures our API integration works before processing all tools.
"""

import sys
import os
import json

# Add the current directory to Python path to import our module
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

try:
    from github_metrics_enhancer import GitHubMetricsCollector
except ImportError as e:
    print(f"Error importing GitHubMetricsCollector: {e}")
    sys.exit(1)

def test_github_api():
    """Test GitHub API with a few sample repositories."""

    print("Testing GitHub API integration...")

    # Initialize collector
    github_token = os.environ.get('GITHUB_TOKEN')
    collector = GitHubMetricsCollector(github_token)

    # Check rate limit
    rate_info = collector.check_rate_limit()
    if rate_info:
        print(f"✅ GitHub API accessible. Rate limit: {rate_info['remaining']}/{rate_info['limit']}")
    else:
        print("⚠️  Warning: Could not check GitHub API rate limit")

    # Test repository extraction
    test_urls = [
        "https://github.com/mtisza1/Cenote-Taker3",
        "https://github.com/bbuchfink/diamond",
        "https://github.com/apcamargo/genomad"
    ]

    print("\nTesting repository URL extraction:")
    for url in test_urls:
        repo_path = collector.extract_github_repo(url)
        print(f"  {url} → {repo_path}")

    # Test metrics collection for one repository
    print(f"\nTesting metrics collection for {test_urls[0]}:")
    repo_path = collector.extract_github_repo(test_urls[0])
    if repo_path:
        metrics = collector.get_repo_metrics(repo_path)
        if 'error' not in metrics:
            print("✅ Successfully collected metrics:")
            print(f"  - Stars: {metrics.get('stars', 'N/A')}")
            print(f"  - Forks: {metrics.get('forks', 'N/A')}")
            print(f"  - Language: {metrics.get('language', 'N/A')}")
            print(f"  - Last updated: {metrics.get('updated_at', 'N/A')}")
            print(f"  - Archived: {metrics.get('archived', 'N/A')}")
            if metrics.get('topics'):
                print(f"  - Topics: {', '.join(metrics['topics'])}")
        else:
            print(f"❌ Error collecting metrics: {metrics['error']}")
            return False

    print("\n✅ GitHub API integration test completed successfully!")
    return True

def test_with_sample_tools():
    """Test enhancement with sample tools from data.json."""

    print("\nTesting with sample tools from data.json...")

    # Load a few sample tools
    try:
        with open('data.json', 'r') as f:
            data = json.load(f)
    except Exception as e:
        print(f"Error loading data.json: {e}")
        return False

    # Find first few GitHub tools
    github_tools = []
    for node in data.get('nodes', []):
        if (node.get('type') == 'tool' and
            'github.com' in node.get('url', '') and
            len(github_tools) < 3):
            github_tools.append(node)

    if not github_tools:
        print("No GitHub tools found in data.json")
        return False

    print(f"Testing enhancement with {len(github_tools)} sample tools:")

    collector = GitHubMetricsCollector(os.environ.get('GITHUB_TOKEN'))

    for tool in github_tools:
        print(f"\nTesting: {tool['name']}")
        print(f"  URL: {tool['url']}")

        enhanced_tool = collector.enhance_tool_metrics(tool)

        if 'github_metrics' in enhanced_tool:
            print("  ✅ Enhancement successful!")
            print(f"    - Stars: {enhanced_tool['stars']}")
            print(f"    - Forks: {enhanced_tool['forks']}")
            print(f"    - Languages: {enhanced_tool.get('all_languages', [])}")
        else:
            print("  ❌ Enhancement failed")
            return False

    print("\n✅ Sample tool enhancement test completed successfully!")
    return True

def main():
    """Run all tests."""
    print("GitHub API Integration Test Suite")
    print("=" * 40)

    # Test 1: Basic API functionality
    if not test_github_api():
        print("❌ Basic API test failed")
        return 1

    # Test 2: Sample tool enhancement
    if not test_with_sample_tools():
        print("❌ Sample tool test failed")
        return 1

    print("\n🎉 All tests passed! Ready for Phase 1 implementation.")
    return 0

if __name__ == "__main__":
    sys.exit(main())