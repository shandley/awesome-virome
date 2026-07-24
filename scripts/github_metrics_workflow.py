#!/usr/bin/env python3
"""
GitHub Metrics Workflow - Production Version for GitHub Actions

A production-ready version of the GitHub metrics enhancer specifically designed
for GitHub Actions workflows. Includes enhanced logging, error handling, and
integration with the existing workflow system.

Features:
- Production-level logging and error reporting
- Integration with existing workflow caching
- GitHub Actions-friendly output formatting
- Automatic retry logic for transient failures
- Comprehensive metrics reporting
"""

import json
import sys
import os
from datetime import datetime
import time

# Import our Phase 1 system
try:
    from github_metrics_enhancer import GitHubMetricsCollector, enhance_bitbucket_tool
except ImportError:
    print("❌ Error: Could not import GitHubMetricsCollector")
    print("   Make sure scripts/github_metrics_enhancer.py is available")
    sys.exit(1)

def setup_workflow_logging():
    """Configure logging for GitHub Actions environment."""
    import logging

    # Configure logging with GitHub Actions-friendly format
    logging.basicConfig(
        level=logging.INFO,
        format='%(levelname)s: %(message)s',
        handlers=[logging.StreamHandler(sys.stdout)]
    )

    return logging.getLogger(__name__)

def github_actions_summary(stats):
    """Create GitHub Actions job summary."""
    summary_file = os.environ.get('GITHUB_STEP_SUMMARY')
    if not summary_file:
        return

    try:
        with open(summary_file, 'a') as f:
            f.write(f"""
## 🎯 Citation System V2 Phase 1 Results

### GitHub Metrics Enhancement Summary
- **Total GitHub tools**: {stats['total_tools']}
- **Successfully enhanced**: {stats['enhanced_count']}
- **Success rate**: {stats['success_rate']:.1f}%
- **Processing time**: {stats['processing_time']:.1f} seconds
- **API calls made**: {stats['api_calls']}

### Key Improvements
{stats['improvements']}

### Rate Limit Status
- **Remaining calls**: {stats['rate_limit_remaining']}
- **Rate limit**: {stats['rate_limit_total']}
""")
    except Exception as e:
        print(f"Warning: Could not write job summary: {e}")

def run_workflow_enhancement():
    """Run GitHub metrics enhancement in workflow environment."""

    logger = setup_workflow_logging()
    start_time = time.time()

    print("🚀 Starting Citation System V2 Phase 1: GitHub Metrics Enhancement")
    print("=" * 70)

    # Check for GitHub token
    github_token = os.environ.get('GITHUB_TOKEN')
    if not github_token:
        print("❌ Error: GITHUB_TOKEN environment variable not found")
        print("   This is required for GitHub Actions workflow")
        return 1

    # Initialize collector
    collector = GitHubMetricsCollector(github_token)

    # Check rate limit
    print("📊 Checking GitHub API rate limits...")
    rate_info = collector.check_rate_limit()
    if rate_info:
        print(f"   Rate limit: {rate_info['remaining']}/{rate_info['limit']} remaining")
        if rate_info['remaining'] < 100:
            print("⚠️  Warning: Low rate limit remaining. Enhancement may be incomplete.")

    # Load existing data
    try:
        with open('data.json', 'r') as f:
            data = json.load(f)
        print("✅ Successfully loaded data.json")
    except Exception as e:
        print(f"❌ Error loading data.json: {e}")
        return 1

    # Find GitHub tools
    github_tools = []
    for node in data.get('nodes', []):
        if node.get('type') == 'tool' and 'github.com' in node.get('url', ''):
            github_tools.append(node)

    print(f"🔍 Found {len(github_tools)} GitHub tools to enhance")

    if not github_tools:
        print("⚠️  No GitHub tools found in data.json")
        return 0

    # Process tools
    enhanced_count = 0
    api_calls = 0
    improvements = []

    print("\n🔄 Processing GitHub tools...")
    print("-" * 40)

    for i, tool in enumerate(github_tools):
        tool_name = tool.get('name', 'Unknown')
        print(f"   {i+1:3d}/{len(github_tools)} {tool_name:<30}", end=" ")

        # Store original values for comparison
        original_stars = tool.get('stars', 0)

        try:
            enhanced_tool = collector.enhance_tool_metrics(tool)

            if 'github_metrics' in enhanced_tool:
                enhanced_count += 1
                api_calls += 2  # Estimate: repo + languages API calls

                # Track improvements
                new_stars = enhanced_tool.get('stars', 0)
                if new_stars != original_stars:
                    improvements.append(f"- **{tool_name}**: {original_stars} → {new_stars} stars")

                # Update the tool in the data structure
                for j, node in enumerate(data['nodes']):
                    if node.get('id') == tool['id']:
                        data['nodes'][j] = enhanced_tool
                        break

                print("✅")
            else:
                print("❌")

        except Exception as e:
            print(f"❌ ({str(e)[:20]}...)")
            logger.warning(f"Failed to enhance {tool_name}: {e}")

        # Progress update every 20 tools
        if (i + 1) % 20 == 0:
            elapsed = time.time() - start_time
            print(f"   Progress: {i+1}/{len(github_tools)} tools ({elapsed:.1f}s elapsed)")

    # Save enhanced data
    print(f"\n💾 Saving enhanced data...")
    try:
        with open('data.json', 'w') as f:
            json.dump(data, f, indent=2, ensure_ascii=False)
        print("✅ Successfully saved enhanced data.json")
    except Exception as e:
        print(f"❌ Error saving data.json: {e}")
        return 1

    # Final statistics
    processing_time = time.time() - start_time
    success_rate = (enhanced_count / len(github_tools)) * 100 if github_tools else 0

    print(f"\n📊 Enhancement Summary:")
    print(f"   Total tools processed: {len(github_tools)}")
    print(f"   Successfully enhanced: {enhanced_count}")
    print(f"   Success rate: {success_rate:.1f}%")
    print(f"   Processing time: {processing_time:.1f} seconds")
    print(f"   API calls made: ~{api_calls}")

    # Check final rate limit
    final_rate_info = collector.check_rate_limit()
    if final_rate_info:
        print(f"   Remaining rate limit: {final_rate_info['remaining']}/{final_rate_info['limit']}")

    # Show improvements
    if improvements:
        print(f"\n⭐ Star count updates detected:")
        for improvement in improvements[:5]:  # Show first 5
            print(f"   {improvement}")
        if len(improvements) > 5:
            print(f"   ... and {len(improvements) - 5} more updates")

    # Create GitHub Actions summary
    stats = {
        'total_tools': len(github_tools),
        'enhanced_count': enhanced_count,
        'success_rate': success_rate,
        'processing_time': processing_time,
        'api_calls': api_calls,
        'improvements': '\n'.join(improvements[:10]) if improvements else 'No star count changes detected',
        'rate_limit_remaining': final_rate_info['remaining'] if final_rate_info else 'Unknown',
        'rate_limit_total': final_rate_info['limit'] if final_rate_info else 'Unknown'
    }
    github_actions_summary(stats)

    print(f"\n🎉 Citation System V2 Phase 1 completed successfully!")

    if enhanced_count == len(github_tools):
        print("✨ Perfect run - all tools enhanced successfully!")
    elif enhanced_count > len(github_tools) * 0.9:
        print("🎯 Excellent run - over 90% success rate!")
    elif enhanced_count > len(github_tools) * 0.8:
        print("👍 Good run - over 80% success rate!")
    else:
        print("⚠️  Some tools failed to enhance - check logs for details")

    return 0

def run_bitbucket_enhancement():
    """Enhance Bitbucket-hosted tools in data.json.

    Runs independently of GitHub auth (Bitbucket's public API needs no token),
    so it happens on every workflow run even when GITHUB_TOKEN is absent.
    """
    print("\n🪣 Starting Bitbucket Metrics Enhancement")
    print("=" * 70)

    try:
        with open('data.json', 'r') as f:
            data = json.load(f)
    except Exception as e:
        print(f"❌ Error loading data.json: {e}")
        return 1

    bitbucket_tools = [
        node for node in data.get('nodes', [])
        if node.get('type') == 'tool' and 'bitbucket.org' in node.get('url', '')
    ]
    print(f"🔍 Found {len(bitbucket_tools)} Bitbucket tools to enhance")

    if not bitbucket_tools:
        return 0

    enhanced_count = 0
    for i, tool in enumerate(bitbucket_tools):
        tool_name = tool.get('name', 'Unknown')
        print(f"   {i+1:3d}/{len(bitbucket_tools)} {tool_name:<30}", end=" ")
        try:
            enhanced_tool = enhance_bitbucket_tool(tool)
            if 'bitbucket_metrics' in enhanced_tool:
                enhanced_count += 1
                for j, node in enumerate(data['nodes']):
                    if node.get('id') == tool['id']:
                        data['nodes'][j] = enhanced_tool
                        break
                print("✅")
            else:
                print("❌ (unverified, left as-is)")
        except Exception as e:
            print(f"❌ ({str(e)[:20]}...)")

    print(f"\n💾 Saving Bitbucket-enhanced data...")
    try:
        with open('data.json', 'w') as f:
            json.dump(data, f, indent=2, ensure_ascii=False)
        print(f"✅ Enhanced {enhanced_count}/{len(bitbucket_tools)} Bitbucket tools")
    except Exception as e:
        print(f"❌ Error saving data.json: {e}")
        return 1

    return 0


def main():
    """Main entry point for workflow."""
    try:
        # Bitbucket first: it needs no auth and must run even without a token,
        # unlike the GitHub pass which returns early when GITHUB_TOKEN is absent.
        run_bitbucket_enhancement()
        return run_workflow_enhancement()
    except KeyboardInterrupt:
        print("\n❌ Enhancement cancelled by user")
        return 1
    except Exception as e:
        print(f"\n❌ Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())