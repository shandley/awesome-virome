#!/usr/bin/env python3
"""
Phase 1 Demo - GitHub Metrics Enhancement

A minimal working example that enhances a small subset of tools
to demonstrate the Phase 1 citation system functionality.

This is safe to run and will only modify a backup copy of data.json.
"""

import json
import os
import sys
from datetime import datetime

# Add the current directory to Python path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

try:
    from github_metrics_enhancer import GitHubMetricsCollector
except ImportError as e:
    print(f"Error importing GitHubMetricsCollector: {e}")
    sys.exit(1)

def create_backup():
    """Create a backup of data.json before modification."""
    backup_filename = f"data_backup_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    try:
        with open('data.json', 'r') as source:
            with open(backup_filename, 'w') as backup:
                backup.write(source.read())
        print(f"✅ Backup created: {backup_filename}")
        return backup_filename
    except Exception as e:
        print(f"❌ Error creating backup: {e}")
        return None

def run_phase1_demo():
    """Run Phase 1 demo with a small subset of tools."""

    print("Phase 1 Demo: GitHub Metrics Enhancement")
    print("=" * 45)

    # Create backup first
    backup_file = create_backup()
    if not backup_file:
        return 1

    # Load data
    try:
        with open('data.json', 'r') as f:
            data = json.load(f)
    except Exception as e:
        print(f"❌ Error loading data.json: {e}")
        return 1

    # Find first 5 GitHub tools for demo
    github_tools = []
    for node in data.get('nodes', []):
        if (node.get('type') == 'tool' and
            'github.com' in node.get('url', '') and
            len(github_tools) < 5):
            github_tools.append(node)

    if not github_tools:
        print("❌ No GitHub tools found in data.json")
        return 1

    print(f"Demo will enhance {len(github_tools)} GitHub tools:")
    for tool in github_tools:
        print(f"  - {tool['name']}: {tool['url']}")

    # Initialize collector
    github_token = os.environ.get('GITHUB_TOKEN')
    if github_token:
        print(f"✅ Using authenticated GitHub API (higher rate limits)")
    else:
        print("⚠️  Using unauthenticated GitHub API (60 requests/hour)")

    collector = GitHubMetricsCollector(github_token)

    # Check rate limit
    rate_info = collector.check_rate_limit()
    if rate_info:
        print(f"📊 Rate limit: {rate_info['remaining']}/{rate_info['limit']} remaining")

    print("\nProcessing tools...")
    print("-" * 30)

    enhanced_count = 0
    for i, tool in enumerate(github_tools):
        print(f"\n{i+1}. Processing: {tool['name']}")

        # Store original values for comparison
        original_stars = tool.get('stars', 'unknown')
        original_forks = tool.get('forks', 'unknown')

        # Enhance metrics
        enhanced_tool = collector.enhance_tool_metrics(tool)

        if 'github_metrics' in enhanced_tool:
            enhanced_count += 1

            # Show enhancement results
            new_stars = enhanced_tool.get('stars', 'unknown')
            new_forks = enhanced_tool.get('forks', 'unknown')

            print(f"   ✅ Enhanced successfully!")
            print(f"   📊 Stars: {original_stars} → {new_stars}")
            print(f"   🍴 Forks: {original_forks} → {new_forks}")

            if enhanced_tool.get('github_metrics', {}).get('language'):
                print(f"   💻 Language: {enhanced_tool['github_metrics']['language']}")

            if enhanced_tool.get('github_metrics', {}).get('topics'):
                topics = enhanced_tool['github_metrics']['topics'][:3]  # Show first 3
                if topics:
                    print(f"   🏷️  Topics: {', '.join(topics)}")

            if enhanced_tool.get('github_metrics', {}).get('recent_activity') is not None:
                activity = "🟢 Active" if enhanced_tool['github_metrics']['recent_activity'] else "🟡 Quiet"
                print(f"   📈 Recent Activity: {activity}")

            # Update the tool in the data structure
            for j, node in enumerate(data['nodes']):
                if node.get('id') == tool['id']:
                    data['nodes'][j] = enhanced_tool
                    break

        else:
            print(f"   ❌ Enhancement failed")

    print(f"\n📊 Summary:")
    print(f"   - Total tools processed: {len(github_tools)}")
    print(f"   - Successfully enhanced: {enhanced_count}")
    print(f"   - Success rate: {enhanced_count/len(github_tools)*100:.1f}%")

    # Save enhanced data
    output_file = f"data_phase1_demo_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    try:
        with open(output_file, 'w') as f:
            json.dump(data, f, indent=2, ensure_ascii=False)
        print(f"✅ Enhanced data saved to: {output_file}")
        print(f"🔒 Original data.json unchanged (backup: {backup_file})")
    except Exception as e:
        print(f"❌ Error saving enhanced data: {e}")
        return 1

    print(f"\n🎉 Phase 1 Demo completed successfully!")
    print(f"💡 To apply changes permanently, replace data.json with {output_file}")

    return 0

def main():
    """Main function."""
    if len(sys.argv) > 1 and sys.argv[1] == '--help':
        print(__doc__)
        return 0

    return run_phase1_demo()

if __name__ == "__main__":
    sys.exit(main())