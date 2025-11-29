#!/usr/bin/env python3
"""
Quick script to verify GitHub token is working properly.
"""

import os
import requests

def verify_github_token():
    """Verify GitHub token configuration."""

    token = os.environ.get('GITHUB_TOKEN')
    if not token:
        print("❌ GITHUB_TOKEN environment variable not found")
        print("\nTo set your token, run:")
        print('export GITHUB_TOKEN="your_token_here"')
        return False

    print("✅ GITHUB_TOKEN environment variable is set")
    print(f"   Token length: {len(token)} characters")
    print(f"   Token prefix: {token[:4]}...")

    # Test the token
    headers = {
        'Authorization': f'token {token}',
        'Accept': 'application/vnd.github.v3+json'
    }

    try:
        response = requests.get('https://api.github.com/rate_limit', headers=headers)

        if response.status_code == 200:
            data = response.json()
            print("✅ Token is valid and working!")
            print(f"   Rate limit: {data['rate']['remaining']}/{data['rate']['limit']} remaining")
            print(f"   Reset time: {data['rate']['reset']}")

            # Check user info
            user_response = requests.get('https://api.github.com/user', headers=headers)
            if user_response.status_code == 200:
                user_data = user_response.json()
                print(f"   Authenticated as: {user_data.get('login', 'unknown')}")

            # Estimate processing time
            tools_count = 164
            calls_per_tool = 2  # repo + languages API calls
            total_calls_needed = tools_count * calls_per_tool
            remaining_calls = data['rate']['remaining']

            print(f"\n📊 Processing Estimates:")
            print(f"   GitHub tools to process: {tools_count}")
            print(f"   API calls needed: ~{total_calls_needed}")
            print(f"   Current remaining calls: {remaining_calls}")

            if remaining_calls >= total_calls_needed:
                print("   ✅ Sufficient rate limit for full processing!")
                print("   🕐 Estimated time: 10-15 minutes")
            else:
                print("   ⚠️  May need to wait for rate limit reset")

            return True

        else:
            print(f"❌ Token test failed: HTTP {response.status_code}")
            if response.status_code == 401:
                print("   This usually means the token is invalid or expired")
            elif response.status_code == 403:
                print("   This usually means the token lacks required permissions")
            print(f"   Response: {response.text[:200]}")
            return False

    except Exception as e:
        print(f"❌ Error testing token: {e}")
        return False

if __name__ == "__main__":
    if verify_github_token():
        print("\n🎉 Ready to run full dataset enhancement!")
    else:
        print("\n🔧 Please fix the token setup before proceeding.")