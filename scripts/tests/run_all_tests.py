#!/usr/bin/env python3
"""
Run all citation workflow tests

This script runs all the tests for the citation workflows to ensure they
work as expected and do not generate synthetic data.
"""

import os
import sys
import subprocess
import argparse
from pathlib import Path

def run_test(test_script, *args, **kwargs):
    """
    Run a test script and report results
    
    Args:
        test_script: Path to the test script
        *args: Additional arguments to pass to the script
        **kwargs: Additional keyword arguments
        
    Returns:
        True if the test passed, False otherwise
    """
    script_name = os.path.basename(test_script)
    print(f"\n{'='*80}")
    print(f"Running test: {script_name}")
    print(f"{'='*80}")
    
    cmd = [sys.executable, test_script]
    cmd.extend(args)
    
    verbose = kwargs.get("verbose", False)
    
    process = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True
    )
    
    stdout, stderr = process.communicate()
    
    if process.returncode == 0:
        print(f"PASSED: {script_name} ✅")
        if verbose and stdout.strip():
            print("\nOutput:")
            print(stdout)
        return True
    else:
        print(f"FAILED: {script_name} ❌")
        print("\nStandard Output:")
        print(stdout)
        print("\nStandard Error:")
        print(stderr)
        return False

def main():
    """Main entry point for the script"""
    parser = argparse.ArgumentParser(description="Run all citation workflow tests")
    parser.add_argument("--verbose", "-v", action="store_true", help="Show verbose output")
    args = parser.parse_args()
    
    # Get the directory of this script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Define the test scripts to run
    test_scripts = [
        os.path.join(script_dir, "test_comprehensive_citation_data.py"),
        os.path.join(script_dir, "test_auto_fix_dois.py"),
        os.path.join(script_dir, "validate_workflows.py"),
        os.path.join(script_dir, "run_workflows_dryrun.py"),
        # Run verify_impact_data.py only if impact_data.json exists
        os.path.join(script_dir, "verify_impact_data.py")
    ]
    
    # Run each test script
    results = {}
    for test_script in test_scripts:
        if os.path.exists(test_script):
            results[test_script] = run_test(test_script, verbose=args.verbose)
        else:
            print(f"WARNING: Test script not found: {test_script}")
            results[test_script] = False
    
    # Print summary
    print("\n\nTest Results Summary:")
    print("=====================")
    
    all_passed = True
    for test_script, passed in results.items():
        script_name = os.path.basename(test_script)
        status = "PASSED" if passed else "FAILED"
        print(f"{script_name}: {status}")
        if not passed:
            all_passed = False
    
    # Exit with appropriate status code
    if all_passed:
        print("\nAll tests passed! ✅")
        print("The citation workflows are functioning as expected and do not generate synthetic data.")
        sys.exit(0)
    else:
        print("\nSome tests failed. ❌")
        print("Please fix the issues before pushing to GitHub.")
        sys.exit(1)

if __name__ == "__main__":
    main()