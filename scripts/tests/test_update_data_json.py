#!/usr/bin/env python3
"""Tests for the language handling in update_data_json.py."""
import sys
from pathlib import Path
from typing import Any, Optional

# update_data_json.py lives at the repository root
sys.path.append(str(Path(__file__).resolve().parents[2]))

from update_data_json import update_data_json  # noqa: E402


def make_tool(name: str, language: Optional[str] = None) -> dict[str, Any]:
    """A tool entry as extract_tools_from_readme() returns it."""
    return {
        "name": name,
        "url": f"https://github.com/example/{name}",
        "description": "Example tool.",
        "category": "Host Prediction",
        "subcategory": None,
        "language": language,
        "package_manager": None,
    }


def run_weekly(existing_nodes: list[dict[str, Any]], tools: list[dict[str, Any]]) -> dict[str, Any]:
    """Run the weekly path (no metadata pass) and return the result."""
    return update_data_json({"nodes": existing_nodes, "links": []}, tools)


def tool_node(data: dict[str, Any], tool_id: str) -> dict[str, Any]:
    return next(node for node in data["nodes"] if node["id"] == tool_id)


def test_untagged_tool_keeps_github_language() -> None:
    existing = [{"id": "tool-alpha", "type": "tool", "language": None,
                 "github_metrics": {"language": "Nextflow"}}]
    data = run_weekly(existing, [make_tool("alpha")])
    assert tool_node(data, "tool-alpha")["language"] == "Nextflow"
    assert "Nextflow" in data["languages"]


def test_untagged_tool_without_github_metrics_keeps_existing_language() -> None:
    existing = [{"id": "tool-beta", "type": "tool", "language": "Perl"}]
    data = run_weekly(existing, [make_tool("beta")])
    assert tool_node(data, "tool-beta")["language"] == "Perl"


def test_readme_tag_takes_precedence() -> None:
    existing = [{"id": "tool-gamma", "type": "tool", "language": "Shell",
                 "github_metrics": {"language": "Shell"}}]
    data = run_weekly(existing, [make_tool("gamma", language="Python")])
    assert tool_node(data, "tool-gamma")["language"] == "Python"


def test_new_untagged_tool_has_no_language() -> None:
    data = run_weekly([], [make_tool("delta")])
    assert tool_node(data, "tool-delta")["language"] is None
