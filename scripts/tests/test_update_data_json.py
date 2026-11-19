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
    existing = [{"id": "tool-alpha", "type": "tool", "url": "https://github.com/example/alpha", "language": None,
                 "github_metrics": {"language": "Nextflow"}}]
    data = run_weekly(existing, [make_tool("alpha")])
    assert tool_node(data, "tool-alpha")["language"] == "Nextflow"
    assert "Nextflow" in data["languages"]


def test_untagged_tool_without_github_metrics_keeps_existing_language() -> None:
    existing = [{"id": "tool-beta", "type": "tool", "url": "https://github.com/example/beta", "language": "Perl"}]
    data = run_weekly(existing, [make_tool("beta")])
    assert tool_node(data, "tool-beta")["language"] == "Perl"


def test_readme_tag_takes_precedence() -> None:
    existing = [{"id": "tool-gamma", "type": "tool", "url": "https://github.com/example/gamma", "language": "Shell",
                 "github_metrics": {"language": "Shell"}}]
    data = run_weekly(existing, [make_tool("gamma", language="Python")])
    assert tool_node(data, "tool-gamma")["language"] == "Python"


def test_new_untagged_tool_has_no_language() -> None:
    data = run_weekly([], [make_tool("delta")])
    assert tool_node(data, "tool-delta")["language"] is None


def test_moved_tool_drops_metrics_from_old_repository() -> None:
    existing = [{"id": "tool-epsilon", "type": "tool", "url": "https://bitbucket.org/example/epsilon",
                 "stars": 2, "color": "#198754", "doi": "10.1000/example", "provider": "bitbucket",
                 "bitbucket_metrics": {"watchers": 2}, "maintenance_status": "active"}]
    node = tool_node(run_weekly(existing, [make_tool("epsilon")]), "tool-epsilon")
    assert "bitbucket_metrics" not in node and "provider" not in node
    assert node["stars"] == 0
    assert node["maintenance_status"] == "unknown"
    assert node["doi"] == "10.1000/example"


def test_trivial_url_edit_keeps_metrics() -> None:
    existing = [{"id": "tool-zeta", "type": "tool", "url": "http://GitHub.com/example/zeta/",
                 "stars": 40, "github_metrics": {"language": "Rust"}}]
    node = tool_node(run_weekly(existing, [make_tool("zeta")]), "tool-zeta")
    assert node["stars"] == 40
    assert node["language"] == "Rust"
