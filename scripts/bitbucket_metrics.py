#!/usr/bin/env python3
"""
Bitbucket Metrics Collector

Collects repository metrics for Bitbucket-hosted tools using the public
Bitbucket Cloud API 2.0 (no authentication required for public repos).

Bitbucket has no concept of "stars"; the watcher count is used as the
stars-equivalent so the dashboard treats Bitbucket and GitHub tools alike.
"""

import re
import time
from typing import Dict, Optional

import requests

API_BASE = "https://api.bitbucket.org/2.0/repositories"
TIMEOUT = 15
MAX_RETRIES = 3  # for transient HTTP 429 rate limiting


def _get(session: requests.Session, url: str, params: Optional[Dict] = None):
    """GET with simple backoff retries on HTTP 429 (Bitbucket rate limiting)."""
    resp = None
    for attempt in range(MAX_RETRIES):
        resp = session.get(url, params=params, timeout=TIMEOUT)
        if resp.status_code != 429:
            return resp
        if attempt < MAX_RETRIES - 1:
            time.sleep(2 ** attempt)  # 1s, 2s
    return resp


def extract_bitbucket_repo(url: str) -> Optional[str]:
    """Extract workspace/repo from a Bitbucket URL.

    Handles URLs with trailing paths such as ``/src/main/`` or ``/src/master/``
    by keeping only the first two path segments after the host.
    """
    match = re.search(r"bitbucket\.org/([^/]+/[^/]+)", url)
    if not match:
        return None
    repo_path = match.group(1)
    # Keep only workspace/repo, dropping any '/src/...' style suffix.
    parts = repo_path.split("/")
    return "/".join(parts[:2])


def get_repo_metrics(repo_path: str) -> Dict:
    """Fetch metrics for a Bitbucket repo.

    Returns a dict with keys: stars (=watchers), forks, language,
    pushed_at (=updated_on), archived (always False; Bitbucket has no archive
    flag). On any failure returns ``{'error': <message>}`` and never fabricates
    numbers.
    """
    session = requests.Session()
    try:
        resp = _get(session, f"{API_BASE}/{repo_path}")
        if resp.status_code != 200:
            return {"error": f"HTTP {resp.status_code}"}
        repo = resp.json()

        # Guard against a redirect/mismatch: the returned slug should match.
        expected_slug = repo_path.split("/")[-1].lower()
        if repo.get("slug", "").lower() != expected_slug:
            return {"error": f"slug mismatch: got {repo.get('slug')!r}"}

        # Coerce '' -> None; capitalize to match GitHub's convention ('Python',
        # not Bitbucket's lowercase 'python') so language grouping stays uniform.
        raw_language = repo.get("language")
        language = raw_language.capitalize() if raw_language else None
        updated_on = repo.get("updated_on")

        watchers = _count(session, repo_path, "watchers")
        forks = _count(session, repo_path, "forks")

        return {
            "stars": watchers,
            "forks": forks,
            "language": language,
            "pushed_at": updated_on,
            "archived": False,
        }
    except requests.RequestException as e:
        return {"error": str(e)}


def _count(session: requests.Session, repo_path: str, endpoint: str) -> Optional[int]:
    """Return the ``size`` (total count) for watchers or forks, or None."""
    try:
        resp = _get(session, f"{API_BASE}/{repo_path}/{endpoint}", params={"pagelen": 0})
        if resp.status_code == 200:
            return resp.json().get("size")
    except requests.RequestException:
        pass
    return None


if __name__ == "__main__":
    import json
    import sys

    for path in sys.argv[1:]:
        rp = extract_bitbucket_repo(path) or path
        print(rp, json.dumps(get_repo_metrics(rp)))
