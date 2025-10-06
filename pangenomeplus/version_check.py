"""Version checking against GitHub releases."""

import json
import urllib.request
from typing import Optional

from pangenomeplus import __version__


def check_for_updates() -> Optional[str]:
    """Check GitHub for newer version of PanGenomePlus.

    Returns:
        Version string of latest release if newer than current, None otherwise.
        Returns None on any error (network, parsing, etc.) - fails silently.
    """
    try:
        # GitHub API endpoint for latest release
        url = "https://api.github.com/repos/mol-evol/pangenomeplus/releases/latest"

        # Non-blocking request with 2 second timeout
        req = urllib.request.Request(url)
        req.add_header("Accept", "application/vnd.github.v3+json")

        with urllib.request.urlopen(req, timeout=2) as response:
            data = json.loads(response.read().decode())

            # Extract version from tag_name (e.g., "v0.2.0" -> "0.2.0")
            latest_tag = data.get("tag_name", "")
            latest_version = latest_tag.lstrip("v")

            # Compare versions (simple string comparison works for semantic versioning)
            if latest_version and latest_version > __version__:
                return latest_version

    except (urllib.error.URLError, urllib.error.HTTPError, json.JSONDecodeError, KeyError, TimeoutError):
        # Fail silently - don't annoy user with network/API errors
        pass

    return None
