#!/usr/bin/env python3
"""prepare_release.py - Prepare a CRISPRme release across all version-pinned files.

Given a target version (e.g. 2.1.12), this script:

  1. Bumps the version string in ``crisprme.py`` (``version = "X.Y.Z"``).
  2. Bumps ``ARG crisprme_version=X.Y.Z`` in the ``Dockerfile``.
  3. Downloads the GitHub *release* tarball for the tag ``vX.Y.Z`` and computes
     its sha256 (used by the Bioconda recipe ``source.sha256``).
  4. Prints the exact ``meta.yaml`` changes needed for the Bioconda recipe.
  5. Prints a ``CHANGELOG.md`` scaffold entry for the new version.

The script is IDEMPOTENT: re-running it for the same version makes no further
edits (it reports "already at X.Y.Z"). It only ever touches ``crisprme.py`` and
``Dockerfile``; it never edits the recipe or the changelog for you -- those are
printed for you to apply/verify manually (the recipe lives in a different repo,
and changelog prose needs a human).

It does NOT create git tags, GitHub releases, or Bioconda PRs.

Usage:
    python scripts/prepare_release.py 2.1.12
    python scripts/prepare_release.py 2.1.12 --no-download   # skip sha256 fetch
    python scripts/prepare_release.py 2.1.12 --repo-root /path/to/CRISPRme
"""
from __future__ import annotations

import argparse
import hashlib
import os
import re
import sys
import urllib.request

REPO = "pinellolab/CRISPRme"
TARBALL_URL = "https://github.com/{repo}/archive/refs/tags/v{version}.tar.gz"

VERSION_RE = re.compile(r"^\d+\.\d+\.\d+$")


def die(msg: str) -> "None":
    sys.stderr.write(f"ERROR: {msg}\n")
    sys.exit(1)


def bump_crisprme_py(path: str, version: str) -> bool:
    """Set ``version = "X.Y.Z"`` in crisprme.py. Returns True if changed."""
    with open(path, "r", encoding="utf-8") as fh:
        text = fh.read()
    pat = re.compile(r'^(version\s*=\s*)"(\d+\.\d+\.\d+)"', re.MULTILINE)
    m = pat.search(text)
    if not m:
        die(f"could not find `version = \"...\"` in {path}")
    current = m.group(2)
    if current == version:
        print(f"  crisprme.py       already at {version} (no change)")
        return False
    new_text = pat.sub(lambda _m: f'{_m.group(1)}"{version}"', text, count=1)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(new_text)
    print(f"  crisprme.py       {current} -> {version}")
    return True


def bump_dockerfile(path: str, version: str) -> bool:
    """Set ``ARG crisprme_version=X.Y.Z`` in Dockerfile. Returns True if changed."""
    with open(path, "r", encoding="utf-8") as fh:
        text = fh.read()
    pat = re.compile(r"^(ARG\s+crisprme_version=)(\d+\.\d+\.\d+)", re.MULTILINE)
    m = pat.search(text)
    if not m:
        die(f"could not find `ARG crisprme_version=...` in {path}")
    current = m.group(2)
    if current == version:
        print(f"  Dockerfile        already at {version} (no change)")
        return False
    new_text = pat.sub(lambda _m: f"{_m.group(1)}{version}", text, count=1)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(new_text)
    print(f"  Dockerfile        {current} -> {version}")
    return True


def compute_sha256(version: str) -> str:
    url = TARBALL_URL.format(repo=REPO, version=version)
    print(f"\nFetching release tarball to compute sha256:\n  {url}")
    try:
        with urllib.request.urlopen(url) as resp:  # noqa: S310 (trusted host)
            data = resp.read()
    except Exception as exc:  # noqa: BLE001
        die(
            f"could not download tarball for v{version}: {exc}\n"
            "       (Has the GitHub release/tag been created yet? "
            "The sha256 can only be computed AFTER the tag exists.)"
        )
    digest = hashlib.sha256(data).hexdigest()
    print(f"  bytes: {len(data)}")
    print(f"  sha256: {digest}")
    return digest


def print_meta_yaml_changes(version: str, sha256: str | None) -> None:
    sha_line = sha256 if sha256 else "<run without --no-download to compute>"
    print(
        "\n"
        "==============================================================\n"
        "  Bioconda recipe changes  (recipes/crisprme/meta.yaml)\n"
        "  Repo: https://github.com/bioconda/bioconda-recipes\n"
        "==============================================================\n"
        f'  {{% set version = "{version}" %}}      # bump version\n'
        f"  source:\n"
        f"    sha256: {sha_line}\n"
        f"  build:\n"
        f"    number: 0                           # RESET to 0 for a new version\n"
        "\n"
        "  NOTE: source.url uses {{ version }}, so it auto-resolves to\n"
        f"        .../archive/refs/tags/v{version}.tar.gz -- do not edit it.\n"
        "  NOTE: if this is a REBUILD of the same version (deps changed only),\n"
        "        keep version the same and INCREMENT build.number instead.\n"
    )


def print_changelog_scaffold(version: str) -> None:
    import datetime

    today = datetime.date.today().isoformat()
    print(
        "\n"
        "==============================================================\n"
        "  CHANGELOG.md scaffold  (move items out of [Unreleased])\n"
        "==============================================================\n"
        f"## [{version}] - {today}\n"
        "### Added\n-\n"
        "### Changed\n-\n"
        "### Fixed\n-\n"
        "\n"
        "  Also update the link-reference footer, e.g.:\n"
        f"  [{version}]: https://github.com/{REPO}/releases/tag/v{version}\n"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("version", help="target version, e.g. 2.1.12 (no leading v)")
    parser.add_argument(
        "--repo-root",
        default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        help="path to the CRISPRme repo root (default: parent of scripts/)",
    )
    parser.add_argument(
        "--no-download",
        action="store_true",
        help="skip downloading the tarball / computing sha256",
    )
    args = parser.parse_args()

    version = args.version.lstrip("v")
    if not VERSION_RE.match(version):
        die(f"version must look like X.Y.Z (got {args.version!r})")

    crisprme_py = os.path.join(args.repo_root, "crisprme.py")
    dockerfile = os.path.join(args.repo_root, "Dockerfile")
    for p in (crisprme_py, dockerfile):
        if not os.path.isfile(p):
            die(f"expected file not found: {p}")

    print(f"Preparing CRISPRme release v{version}\n")
    print("Bumping in-repo version pins:")
    changed = False
    changed |= bump_crisprme_py(crisprme_py, version)
    changed |= bump_dockerfile(dockerfile, version)
    if not changed:
        print("  (all in-repo files already at target version)")

    sha256 = None if args.no_download else compute_sha256(version)

    print_meta_yaml_changes(version, sha256)
    print_changelog_scaffold(version)

    print(
        "Next steps (see docs/RELEASING.md or the release-crisprme skill):\n"
        "  1. Review the git diff of crisprme.py + Dockerfile.\n"
        "  2. Move [Unreleased] changelog items into the new version section.\n"
        "  3. Commit, then create the GitHub release/tag vX.Y.Z.\n"
        "  4. Update the Bioconda recipe (autobump PR or manual PR) with the sha256 above.\n"
    )


if __name__ == "__main__":
    main()
