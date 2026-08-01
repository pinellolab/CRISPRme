---
name: release-crisprme
description: Cut a new CRISPRme release and keep GitHub, Docker, and Bioconda in sync. Use when the user wants to release, publish, tag, or bump the version of CRISPRme; update the Bioconda recipe; or update the changelog for a new version. Covers the version-consistency pre-flight, CHANGELOG.md, the GitHub release/tag, the Bioconda recipe update (autobump or manual PR), and final version verification.
---

# Release CRISPRme

Cut a CRISPRme release and keep three artifacts in sync:

```
git tag vX.Y.Z ──▶ GitHub release tarball ──▶ Bioconda meta.yaml (sha256 of that tarball)
        │                                              │
        └── version pinned in crisprme.py + Dockerfile ┘
```

Version-sync points (ALL must equal the same `X.Y.Z`):
- `crisprme.py` → `version = "X.Y.Z"` (also sets `__version__`)
- `Dockerfile` → `ARG crisprme_version=X.Y.Z`
- git tag / GitHub release → `vX.Y.Z`
- Bioconda `recipes/crisprme/meta.yaml` → `{% set version = "X.Y.Z" %}` + fresh `sha256` + `build.number: 0`

Bioconda source URL is `https://github.com/pinellolab/CRISPRme/archive/refs/tags/v{{ version }}.tar.gz`,
so the recipe pulls the **GitHub release tarball** — its sha256 is the hash of that exact tarball.

Run all commands from the CRISPRme repo root. Replace `X.Y.Z` with the real version (no leading `v` except in tags/URLs).

---

## Step 0 — Pre-flight: version-consistency check

Confirm where the repo currently stands and that nothing is half-bumped.

```bash
echo "crisprme.py: $(grep -E '^version *= *' crisprme.py)"
echo "Dockerfile:  $(grep -E '^ARG crisprme_version=' Dockerfile)"
echo "latest tag:  $(git tag --sort=-v:refname | head -1)"
echo "latest release: $(gh release list --repo pinellolab/CRISPRme -L1)"
curl -sL https://raw.githubusercontent.com/bioconda/bioconda-recipes/master/recipes/crisprme/meta.yaml | grep -E 'set version|sha256|number:'
```

Expected before a release: `crisprme.py`, `Dockerfile`, latest tag, and the recipe all show the **previous** version, consistent with each other. If they disagree (e.g. Dockerfile lags crisprme.py), fix the lag first. **Do not proceed with inconsistent state.**

## Step 1 — Bump in-repo version pins (helper script)

The helper edits ONLY `crisprme.py` and `Dockerfile`, is idempotent, and prints the recipe + changelog scaffolds. The sha256 fetch only works AFTER the tag exists, so run it now with `--no-download`, then again after Step 3 to get the hash.

```bash
python scripts/prepare_release.py X.Y.Z --no-download
git --no-pager diff --stat        # must show ONLY crisprme.py + Dockerfile
git --no-pager diff crisprme.py Dockerfile
```

Rollback if wrong: `git checkout crisprme.py Dockerfile`.

## Step 2 — Update CHANGELOG.md

Move items from `## [Unreleased]` into a new dated section, using the script's scaffold as a template. Keep-a-Changelog categories: Added / Changed / Deprecated / Removed / Fixed / Security.

```bash
python scripts/prepare_release.py X.Y.Z --no-download | sed -n '/CHANGELOG.md scaffold/,/releases\/tag/p'
```

Edit `CHANGELOG.md`:
- Add `## [X.Y.Z] - YYYY-MM-DD` with the curated entries.
- Leave a fresh empty `## [Unreleased]` at the top.
- Update the link-reference footer:
  `[X.Y.Z]: https://github.com/pinellolab/CRISPRme/releases/tag/vX.Y.Z`
  and repoint `[Unreleased]` to `.../compare/vX.Y.Z...HEAD`.

## Step 3 — Commit, then create the GitHub release + tag

Commit the version bump + changelog on the release branch/`main` (per repo policy), push, then create the release. `gh release create` creates the annotated tag `vX.Y.Z` for you.

```bash
git add crisprme.py Dockerfile CHANGELOG.md
git commit -m "Release vX.Y.Z"
git push origin HEAD

# Notes come straight from the changelog section you just wrote:
awk '/^## \[X.Y.Z\]/{f=1;next} /^## \[/{f=0} f' CHANGELOG.md > /tmp/crisprme_notes.md

gh release create vX.Y.Z \
  --repo pinellolab/CRISPRme \
  --title "CRISPRme vX.Y.Z" \
  --notes-file /tmp/crisprme_notes.md \
  --target main            # or the branch you released from
```

Verify the tarball is now published:

```bash
curl -sIL https://github.com/pinellolab/CRISPRme/archive/refs/tags/vX.Y.Z.tar.gz | head -1  # expect 200
```

Rollback: `gh release delete vX.Y.Z --repo pinellolab/CRISPRme --cleanup-tag` (only if the release must be withdrawn before Bioconda picks it up).

## Step 4 — Update the Bioconda recipe

First compute the sha256 of the now-published tarball (this is the recipe's `source.sha256`):

```bash
python scripts/prepare_release.py X.Y.Z    # no --no-download; prints sha256 + meta.yaml diff
# or directly:
curl -sL https://github.com/pinellolab/CRISPRme/archive/refs/tags/vX.Y.Z.tar.gz | shasum -a 256
```

You have two routes. **Prefer autobump; fall back to a manual PR.**

### Route A — Autobump bot (preferred, low-effort)

Bioconda runs an autobump bot that periodically scans GitHub releases and opens
`Update crisprme to X.Y.Z by BiocondaBot` PRs on its own. Wait ~a day, then check:

```bash
gh pr list --repo bioconda/bioconda-recipes --search "crisprme in:title" --state all -L 5
```

If a PR exists, review that it has the right version + sha256 + `build.number: 0`, wait for CI/lint (green), then a Bioconda member or you (once approved) comment on the PR:

```
@BiocondaBot please merge
```

To nudge the bot to scan immediately instead of waiting, comment on any recipe PR/issue:

```
@BiocondaBot autobump crisprme
```

### Route B — Manual PR (reliable, when you don't want to wait)

```bash
# one-time: fork bioconda/bioconda-recipes to your account
gh repo fork bioconda/bioconda-recipes --clone --remote
cd bioconda-recipes
git checkout master && git pull upstream master
git checkout -b crisprme-X.Y.Z
```

Edit `recipes/crisprme/meta.yaml`:
- `{% set version = "X.Y.Z" %}`  ← bump
- `source.sha256:` ← the hash from above
- `build.number: 0`  ← RESET to 0 (only INCREMENT, keeping version, for a same-version rebuild)
- Do NOT edit `source.url` — it uses `{{ version }}` and resolves automatically.

```bash
git commit -am "Update crisprme to X.Y.Z"
git push -u origin crisprme-X.Y.Z
gh pr create --repo bioconda/bioconda-recipes \
  --title "Update crisprme to X.Y.Z" \
  --body "Bump crisprme to vX.Y.Z. sha256 recomputed from the GitHub release tarball. build.number reset to 0."
```

Useful @BiocondaBot PR comments:
- `@BiocondaBot please update` — sync your PR branch with upstream master.
- `@BiocondaBot please add label` — request the `please review & merge` label.
- `@BiocondaBot please fetch artifacts` — links to the built test packages/containers.
- `@BiocondaBot please merge` — after lint + CI are green and the PR is approved, uploads packages and merges.

Wait for CI/lint to pass. Fix any lint findings the bot reports.

## Step 5 — Verify GitHub ↔ Bioconda ↔ Docker match

Do this AFTER the Bioconda PR merges and the package is on the channel (can take an hour+).

```bash
# GitHub / in-repo
grep -E '^version *= *' crisprme.py
grep -E '^ARG crisprme_version=' Dockerfile
git tag --sort=-v:refname | head -1

# Bioconda published package
conda search -c bioconda crisprme | tail -3          # or: mamba search -c bioconda crisprme

# Optional: build the Docker image and check it installs the new version
# (Dockerfile installs crisprme=$crisprme_version from bioconda, so this fails
#  until Bioconda has published X.Y.Z)
docker build -t crisprme:X.Y.Z . && \
  docker run --rm crisprme:X.Y.Z crisprme.py --version   # expect vX.Y.Z
```

All four (crisprme.py, Dockerfile ARG, git tag, bioconda version) must read `X.Y.Z`. Done.

---

## Failure / rollback quick reference
- Wrong in-repo bump: `git checkout crisprme.py Dockerfile CHANGELOG.md`.
- Wrong/premature GitHub release: `gh release delete vX.Y.Z --repo pinellolab/CRISPRme --cleanup-tag` and redo.
- sha256 mismatch on the Bioconda PR: GitHub can regenerate archive tarballs; recompute the sha256 from the live tarball and update `meta.yaml`, then push.
- Deps changed but version unchanged: keep `version`, INCREMENT `build.number` (do not reset to 0).
- Docker build fails on `crisprme=$crisprme_version`: the Bioconda package for that version isn't published yet — wait, then rebuild.
