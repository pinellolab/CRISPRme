# Releasing CRISPRme

This guide explains how a maintainer cuts a new CRISPRme release and keeps the
three published artifacts — the **GitHub release**, the **Docker image**, and the
**Bioconda package** — in sync. For an automated, step-by-step assistant version,
see the `release-crisprme` Claude Code skill (`.claude/skills/release-crisprme/SKILL.md`);
this document is the human-readable counterpart.

## Using the release skill (Claude Code)

Day to day, the easiest way to cut a release is the `release-crisprme` skill from
a Claude Code session opened at the repo root:

1. **Prerequisites (one-time / per-release):** `main` is green and up to date,
   `gh` is authenticated, and the `DOCKERHUB_USERNAME` / `DOCKERHUB_TOKEN`
   repository secrets are set (already configured for `pinellolab/crisprme`, so
   the multi-arch image publishes on the release tag).
2. **Invoke it:** type `/release-crisprme`, or just ask in plain language, e.g.
   *"release CRISPRme 2.1.13"*. The skill then walks the whole flow:
   - pre-flight version-consistency check across `crisprme.py`, `Dockerfile`, the
     git tag, and the Bioconda recipe;
   - move `CHANGELOG.md` `[Unreleased]` into the new version;
   - bump the version everywhere (`scripts/prepare_release.py X.Y.Z`);
   - create the GitHub release/tag `vX.Y.Z` — this triggers the multi-arch Docker
     publish automatically;
   - update the Bioconda recipe (autobump PR or a manual PR) with the tarball
     `sha256` and `build.number: 0`;
   - verify GitHub, Docker, and Bioconda all report `X.Y.Z`.
3. It **pauses before anything irreversible** (creating the public tag/release)
   so you can confirm first.

The rest of this document is the manual, step-by-step version of exactly what the
skill automates — use it if you prefer to run the commands yourself, or if the
skill isn't available.

## The sync model

There is one source of truth for a version — the number `X.Y.Z` — and it must be
identical in four places. The Bioconda recipe additionally pins the **sha256 of the
GitHub release tarball**, which is why the order of operations matters: you cannot
compute that hash until the GitHub release (and therefore the tarball) exists.

```
                    ┌──────────────────────────────────────────────┐
                    │             version  X.Y.Z                    │
                    └──────────────────────────────────────────────┘
                                        │
        ┌───────────────────────────────┼───────────────────────────────┐
        ▼                               ▼                               ▼
  crisprme.py                      Dockerfile                       git tag vX.Y.Z
  version = "X.Y.Z"                ARG crisprme_version=X.Y.Z             │
                                                                         ▼
                                                          GitHub release tarball
                                                 .../archive/refs/tags/vX.Y.Z.tar.gz
                                                                         │
                                                             sha256 of that tarball
                                                                         │
                                                                         ▼
                                              Bioconda recipes/crisprme/meta.yaml
                                                {% set version = "X.Y.Z" %}
                                                source.sha256: <hash>
                                                build.number: 0
```

Note the dependency chain: the Dockerfile installs `crisprme=$crisprme_version`
from the Bioconda channel, so a freshly built Docker image only works once the
Bioconda package for that version has been published. The typical timeline is:

1. Bump in-repo files → 2. GitHub release/tag → 3. Bioconda recipe update (merged)
→ 4. Docker image rebuild works.

## The four version-sync points

| Location | What to change |
|---|---|
| `crisprme.py` | `version = "X.Y.Z"` (this also sets `__version__`, surfaced by `crisprme.py --version`) |
| `Dockerfile` | `ARG crisprme_version=X.Y.Z` |
| Git tag / GitHub release | tag `vX.Y.Z`, notes from the changelog |
| `recipes/crisprme/meta.yaml` (in `bioconda/bioconda-recipes`) | `{% set version = "X.Y.Z" %}`, `source.sha256`, `build.number: 0` |

The Bioconda `source.url` is templated
(`.../archive/refs/tags/v{{ version }}.tar.gz`) and must not be edited by hand —
bumping `version` repoints it automatically.

## Step-by-step

### 1. Pre-flight consistency check

Before starting, confirm the repo, the latest tag, and the live recipe all agree
on the *current* (previous) version. If, for example, the Dockerfile lags behind
`crisprme.py` (this has happened historically), fix that lag first so you start
from a consistent state.

```bash
grep -E '^version *= *' crisprme.py
grep -E '^ARG crisprme_version=' Dockerfile
git tag --sort=-v:refname | head -1
curl -sL https://raw.githubusercontent.com/bioconda/bioconda-recipes/master/recipes/crisprme/meta.yaml \
  | grep -E 'set version|sha256|number:'
```

### 2. Bump the in-repo version pins

Use the helper script, which edits only `crisprme.py` and `Dockerfile`, is
idempotent, and prints the Bioconda-recipe and changelog scaffolds for you:

```bash
python scripts/prepare_release.py X.Y.Z --no-download
git diff crisprme.py Dockerfile     # sanity-check the two-line change
```

(The `--no-download` flag skips the sha256 fetch, which cannot succeed until the
tag exists. You will re-run the script without it in step 5.)

### 3. Update the changelog

CRISPRme follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
Move the accumulated entries from the `[Unreleased]` section into a new
`## [X.Y.Z] - YYYY-MM-DD` section, grouped under Added / Changed / Deprecated /
Removed / Fixed / Security. Leave a fresh empty `[Unreleased]` at the top and
update the link-reference footer at the bottom of the file.

### 4. Create the GitHub release and tag

Commit the version bump and changelog, push, then create the release. `gh`
creates the annotated tag `vX.Y.Z` as part of the release.

```bash
git add crisprme.py Dockerfile CHANGELOG.md
git commit -m "Release vX.Y.Z"
git push origin HEAD

gh release create vX.Y.Z \
  --repo pinellolab/CRISPRme \
  --title "CRISPRme vX.Y.Z" \
  --notes-file <changelog-section-for-X.Y.Z>
```

Confirm the tarball is live (HTTP 200) before touching Bioconda:

```bash
curl -sIL https://github.com/pinellolab/CRISPRme/archive/refs/tags/vX.Y.Z.tar.gz | head -1
```

### 5. Update the Bioconda recipe

Recompute the sha256 from the *published* tarball:

```bash
python scripts/prepare_release.py X.Y.Z        # prints sha256 + the exact meta.yaml diff
# equivalently:
curl -sL https://github.com/pinellolab/CRISPRme/archive/refs/tags/vX.Y.Z.tar.gz | shasum -a 256
```

There are two ways to get this into `bioconda/bioconda-recipes`:

**Autobump (preferred).** Bioconda operates an autobump bot that periodically
watches GitHub releases and opens `Update crisprme to X.Y.Z by BiocondaBot` pull
requests automatically — often within a day of the release. Watch for it:

```bash
gh pr list --repo bioconda/bioconda-recipes --search "crisprme in:title" --state all -L5
```

Verify the auto-PR carries the correct version, sha256, and `build.number: 0`,
wait for CI/lint to go green, and once it is approved, comment
`@BiocondaBot please merge`. To ask the bot to scan immediately rather than wait,
comment `@BiocondaBot autobump crisprme` on a recipe PR/issue.

**Manual PR (when you don't want to wait, or the auto-PR needs correcting).**
Fork `bioconda/bioconda-recipes`, branch from an up-to-date `master`, edit
`recipes/crisprme/meta.yaml` (bump `version`, set the new `sha256`, reset
`build.number` to `0`), and open a PR titled `Update crisprme to X.Y.Z`. CI builds
and tests the recipe automatically. Handy PR comments: `@BiocondaBot please update`
(sync your branch), `@BiocondaBot please add label` (request review), and
`@BiocondaBot please merge` (after green CI + approval).

> Build-number rule: reset `build.number` to `0` for a **new version**. Only
> **increment** it (keeping the same version) when you rebuild an existing version
> because its dependencies changed.

### 6. Verify everything matches

After the Bioconda PR merges and the package publishes (this can take an hour or
more), confirm all four sources report `X.Y.Z`:

```bash
grep -E '^version *= *' crisprme.py
grep -E '^ARG crisprme_version=' Dockerfile
git tag --sort=-v:refname | head -1
conda search -c bioconda crisprme | tail -3
```

Optionally rebuild the Docker image to confirm the full chain:

```bash
docker build -t crisprme:X.Y.Z . && docker run --rm crisprme:X.Y.Z crisprme.py --version
```

## Rollback / troubleshooting

- **Bad in-repo bump:** `git checkout crisprme.py Dockerfile CHANGELOG.md`.
- **Premature/incorrect GitHub release:** `gh release delete vX.Y.Z --repo pinellolab/CRISPRme --cleanup-tag`, then redo.
- **sha256 mismatch in the Bioconda PR:** GitHub occasionally regenerates archive tarballs; recompute the hash from the live tarball and update `meta.yaml`.
- **Docker build fails on `crisprme=$crisprme_version`:** the Bioconda package for that version has not published yet — wait and rebuild.
- **Version lag between files:** always run the step-1 pre-flight; the most common historical drift is the `Dockerfile` ARG trailing `crisprme.py`.
