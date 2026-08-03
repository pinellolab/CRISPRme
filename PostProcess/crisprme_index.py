#!/usr/bin/env python
"""crisprme_index.py — build, upload, download and list precomputed CRISPRme
indexes on HuggingFace.

CRISPRme's two slowest one-time stages for a given (genome, PAM, VCF-set) are:

  1. Genome enrichment  — ``crispritz.py add-variants`` produces the
     variant-enriched FASTAs, the per-chromosome SNP dictionaries
     (``my_dict_<chr>.json``) and the fake-indel FASTAs.
  2. TST index build    — ``crispritz.py index-genome`` produces the
     ``genome_library/<true_pam>_<bMax>_<name>`` trie indexes (the pipeline
     builds the index for bMax, bMax+1 and 1).

These artifacts are fully reusable: once built for the default reference they
never change. This tool packages them and ships them through the CRISPRme
HuggingFace dataset so end users can download a precomputed index and skip both
stages entirely.

HuggingFace layout
------------------
Every index lives under the ``indexes/<name>/`` prefix of the HF dataset
(default ``lucapinello/crisprme-data``, overridable via
``CRISPRME_HF_DATA_REPO``)::

    indexes/<name>/
        manifest.json              # provenance: genome, pam, bMax, vcf, ...
        genome_library/<true_pam>_<bMax>_<name>/...        # reference TST index
        genome_library/<true_pam>_<bMax>_<name>+<vcf>/...  # variant TST index (if --vcf)
        Genomes/<name>+<vcf>/...                            # enriched genome (if --vcf)
        Genomes/<name>+<vcf>_INDELS/...                     # fake-indel FASTAs (if --vcf)
        Dictionaries/...                                    # variant dictionaries (if --vcf)

Download target layout (for a CRISPRme working directory)
---------------------------------------------------------
After ``download``, the ``genome_library/``, ``Genomes/`` and ``Dictionaries/``
folders are placed directly at the root of the destination directory so a
subsequent ``crisprme.py complete-search`` (run with that directory as the
current working directory) finds them without rebuilding.

Subcommands
-----------
  build     orchestrate ``crispritz.py`` to generate an index staging folder
  upload    push a staging folder to ``indexes/<name>/`` on the HF dataset
  download  pull ``indexes/<name>/`` and lay it out for a CRISPRme search
  list      list the precomputed indexes available on the HF dataset
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from typing import List, Optional

# Reuse the HF dataset plumbing shipped in PostProcess/utils.py (stacked on the
# hf-downloads branch). HF_DATA_REPO honors the CRISPRME_HF_DATA_REPO env var.
try:  # when imported as part of the PostProcess package
    from .utils import HF_DATA_REPO
except Exception:  # when run as a standalone script (python crisprme_index.py)
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from utils import HF_DATA_REPO  # type: ignore

# prefix under which all precomputed indexes live on the HF dataset
HF_INDEX_PREFIX = "indexes"


# --------------------------------------------------------------------------- #
# lazy huggingface_hub import
# --------------------------------------------------------------------------- #
def _require_hf():
    """Import ``huggingface_hub`` lazily, raising a clear install hint if missing."""
    try:
        import huggingface_hub  # noqa: F401
    except ImportError as e:  # pragma: no cover - trivial guard
        raise SystemExit(
            "huggingface_hub is required for this command.\n"
            "Install it with:  pip install huggingface_hub"
        ) from e
    return huggingface_hub


# --------------------------------------------------------------------------- #
# helpers shared with the crispritz pipeline
# --------------------------------------------------------------------------- #
def read_true_pam(pam_file: str) -> str:
    """Extract the ``true_pam`` string from a CRISPRme PAM file.

    A PAM file has a single line ``<full_seq> <pos>`` where ``full_seq`` is the
    guide+PAM sequence and ``pos`` is the signed PAM length. This mirrors the
    bash logic in ``submit_job_automated_new_multiple_vcfs.sh``::

        pos > 0  -> last  ``pos``  characters (e.g. NGG for SpCas9)
        pos < 0  -> first ``-pos`` characters (e.g. TTTV for Cas12a)

    Args:
        pam_file: Path to the PAM definition file.

    Returns:
        The true PAM substring (used in the index folder name).
    """
    with open(pam_file) as fh:
        line = fh.readline().strip()
    fullseqpam, pos_str = line.split()
    pos = int(pos_str)
    if pos > 0:
        return fullseqpam[len(fullseqpam) - pos :]
    # pos < 0: first -pos characters (bash ${fullseqpam:0:-$pos})
    return fullseqpam[:-pos] if pos < 0 else fullseqpam


def index_genome_command(
    ref_name: str,
    ref_dir: str,
    pam_file: str,
    bmax: int,
    thread: int,
    crispritz: str = "crispritz.py",
) -> List[str]:
    """Build the ``crispritz.py index-genome`` argv for a single bMax.

    Mirrors the pipeline call::

        crispritz.py index-genome <ref_name> <ref_dir>/ <pam_file> -bMax <bMax> -th <n>

    Args:
        ref_name: Genome name (index folder is ``<true_pam>_<bMax>_<ref_name>``).
        ref_dir: Directory holding the per-chromosome FASTAs (trailing ``/`` added).
        pam_file: PAM definition file.
        bmax: Maximum bulge size for this index.
        thread: Number of threads.
        crispritz: crispritz launcher (overridable for testing).

    Returns:
        The command as an argv list.
    """
    ref_dir = ref_dir.rstrip("/") + "/"  # crispritz expects a trailing slash
    return [
        crispritz,
        "index-genome",
        ref_name,
        ref_dir,
        pam_file,
        "-bMax",
        str(bmax),
        "-th",
        str(thread),
    ]


def add_variants_command(
    vcf_dir: str,
    ref_dir: str,
    crispritz: str = "crispritz.py",
) -> List[str]:
    """Build the ``crispritz.py add-variants`` argv.

    Mirrors the pipeline call::

        crispritz.py add-variants <vcf_dir>/ <ref_dir>/ true

    Args:
        vcf_dir: Directory of bgzipped VCFs (trailing ``/`` added).
        ref_dir: Reference genome FASTA directory (trailing ``/`` added).
        crispritz: crispritz launcher (overridable for testing).

    Returns:
        The command as an argv list.
    """
    vcf_dir = vcf_dir.rstrip("/") + "/"
    ref_dir = ref_dir.rstrip("/") + "/"
    return [crispritz, "add-variants", vcf_dir, ref_dir, "true"]


def _crispritz_available(crispritz: str = "crispritz.py") -> bool:
    """Return True if the crispritz launcher is on PATH."""
    return shutil.which(crispritz) is not None


def _crispritz_version(crispritz: str = "crispritz.py") -> Optional[str]:
    """Best-effort crispritz version string (None if unavailable)."""
    if not _crispritz_available(crispritz):
        return None
    try:
        out = subprocess.run(
            [crispritz, "--version"],
            capture_output=True,
            text=True,
            timeout=30,
        )
        text = (out.stdout or out.stderr).strip()
        # crispritz.py rejects --version; don't record that error as a "version"
        if not text or "not allowed" in text.lower() or text.lower().startswith(
            "error"
        ):
            return None
        return text
    except Exception:
        return None


def _run(cmd: List[str], cwd: Optional[str] = None) -> None:
    """Log and execute a command, raising SystemExit on a clear failure."""
    print("[crisprme-index] $ " + " ".join(cmd) + (f"   (cwd={cwd})" if cwd else ""))
    if not _crispritz_available(cmd[0]):
        raise SystemExit(
            f"'{cmd[0]}' not found on PATH. crispritz must be installed to run "
            f"'build' (it ships in the CRISPRme conda/Docker environment). "
            f"The command that would have run is:\n    {' '.join(cmd)}"
        )
    result = subprocess.run(cmd, cwd=cwd)
    if result.returncode != 0:
        raise SystemExit(
            f"Command failed (exit {result.returncode}): {' '.join(cmd)}"
        )


# --------------------------------------------------------------------------- #
# build
# --------------------------------------------------------------------------- #
def build(args: argparse.Namespace) -> None:
    """Orchestrate crispritz to produce an index staging folder.

    Populates ``<work-dir>/indexes-build/<name>/`` with ``genome_library/*``,
    the enriched genome + ``Dictionaries/*`` (if ``--vcf``) and a
    ``manifest.json``. crispritz is invoked via subprocess; each command is
    logged. Command construction is factored into the unit-testable helpers
    :func:`index_genome_command` / :func:`add_variants_command`.
    """
    genome_dir = os.path.abspath(args.genome)
    pam_file = os.path.abspath(args.pam)
    name = args.name
    bmax = args.bmax
    thread = args.thread
    work_dir = os.path.abspath(args.work_dir or os.getcwd())
    crispritz = args.crispritz

    if not os.path.isdir(genome_dir):
        raise SystemExit(f"--genome directory not found: {genome_dir}")
    if not os.path.isfile(pam_file):
        raise SystemExit(f"--pam file not found: {pam_file}")

    true_pam = read_true_pam(pam_file)
    pam_name = os.path.basename(pam_file)

    stage = os.path.join(work_dir, "indexes-build", name)
    stage_genome_library = os.path.join(stage, "genome_library")
    stage_genomes = os.path.join(stage, "Genomes")
    stage_dicts = os.path.join(stage, "Dictionaries")
    for d in (stage, stage_genome_library, stage_genomes, stage_dicts):
        os.makedirs(d, exist_ok=True)

    # crispritz writes genome_library/ (and, for add-variants, variants_genome/)
    # relative to the working directory. Run everything from a scratch build
    # root so we don't pollute the caller's cwd, then collect the outputs.
    build_root = os.path.join(stage, "_build")
    os.makedirs(build_root, exist_ok=True)

    # ---- STEP 1: reference TST indexes (bMax, bMax+1, 1) -------------------- #
    # matches idx_folder1/2/3 in the pipeline
    ref_bmax_set = sorted({bmax, bmax + 1, 1})
    print(f"[crisprme-index] building reference indexes for bMax in {ref_bmax_set}")
    for b in ref_bmax_set:
        cmd = index_genome_command(name, genome_dir, pam_file, b, thread, crispritz)
        _run(cmd, cwd=build_root)

    enriched_dir = None
    vcf_name = None
    if args.vcf:
        vcf_dir = os.path.abspath(args.vcf)
        if not os.path.isdir(vcf_dir):
            raise SystemExit(f"--vcf directory not found: {vcf_dir}")
        vcf_name = os.path.basename(vcf_dir.rstrip("/"))

        # ---- STEP 2: enrich genome ---------------------------------------- #
        print("[crisprme-index] enriching genome (add-variants)")
        _run(add_variants_command(vcf_dir, genome_dir, crispritz), cwd=build_root)

        # crispritz's add-variants names the enriched dir after the GENOME DIR
        # basename (not --name): "<genome_basename>_enriched".
        genome_basename = os.path.basename(os.path.normpath(genome_dir))
        snps_genome = os.path.join(build_root, "variants_genome", "SNPs_genome")
        enriched_src = os.path.join(snps_genome, f"{genome_basename}_enriched")
        enriched_name = f"{name}+{vcf_name}"
        enriched_dir = os.path.join(build_root, "Genomes", enriched_name)
        os.makedirs(os.path.dirname(enriched_dir), exist_ok=True)
        if os.path.isdir(enriched_src):
            shutil.move(enriched_src, enriched_dir)
        # collect SNP dictionaries + fake-indel FASTAs
        if os.path.isdir(snps_genome):
            for entry in os.listdir(snps_genome):
                src = os.path.join(snps_genome, entry)
                if entry.endswith(".json"):
                    shutil.move(src, os.path.join(stage_dicts, entry))
                elif entry.startswith("fake"):
                    indels_dir = os.path.join(
                        build_root, "Genomes", f"{enriched_name}_INDELS"
                    )
                    os.makedirs(indels_dir, exist_ok=True)
                    shutil.move(src, os.path.join(indels_dir, entry))

        # ---- STEP 3: index the enriched genome ---------------------------- #
        if os.path.isdir(enriched_dir):
            print("[crisprme-index] building variant indexes")
            for b in ref_bmax_set:
                cmd = index_genome_command(
                    enriched_name, enriched_dir, pam_file, b, thread, crispritz
                )
                _run(cmd, cwd=build_root)

    # ---- collect artifacts into the staging layout ------------------------- #
    src_genome_library = os.path.join(build_root, "genome_library")
    if os.path.isdir(src_genome_library):
        for entry in os.listdir(src_genome_library):
            src = os.path.join(src_genome_library, entry)
            dst = os.path.join(stage_genome_library, entry)
            if not os.path.exists(dst):
                shutil.move(src, dst)
    src_genomes = os.path.join(build_root, "Genomes")
    if os.path.isdir(src_genomes):
        for entry in os.listdir(src_genomes):
            src = os.path.join(src_genomes, entry)
            dst = os.path.join(stage_genomes, entry)
            if not os.path.exists(dst):
                shutil.move(src, dst)

    # drop empty staging folders so we don't upload empty dirs
    for d in (stage_genomes, stage_dicts):
        if os.path.isdir(d) and not os.listdir(d):
            os.rmdir(d)
    # remove the crispritz scratch working dir entirely — the real artifacts have
    # been moved out to the staging layout above, and it must not be uploaded
    shutil.rmtree(build_root, ignore_errors=True)

    # ---- manifest ---------------------------------------------------------- #
    manifest = {
        "name": name,
        "genome": os.path.basename(genome_dir.rstrip("/")),
        "pam": pam_name,
        "pam_name": true_pam,
        "bMax": bmax,
        "vcf": vcf_name,
        "thread": thread,
        "crispritz_version": _crispritz_version(crispritz),
        "created_by": os.environ.get("USER") or os.environ.get("USERNAME"),
        "created_at": datetime.now(timezone.utc).isoformat(),
        "note": args.note,
    }
    manifest_path = os.path.join(stage, "manifest.json")
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2)
    print(f"[crisprme-index] wrote manifest -> {manifest_path}")
    print(f"[crisprme-index] staging folder ready -> {stage}")
    print(f"[crisprme-index] next: crisprme_index.py upload --name {name}")


# --------------------------------------------------------------------------- #
# upload
# --------------------------------------------------------------------------- #
def _mirror_tree(src: str, dst: str) -> None:
    """Recreate ``src`` at ``dst`` using hardlinks (no data copy) where possible.

    Falls back to a real copy when hardlinking is unsupported (e.g. across
    filesystems). Used to give ``upload_large_folder`` a scratch tree rooted at
    ``indexes/<name>/`` without duplicating multi-GB index files on disk.
    """
    os.makedirs(dst, exist_ok=True)
    for entry in os.listdir(src):
        s = os.path.join(src, entry)
        d = os.path.join(dst, entry)
        if os.path.isdir(s):
            _mirror_tree(s, d)
        else:
            try:
                os.link(s, d)
            except OSError:
                shutil.copyfile(s, d)


def upload(args: argparse.Namespace) -> None:
    """Upload a staging folder to ``indexes/<name>/`` on the HF dataset."""
    hf = _require_hf()
    repo = args.repo or HF_DATA_REPO
    path = os.path.abspath(
        args.path or os.path.join(os.getcwd(), "indexes-build", args.name)
    )
    if not os.path.isdir(path):
        raise SystemExit(f"index folder not found: {path}")
    path_in_repo = f"{HF_INDEX_PREFIX}/{args.name}"

    print(f"[crisprme-index] uploading {path} -> {repo}:{path_in_repo}")
    # upload_large_folder is the preferred path for big indexes (multi-worker,
    # resumable). Newer hub versions accept path_in_repo directly; older ones
    # upload folder_path to the repo *root*. To target indexes/<name>/ portably,
    # mirror the staging tree into a scratch root laid out as indexes/<name>/...
    # using hardlinks (same filesystem => no data copy; falls back to copy across
    # filesystems). If upload_large_folder is unavailable or errors, fall back to
    # upload_folder, which supports path_in_repo directly.
    uploaded = False
    if hasattr(hf, "upload_large_folder"):
        import tempfile

        # place scratch alongside the source so hardlinks stay on one filesystem
        scratch = tempfile.mkdtemp(
            prefix=".crisprme_index_upload_", dir=os.path.dirname(path)
        )
        try:
            mirror = os.path.join(scratch, HF_INDEX_PREFIX, args.name)
            os.makedirs(os.path.dirname(mirror), exist_ok=True)
            _mirror_tree(path, mirror)  # hardlink files under indexes/<name>/
            hf.upload_large_folder(
                repo_id=repo,
                repo_type="dataset",
                folder_path=scratch,
            )
            uploaded = True
        except Exception as e:  # pragma: no cover - depends on hub version
            print(
                f"[crisprme-index] upload_large_folder failed ({e}); "
                "falling back to upload_folder"
            )
        finally:
            shutil.rmtree(scratch, ignore_errors=True)
    if not uploaded:
        hf.upload_folder(
            repo_id=repo,
            repo_type="dataset",
            folder_path=path,
            path_in_repo=path_in_repo,
        )
    tree_url = f"https://huggingface.co/datasets/{repo}/tree/main/{path_in_repo}"
    print(f"[crisprme-index] uploaded. Browse: {tree_url}")


# --------------------------------------------------------------------------- #
# download
# --------------------------------------------------------------------------- #
def download(args: argparse.Namespace) -> None:
    """Download ``indexes/<name>/`` and lay it out for a CRISPRme search.

    Files are snapshot-downloaded then the ``genome_library/``, ``Genomes/`` and
    ``Dictionaries/`` folders are placed at the root of ``--dest`` so a CRISPRme
    search run from that directory finds them without rebuilding.
    """
    hf = _require_hf()
    repo = args.repo or HF_DATA_REPO
    dest = os.path.abspath(args.dest or os.getcwd())
    os.makedirs(dest, exist_ok=True)
    path_in_repo = f"{HF_INDEX_PREFIX}/{args.name}"

    print(f"[crisprme-index] downloading {repo}:{path_in_repo}/* -> {dest}")
    snapshot = hf.snapshot_download(
        repo_id=repo,
        repo_type="dataset",
        allow_patterns=f"{path_in_repo}/*",
        local_dir=dest,
    )
    downloaded_root = os.path.join(snapshot, HF_INDEX_PREFIX, args.name)
    if not os.path.isdir(downloaded_root):
        raise SystemExit(
            f"nothing downloaded for index '{args.name}'. "
            f"Run 'crisprme_index.py list' to see available indexes."
        )

    # promote genome_library/ Genomes/ Dictionaries/ to the dest root so a
    # CRISPRme search (running from dest) finds them.
    for folder in ("genome_library", "Genomes", "Dictionaries"):
        src = os.path.join(downloaded_root, folder)
        if not os.path.isdir(src):
            continue
        dst = os.path.join(dest, folder)
        os.makedirs(dst, exist_ok=True)
        for entry in os.listdir(src):
            s = os.path.join(src, entry)
            d = os.path.join(dst, entry)
            if os.path.exists(d):
                print(f"[crisprme-index] skip existing {d}")
                continue
            # copy (snapshot may be symlinks into the HF cache); dereference
            if os.path.isdir(s):
                shutil.copytree(s, d, symlinks=False)
            else:
                shutil.copyfile(s, d)
    manifest = os.path.join(downloaded_root, "manifest.json")
    if os.path.isfile(manifest):
        shutil.copyfile(manifest, os.path.join(dest, f"manifest.{args.name}.json"))
    print(f"[crisprme-index] index '{args.name}' ready under {dest}")
    print(
        "[crisprme-index] run your CRISPRme search with this directory as the "
        "current working directory to skip enrichment + indexing."
    )


# --------------------------------------------------------------------------- #
# list
# --------------------------------------------------------------------------- #
def list_indexes(args: argparse.Namespace) -> None:
    """List the ``indexes/<name>/`` prefixes present on the HF dataset."""
    hf = _require_hf()
    repo = args.repo or HF_DATA_REPO
    api = hf.HfApi()
    try:
        files = api.list_repo_files(repo_id=repo, repo_type="dataset")
    except Exception as e:
        raise SystemExit(f"failed to list repo '{repo}': {e}")

    # collect the top-level names under indexes/
    names = set()
    for f in files:
        parts = f.split("/")
        if len(parts) >= 2 and parts[0] == HF_INDEX_PREFIX:
            names.add(parts[1])
    if not names:
        print(f"[crisprme-index] no precomputed indexes found in {repo}")
        return

    # best-effort sizes (cheap: from repo tree metadata if available)
    sizes = {}
    try:
        for item in api.list_repo_tree(
            repo_id=repo, repo_type="dataset", path_in_repo=HF_INDEX_PREFIX,
            recursive=True,
        ):
            size = getattr(item, "size", None)
            path = getattr(item, "path", "")
            parts = path.split("/")
            if size and len(parts) >= 2 and parts[0] == HF_INDEX_PREFIX:
                sizes[parts[1]] = sizes.get(parts[1], 0) + size
    except Exception:
        pass  # sizes are best-effort only

    print(f"Precomputed indexes on {repo} (indexes/<name>/):")
    for name in sorted(names):
        if name in sizes:
            print(f"  {name}    ({_human_size(sizes[name])})")
        else:
            print(f"  {name}")


def _human_size(n: int) -> str:
    """Format a byte count as a human-readable string."""
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024 or unit == "TB":
            return f"{n:.1f} {unit}" if unit != "B" else f"{n} B"
        n /= 1024.0
    return f"{n} B"


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def build_parser() -> argparse.ArgumentParser:
    """Construct the argparse CLI."""
    parser = argparse.ArgumentParser(
        prog="crisprme_index.py",
        description=(
            "Build, upload, download and list precomputed CRISPRme indexes on "
            "HuggingFace. Precomputed indexes let end users skip the two slowest "
            "one-time stages: genome enrichment (add-variants) and TST index "
            "building (index-genome).\n\n"
            "HF layout: indexes/<name>/{manifest.json, genome_library/, "
            "Genomes/, Dictionaries/}.\n"
            "download target layout: genome_library/, Genomes/ and Dictionaries/ "
            "are placed at the root of --dest so a CRISPRme search run from that "
            "directory finds them."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # build
    p_build = sub.add_parser(
        "build",
        help="orchestrate crispritz to generate an index staging folder",
        description=(
            "Build a precomputed index into <work-dir>/indexes-build/<name>/. "
            "Runs 'crispritz.py index-genome' for bMax, bMax+1 and 1; if --vcf "
            "is given, also runs 'crispritz.py add-variants' and indexes the "
            "enriched genome. Writes a manifest.json."
        ),
    )
    p_build.add_argument("--genome", required=True, help="reference FASTA directory")
    p_build.add_argument("--pam", required=True, help="PAM definition file")
    p_build.add_argument("--name", required=True, help="index / genome name")
    p_build.add_argument("--vcf", help="directory of bgzipped VCFs (enable enrichment)")
    p_build.add_argument("--bmax", type=int, default=2, help="max bulge size (default 2)")
    p_build.add_argument("--thread", type=int, default=4, help="threads (default 4)")
    p_build.add_argument("--work-dir", help="staging root (default cwd)")
    p_build.add_argument("--note", default="", help="free-text note for the manifest")
    p_build.add_argument(
        "--crispritz", default="crispritz.py", help="crispritz launcher (default crispritz.py)"
    )
    p_build.set_defaults(func=build)

    # upload
    p_up = sub.add_parser(
        "upload", help="upload a staging folder to indexes/<name>/ on the HF dataset"
    )
    p_up.add_argument("--name", required=True, help="index name (indexes/<name>/)")
    p_up.add_argument("--path", help="folder to upload (default <cwd>/indexes-build/<name>/)")
    p_up.add_argument("--repo", help=f"HF dataset repo (default {HF_DATA_REPO})")
    p_up.set_defaults(func=upload)

    # download
    p_dl = sub.add_parser(
        "download", help="download indexes/<name>/ and lay it out for a CRISPRme search"
    )
    p_dl.add_argument("--name", required=True, help="index name (indexes/<name>/)")
    p_dl.add_argument("--dest", help="destination directory (default cwd)")
    p_dl.add_argument("--repo", help=f"HF dataset repo (default {HF_DATA_REPO})")
    p_dl.set_defaults(func=download)

    # list
    p_ls = sub.add_parser("list", help="list precomputed indexes on the HF dataset")
    p_ls.add_argument("--repo", help=f"HF dataset repo (default {HF_DATA_REPO})")
    p_ls.set_defaults(func=list_indexes)

    return parser


def main(argv: Optional[List[str]] = None) -> None:
    """CLI entrypoint."""
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
