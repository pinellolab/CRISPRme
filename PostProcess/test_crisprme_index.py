#!/usr/bin/env python
"""Tests for PostProcess/crisprme_index.py.

Two groups:

1. Command-construction unit tests (no execution, no network): verify the
   crispritz argv assembled by ``index_genome_command`` / ``add_variants_command``,
   ``read_true_pam`` PAM parsing, and that a missing ``crispritz.py`` yields a
   clean SystemExit from the ``build`` command.

2. A LIVE upload/download/list round-trip against the public HF dataset, isolated
   in a throwaway ``indexes/_smoke_<rand>/`` prefix and cleaned up afterwards.
   The random suffix is derived from os.urandom / os.getpid so parallel runs do
   not collide. If HF write access is unavailable, the round-trip is skipped
   (not failed).

Run all:   python -m pytest PostProcess/test_crisprme_index.py -v
Skip live: CRISPRME_SKIP_HF_TESTS=1 python -m pytest PostProcess/test_crisprme_index.py -v
"""

import binascii
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import crisprme_index as ci  # noqa: E402


# --------------------------------------------------------------------------- #
# command construction (no execution)
# --------------------------------------------------------------------------- #
def test_read_true_pam_ngg(tmp_path):
    pam = tmp_path / "ngg.txt"
    pam.write_text("NNNNNNNNNNNNNNNNNNNNNGG 3\n")
    assert ci.read_true_pam(str(pam)) == "NGG"


def test_read_true_pam_tttv(tmp_path):
    pam = tmp_path / "tttv.txt"
    pam.write_text("TTTVNNNNNNNNNNNNNNNNNNNN -4\n")
    assert ci.read_true_pam(str(pam)) == "TTTV"


def test_index_genome_command():
    cmd = ci.index_genome_command(
        "hg38", "/data/Genomes/hg38", "/data/PAMs/ngg.txt", 2, 8
    )
    assert cmd == [
        "crispritz.py",
        "index-genome",
        "hg38",
        "/data/Genomes/hg38/",  # trailing slash added
        "/data/PAMs/ngg.txt",
        "-bMax",
        "2",
        "-th",
        "8",
    ]


def test_index_genome_command_no_double_slash():
    # trailing slash on input dir must not become a double slash
    cmd = ci.index_genome_command(
        "hg38", "/data/Genomes/hg38/", "/data/PAMs/ngg.txt", 1, 4
    )
    assert cmd[3] == "/data/Genomes/hg38/"


def test_index_genome_command_custom_launcher():
    cmd = ci.index_genome_command(
        "hg38", "/g", "/p.txt", 3, 2, crispritz="/opt/crispritz.py"
    )
    assert cmd[0] == "/opt/crispritz.py"


def test_add_variants_command():
    cmd = ci.add_variants_command("/data/VCFs/1000G", "/data/Genomes/hg38")
    assert cmd == [
        "crispritz.py",
        "add-variants",
        "/data/VCFs/1000G/",
        "/data/Genomes/hg38/",
        "true",
    ]


def test_build_missing_crispritz_clean_error(tmp_path, monkeypatch):
    """build must fail with a clear SystemExit when crispritz.py is absent."""
    genome = tmp_path / "hg38"
    genome.mkdir()
    (genome / "chr22.fa").write_text(">chr22\nACGT\n")
    pam = tmp_path / "ngg.txt"
    pam.write_text("NNNNNNNNNNNNNNNNNNNNNGG 3\n")

    # ensure the launcher is reported as missing regardless of the host machine
    monkeypatch.setattr(ci, "_crispritz_available", lambda c="crispritz.py": False)

    args = ci.build_parser().parse_args(
        [
            "build",
            "--genome", str(genome),
            "--pam", str(pam),
            "--name", "hg38",
            "--work-dir", str(tmp_path),
        ]
    )
    with pytest.raises(SystemExit) as exc:
        args.func(args)
    msg = str(exc.value)
    assert "not found on PATH" in msg
    assert "crispritz" in msg


def test_human_size():
    assert ci._human_size(512) == "512 B"
    assert ci._human_size(2048) == "2.0 KB"
    assert ci._human_size(5 * 1024 * 1024) == "5.0 MB"


# --------------------------------------------------------------------------- #
# LIVE round-trip against the public HF dataset (isolated + cleaned up)
# --------------------------------------------------------------------------- #
def _smoke_suffix() -> str:
    """Random-ish suffix from env override / urandom / pid to avoid collisions."""
    env = os.environ.get("CRISPRME_SMOKE_SUFFIX")
    if env:
        return env
    rnd = binascii.hexlify(os.urandom(4)).decode()
    return f"{os.getpid()}_{rnd}"


@pytest.mark.skipif(
    os.environ.get("CRISPRME_SKIP_HF_TESTS") == "1",
    reason="CRISPRME_SKIP_HF_TESTS=1",
)
def test_hf_upload_download_list_roundtrip(tmp_path):
    hf = pytest.importorskip("huggingface_hub")
    from huggingface_hub import HfApi
    from huggingface_hub.utils import HfHubHTTPError

    repo = ci.HF_DATA_REPO
    api = HfApi()

    # verify we have a usable token with write access; otherwise skip
    try:
        api.whoami()
    except Exception as e:
        pytest.skip(f"no usable HF token: {e}")

    suffix = _smoke_suffix()
    name = f"_smoke_{suffix}"
    path_in_repo = f"{ci.HF_INDEX_PREFIX}/{name}"

    # ---- create a tiny staging folder ------------------------------------- #
    src = tmp_path / "indexes-build" / name
    (src / "genome_library").mkdir(parents=True)
    files = {
        "manifest.json": ('{"name": "%s", "smoke": true}\n' % name).encode(),
        "genome_library/a.txt": b"hello index\n",
        "genome_library/b.bin": bytes(range(256)) + os.urandom(64),
    }
    for rel, content in files.items():
        p = src / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(content)

    class NS:
        pass

    up = NS()
    up.name, up.path, up.repo = name, str(src), repo

    # ---- upload (skip-with-message if write denied) ----------------------- #
    try:
        ci.upload(up)
    except (HfHubHTTPError, PermissionError) as e:
        pytest.skip(f"HF write access unavailable: {e}")
    except Exception as e:
        if "403" in str(e) or "401" in str(e) or "authorized" in str(e).lower():
            pytest.skip(f"HF write access unavailable: {e}")
        raise

    try:
        # ---- list: our smoke prefix must appear --------------------------- #
        repo_files = api.list_repo_files(repo_id=repo, repo_type="dataset")
        assert any(f.startswith(path_in_repo + "/") for f in repo_files), (
            f"{path_in_repo} not found after upload"
        )

        # ---- download + byte-identical round-trip ------------------------- #
        dest = tmp_path / "dl"
        dl = NS()
        dl.name, dl.dest, dl.repo = name, str(dest), repo
        ci.download(dl)

        # genome_library/ must be promoted to the dest root; bytes must match
        got_a = (dest / "genome_library" / "a.txt").read_bytes()
        got_b = (dest / "genome_library" / "b.bin").read_bytes()
        assert got_a == files["genome_library/a.txt"]
        assert got_b == files["genome_library/b.bin"]
    finally:
        # ---- cleanup: always delete the throwaway prefix ------------------ #
        try:
            api.delete_folder(
                path_in_repo=path_in_repo, repo_id=repo, repo_type="dataset"
            )
            print(f"[test] deleted smoke prefix {path_in_repo}")
        except Exception as e:
            print(f"[test] WARNING: failed to delete {path_in_repo}: {e}")

    # confirm deletion took effect
    repo_files_after = api.list_repo_files(repo_id=repo, repo_type="dataset")
    assert not any(f.startswith(path_in_repo + "/") for f in repo_files_after), (
        f"smoke prefix {path_in_repo} still present after cleanup"
    )
