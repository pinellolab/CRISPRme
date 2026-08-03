"""Lightweight tests for the HuggingFace-first reference-data download path.

These tests exercise the NEW ``hf_fetch`` helper and the HF branch of the
CRISPRme download functions end-to-end against the PUBLIC HF dataset
(``lucapinello/crisprme-data``). To keep the suite cheap they only touch small
resources (sample-ID files, a tiny genome FASTA) — never the full genome or the
multi-GB 1000G/HGDP VCF sets.

Run with (network required — dataset is public, no token needed):

    <env>/bin/python3 -m unittest discover -s PostProcess -p 'test_hf_downloads.py' -v

Set ``CRISPRME_SKIP_HF_TESTS=1`` to skip the network-dependent cases.
"""

import importlib
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import utils  # noqa: E402
from utils import compute_md5, hf_fetch, MD5SAMPLES  # noqa: E402

SKIP_NET = os.environ.get("CRISPRME_SKIP_HF_TESTS") == "1"
NET_REASON = "network-dependent HF test skipped (CRISPRME_SKIP_HF_TESTS=1)"

# a tiny genome FASTA in the HF dataset — cheap to fetch for the placement test
SMALL_GENOME_REPO_FILE = "genomes/hg38/chrM.fa"


@unittest.skipIf(SKIP_NET, NET_REASON)
class TestHfFetchSamplesIds(unittest.TestCase):
    """hf_fetch lands the file at the exact CRISPRme target path with matching MD5."""

    def test_1000g_samplesid_lands_and_md5_matches(self):
        with tempfile.TemporaryDirectory() as td:
            target = os.path.join(td, "samplesIDs", "hg38_1000G.samplesID.txt")
            out = hf_fetch("samplesIDs/hg38_1000G.samplesID.txt", target)
            # (a) file lands at the correct CRISPRme target path
            self.assertEqual(out, target)
            self.assertTrue(os.path.isfile(target))
            # (b) md5 matches the expected MD5SAMPLES entry (byte-identical file)
            self.assertEqual(
                compute_md5(target), MD5SAMPLES["samplesIDs.1000G.txt"]
            )

    def test_hgdp_samplesid_lands_and_md5_matches(self):
        with tempfile.TemporaryDirectory() as td:
            target = os.path.join(td, "samplesIDs", "hg38_HGDP.samplesID.txt")
            out = hf_fetch("samplesIDs/hg38_HGDP.samplesID.txt", target)
            self.assertEqual(out, target)
            self.assertTrue(os.path.isfile(target))
            self.assertEqual(
                compute_md5(target), MD5SAMPLES["samplesIDs.HGDP.txt"]
            )


@unittest.skipIf(SKIP_NET, NET_REASON)
class TestHfFetchGenomePlacement(unittest.TestCase):
    """A small genome FASTA lands directly at the target .fa path (no untar)."""

    def test_small_genome_fasta_lands_at_target(self):
        with tempfile.TemporaryDirectory() as td:
            genome_dir = os.path.join(td, "Genomes", "hg38_chrM")
            target = os.path.join(genome_dir, "chrM.fa")
            out = hf_fetch(SMALL_GENOME_REPO_FILE, target)
            self.assertEqual(out, target)
            self.assertTrue(os.path.isfile(target))
            # extracted FASTA, non-empty, starts with a FASTA header
            self.assertGreater(os.path.getsize(target), 0)
            with open(target) as fh:
                self.assertTrue(fh.readline().startswith(">"))


class TestHfFetchFallbackTrigger(unittest.TestCase):
    """When the repo does not exist, hf_fetch raises so callers can fall back."""

    def setUp(self):
        self._orig = os.environ.get("CRISPRME_HF_DATA_REPO")
        os.environ["CRISPRME_HF_DATA_REPO"] = "lucapinello/does-not-exist-xyz-000"
        importlib.reload(utils)

    def tearDown(self):
        if self._orig is None:
            os.environ.pop("CRISPRME_HF_DATA_REPO", None)
        else:
            os.environ["CRISPRME_HF_DATA_REPO"] = self._orig
        importlib.reload(utils)

    @unittest.skipIf(SKIP_NET, NET_REASON)
    def test_nonexistent_repo_raises(self):
        # env override is picked up as the module constant
        self.assertEqual(utils.HF_DATA_REPO, "lucapinello/does-not-exist-xyz-000")
        with tempfile.TemporaryDirectory() as td:
            with self.assertRaises(Exception):
                utils.hf_fetch(
                    "samplesIDs/hg38_1000G.samplesID.txt",
                    os.path.join(td, "x.txt"),
                )


@unittest.skipIf(SKIP_NET, NET_REASON)
class TestDownloadSamplesIdsEndToEnd(unittest.TestCase):
    """The real download_samples_ids_data() uses the HF path end-to-end."""

    def test_complete_test_samples_ids_via_hf(self):
        # complete_test imports sibling modules by bare name; ensure importable
        import complete_test as ct  # noqa: E402

        with tempfile.TemporaryDirectory() as td:
            cwd = os.getcwd()
            os.chdir(td)
            try:
                ct.download_samples_ids_data("1000G")
                landed = os.path.join(
                    td, "samplesIDs", "hg38_1000G.samplesID.txt"
                )
                self.assertTrue(os.path.isfile(landed))
                self.assertEqual(
                    compute_md5(landed), MD5SAMPLES["samplesIDs.1000G.txt"]
                )
            finally:
                os.chdir(cwd)


if __name__ == "__main__":
    unittest.main(verbosity=2)
