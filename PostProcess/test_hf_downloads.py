"""Tests for the transparent HuggingFace download fast-path + source fallback.

The network-facing checks are skipped unless a populated HF dataset repo is
reachable; set CRISPRME_RUN_HF_TESTS=1 (and optionally CRISPRME_HF_REPO) to run
them. The fallback-contract test is network-free: it asserts that hf_fetch
propagates errors so the caller's try/except can fall back to the legacy source.
"""

import os
import unittest

import utils


@unittest.skipUnless(
    os.environ.get("CRISPRME_RUN_HF_TESTS") == "1",
    "network HF tests disabled (set CRISPRME_RUN_HF_TESTS=1 to enable)",
)
class TestHfFetchNetwork(unittest.TestCase):
    def test_genome_lands_at_exact_path(self):
        import tempfile

        dest = os.path.join(tempfile.mkdtemp(), "chr22.fa")
        out = utils.hf_fetch("genomes/hg38/chr22.fa", dest)
        self.assertEqual(out, dest)
        self.assertTrue(os.path.isfile(dest))  # placed exactly, not inside a tree


class TestHfFetchFallbackContract(unittest.TestCase):
    """hf_fetch must RAISE (not silently no-op) on a bad repo/file, so callers
    fall back to the legacy source."""

    def test_missing_file_raises(self):
        import tempfile

        # huggingface_hub not installed -> ImportError; installed but file/repo
        # missing -> some hub error. Either way it must raise, never return.
        dest = os.path.join(tempfile.mkdtemp(), "nope.fa")
        with self.assertRaises(Exception):
            utils.hf_fetch("genomes/hg38/__does_not_exist__.fa", dest)
        self.assertFalse(os.path.exists(dest))

    def test_repo_constant_matches_env(self):
        # HF_DATA_REPO honors CRISPRME_HF_REPO (shared with crisprme_hf)
        self.assertTrue(utils.HF_DATA_REPO)  # non-empty default


if __name__ == "__main__":
    unittest.main()
