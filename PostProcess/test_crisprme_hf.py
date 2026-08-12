"""Unit tests for the pure (network-free) logic of crisprme_hf.

These exercise repo/token resolution, gzip decompression, the index tar
round-trip, and argument validation without touching HuggingFace. The
network-facing calls (snapshot_download / upload_file) are covered separately by
an end-to-end test against a populated repo once one exists.
"""

import contextlib
import gzip
import importlib.util
import os
import tarfile
import unittest
from unittest import mock

import crisprme_hf as hf

# The network-free CI deliberately does NOT install huggingface_hub (see
# unit-tests.yml). Patching "huggingface_hub.get_token" is only possible when the
# module is importable; when it is absent, resolve_token already yields None via
# its own ImportError fallback, so no patch is needed.
_HAS_HF = importlib.util.find_spec("huggingface_hub") is not None


@contextlib.contextmanager
def _no_cached_hf_token():
    """Neutralize a cached ``huggingface-cli login`` token for the no-token path.

    A no-op when huggingface_hub is absent (resolve_token returns None on its own);
    patches get_token to None when it is present so a logged-in dev machine can't
    leak a real token into these tests.
    """
    if _HAS_HF:
        with mock.patch("huggingface_hub.get_token", return_value=None):
            yield
    else:
        yield


class TestResolveRepo(unittest.TestCase):
    def setUp(self):
        self._saved = {k: os.environ.get(k) for k in ("CRISPRME_HF_REPO", "HF_TOKEN",
                                                       "HUGGING_FACE_HUB_TOKEN")}
        for k in self._saved:
            os.environ.pop(k, None)

    def tearDown(self):
        for k, v in self._saved.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v

    def test_cli_wins(self):
        os.environ["CRISPRME_HF_REPO"] = "env/repo"
        self.assertEqual(hf.resolve_repo("cli/repo"), "cli/repo")

    def test_env_over_default(self):
        os.environ["CRISPRME_HF_REPO"] = "env/repo"
        self.assertEqual(hf.resolve_repo(None), "env/repo")

    def test_default(self):
        self.assertEqual(hf.resolve_repo(None), hf.DEFAULT_HF_REPO)

    def test_token_precedence(self):
        os.environ["HF_TOKEN"] = "tok_env"
        self.assertEqual(hf.resolve_token("tok_cli"), "tok_cli")
        self.assertEqual(hf.resolve_token(None), "tok_env")
        os.environ.pop("HF_TOKEN")
        os.environ["HUGGING_FACE_HUB_TOKEN"] = "tok_alt"
        self.assertEqual(hf.resolve_token(None), "tok_alt")

    def test_token_none(self):
        # neutralize a cached `huggingface-cli login` token so the no-token path
        # is exercised even on a machine that is logged into HuggingFace
        with _no_cached_hf_token():
            self.assertIsNone(hf.resolve_token(None))


class TestDecompress(unittest.TestCase):
    def test_gunzip(self):
        import tempfile

        d = tempfile.mkdtemp()
        raw = os.path.join(d, "chr1.fa")
        with open(raw, "w") as f:
            f.write(">chr1\nACGT\n")
        gzp = raw + ".gz"
        with open(raw, "rb") as src, gzip.open(gzp, "wb") as dst:
            dst.write(src.read())
        os.remove(raw)
        out = hf.decompress_gz(gzp)
        self.assertEqual(out, raw)
        self.assertTrue(os.path.exists(raw))
        self.assertFalse(os.path.exists(gzp))  # source removed by default
        with open(raw) as f:
            self.assertIn("ACGT", f.read())

    def test_noop_on_plain(self):
        self.assertEqual(hf.decompress_gz("/x/y/chr1.fa"), "/x/y/chr1.fa")


class TestComponentValidation(unittest.TestCase):
    def test_unknown_component(self):
        with self.assertRaises(ValueError):
            hf.download_component("nonsense", "/tmp")

    def test_vcf_requires_dataset(self):
        with self.assertRaises(ValueError):
            hf.download_component("vcf", "/tmp", dataset=None)

    def test_index_requires_name(self):
        with self.assertRaises(ValueError):
            hf.download_component("index", "/tmp", index_name=None)

    def test_component_maps_are_consistent(self):
        self.assertEqual(
            set(hf._COMPONENT_PREFIXES), set(hf._COMPONENT_LOCALDIR)
        )


class TestPublishValidation(unittest.TestCase):
    def test_missing_dir(self):
        with self.assertRaises(ValueError):
            hf.publish_index("/no/such/index/dir", token="x")

    def test_missing_token(self):
        import tempfile

        d = tempfile.mkdtemp()
        idx = os.path.join(d, "NGG_2_hg38")
        os.makedirs(idx)
        open(os.path.join(idx, "chr1.bin"), "wb").close()
        # no token in env, none passed, and no cached login -> must raise before
        # any network call (patch get_token so a logged-in machine can't cause a
        # real upload during the test)
        for k in ("HF_TOKEN", "HUGGING_FACE_HUB_TOKEN"):
            os.environ.pop(k, None)
        with _no_cached_hf_token():
            with self.assertRaises(ValueError):
                hf.publish_index(idx, token=None)


class TestIndexTarRoundTrip(unittest.TestCase):
    """The tar layout must unpack to genome_library/<name>/ (arcname=name)."""

    def test_arcname_layout(self):
        import tempfile

        d = tempfile.mkdtemp()
        idx = os.path.join(d, "genome_library", "NGG_2_hg38")
        os.makedirs(idx)
        with open(os.path.join(idx, "chr1.bin"), "wb") as f:
            f.write(b"\x00\x01")
        tarball = os.path.join(d, "NGG_2_hg38.tar.gz")
        with tarfile.open(tarball, "w:gz") as tf:
            tf.add(idx, arcname="NGG_2_hg38")
        out = os.path.join(d, "unpack")
        os.makedirs(out)
        with tarfile.open(tarball) as tf:
            names = tf.getnames()
            tf.extractall(out)
        self.assertIn("NGG_2_hg38", names[0])
        self.assertTrue(os.path.exists(os.path.join(out, "NGG_2_hg38", "chr1.bin")))


if __name__ == "__main__":
    unittest.main()
