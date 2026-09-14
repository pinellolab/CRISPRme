"""Regression test for GitHub issue #172: silent loss of ALL indel off-targets.

post_analisi_indel.sh subsets the per-chromosome indel targets out of the
combined CRISPRitz targets file. The indel search runs on the per-chromosome
*fake* genome (see pool_search_indels.py), so every target row's Chromosome
column is named ``fake${chrom}`` (e.g. ``fakechr22``), NOT ``${chrom}``.

A refactor changed the subset from
    awk -v fakechrom="$fakechrom" '$0 ~ fakechrom {print}'
to
    LC_ALL=C grep -F -w $chrom ...
i.e. it matched on ``$chrom`` ("chr22") instead of ``$fakechrom``
("fakechr22"). ``grep -F -w chr22`` can NEVER match ``fakechr22`` because the
left word boundary fails ("chr22" is preceded by the word char "e"), so the
subset came back empty, the entire indel post-analysis processed nothing, and
100% of indel off-targets were silently dropped (exit 0, clean logs).

These tests assert, against a tiny ``fakechr22``-columned fixture:
  * the FIXED command (grep -F -w "$fakechrom") returns the target rows, and
  * the BUGGY command (grep -F -w "$chrom") returns nothing,
  * ``-w`` still prevents ``fakechr2`` from matching ``fakechr22`` (guards the
    fix against over-matching), and
  * the shipped post_analisi_indel.sh actually uses ``$fakechrom`` (not the
    bare ``$chrom``) on both grep lines.

Run with:
    <env>/bin/python3 -m unittest discover -s PostProcess -p 'test_indel_chrom_subset.py' -v
"""
import os
import re
import subprocess
import tempfile
import unittest

PP = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(PP, "post_analisi_indel.sh")

# A minimal, well-formed CRISPRitz .targets.txt fragment. The Chromosome column
# (col 4) is "fakechr22" for the real indel rows, exactly as the indel search
# writes them. Each row has >=10 tab-separated columns so it survives the
# NF>=10 malformed-line guard. We include a decoy "fakechr2" row to prove -w
# does not let chr2 bleed into chr22.
HEADER = (
    "#Bulge_type\tcrRNA\tDNA\tChromosome\tPosition\tCluster_Position\t"
    "Direction\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples"
)
ROW_CHR22_A = (
    "X\tGACTGACTGACTGACTGACTNNN\tGACTGACTGACTGACTGACTAGG\tfakechr22\t"
    "1000\t1000\t+\t2\t0\t2\tNNN\ty\tNA"
)
ROW_CHR22_B = (
    "X\tGACTGACTGACTGACTGACTNNN\tGACTGACTGACTGACTGACTAGG\tfakechr22\t"
    "2000\t2000\t-\t1\t0\t1\tNNN\ty\tNA"
)
ROW_CHR2_DECOY = (
    "X\tGACTGACTGACTGACTGACTNNN\tGACTGACTGACTGACTGACTAGG\tfakechr2\t"
    "3000\t3000\t+\t0\t0\t0\tNNN\ty\tNA"
)


def _targets_file(tmpdir):
    path = os.path.join(tmpdir, "targets.txt")
    with open(path, "w") as fh:
        fh.write(HEADER + "\n")
        fh.write(ROW_CHR22_A + "\n")
        fh.write(ROW_CHR2_DECOY + "\n")
        fh.write(ROW_CHR22_B + "\n")
    return path


class IndelChromSubset(unittest.TestCase):
    def test_fixed_grep_returns_indel_targets(self):
        """grep -F -w "$fakechrom" recovers the fakechr22 rows (the fix)."""
        with tempfile.TemporaryDirectory() as td:
            targets = _targets_file(td)
            r = subprocess.run(
                ["bash", "-c",
                 'LC_ALL=C grep -F -w "$fakechrom" "$1"', "_", targets],
                capture_output=True, text=True,
                env={**os.environ, "fakechrom": "fakechr22"},
            )
            out = r.stdout.splitlines()
            self.assertEqual(len(out), 2, f"expected 2 fakechr22 rows, got: {out!r}")
            self.assertTrue(all("fakechr22" in ln for ln in out))
            self.assertFalse(any("fakechr2\t" in ln for ln in out),
                             "the fakechr2 decoy must NOT be in a chr22 subset")

    def test_buggy_grep_on_chrom_returns_nothing(self):
        """grep -F -w "$chrom" (the bug) yields an EMPTY subset."""
        with tempfile.TemporaryDirectory() as td:
            targets = _targets_file(td)
            r = subprocess.run(
                ["bash", "-c",
                 'LC_ALL=C grep -F -w "$chrom" "$1"', "_", targets],
                capture_output=True, text=True,
                env={**os.environ, "chrom": "chr22"},
            )
            self.assertEqual(r.stdout, "",
                             "buggy grep -w chr22 must not match fakechr22 -- "
                             "this empty subset is the silent-drop bug (#172)")

    def test_w_boundary_does_not_overmatch_chr2_into_chr22(self):
        """grep -F -w "$fakechrom"=fakechr2 must NOT catch fakechr22 rows."""
        with tempfile.TemporaryDirectory() as td:
            targets = _targets_file(td)
            r = subprocess.run(
                ["bash", "-c",
                 'LC_ALL=C grep -F -w "$fakechrom" "$1"', "_", targets],
                capture_output=True, text=True,
                env={**os.environ, "fakechrom": "fakechr2"},
            )
            out = r.stdout.splitlines()
            self.assertEqual(len(out), 1, f"expected only the fakechr2 decoy, got: {out!r}")
            self.assertIn("fakechr2\t", out[0])
            self.assertFalse(any("fakechr22" in ln for ln in out))

    def test_shipped_script_matches_on_fakechrom(self):
        """post_analisi_indel.sh must grep on $fakechrom, not bare $chrom."""
        with open(SCRIPT) as fh:
            text = fh.read()
        grep_lines = [ln.strip() for ln in text.splitlines()
                      if re.search(r"grep\b.*targets_tsv_(ref|alt)\b", ln)]
        self.assertEqual(len(grep_lines), 2,
                         f"expected 2 grep subset lines, found: {grep_lines!r}")
        for ln in grep_lines:
            self.assertIn('"$fakechrom"', ln,
                          f"subset grep must match on $fakechrom: {ln!r}")
            self.assertNotRegex(
                ln, r'grep\b[^\n]*-w\s+"?\$chrom"?\s',
                f"subset grep must NOT match on bare $chrom (bug #172): {ln!r}")


if __name__ == "__main__":
    unittest.main()
