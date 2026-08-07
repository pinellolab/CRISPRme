"""Regression test for CRISPRme's post-search annotation step
(``PostProcess/annotation.py``).

Run with:

    <env>/bin/python3 -m unittest discover -s PostProcess -p 'test_annotation.py' -v

``annotate_offtargets()`` used to split each off-target row with the
default, whitespace-collapsing ``str.split()`` instead of an explicit
``\\t``-split. Any row with a genuinely empty tab-separated column (e.g. a
blank AF value, reachable once CRISPRitz PR #36 lets AF-not-found return
``""`` instead of erroring) silently lost that column, shifting every later
column left by one and breaking the fixed-index readers downstream
(``add_risk_score.py``, ``resultIntegrator.py``) with an ``IndexError``.
"""

import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from PostProcess import annotation as an  # noqa: E402


# 22-column header used by the primary bestCFD/bestmmblg/bestCRISTA
# intermediates (submit_job_automated_new_multiple_vcfs.sh).
HEADER = (
    "#Bulge_type\tcrRNA\tDNA\tReference\tChromosome\tPosition\tCluster_Position\t"
    "Direction\tMismatches\tBulge_Size\tTotal\tPAM_gen\tVar_uniq\tSamples\t"
    "Annotation_Type\tReal_Guide\trsID\tAF\tSNP\t#Seq_in_cluster\tCFD\tCFD_ref"
)


def _row(af=""):
    """A realistic 22-field off-target row; AF (index 17) blank by default,
    matching the "no AF field anywhere in the VCF" / AF-bounds-overflow case.
    """
    fields = [
        "DNA", "CTAAC-AGTTGCTTTTATCACNNN", "CTccCAAGTTGCTgTgATCACAGG",
        "CTccCAAGTTGCTgggATCACAGG", "chr22", "10730842", "10730843", "+",
        "4", "1", "5", "n", "y", "HG04185,HG04171", "n",
        "CTAACAGTTGCTTTTATCACNNN", ".", af, "chr22_10730857_G_T", "0",
        "0.0", "0.0",
    ]
    assert len(fields) == 22
    return "\t".join(fields)


class AnnotateOfftargetsTest(unittest.TestCase):
    def _run(self, rows):
        with tempfile.TemporaryDirectory() as tmp:
            infile = os.path.join(tmp, "in.txt")
            outfile = os.path.join(tmp, "out.txt")
            with open(infile, "w") as fh:
                fh.write(HEADER + "\n")
                for row in rows:
                    fh.write(row + "\n")
            an.annotate_offtargets(infile, None, outfile)  # annotation=None: vuoto.txt path
            with open(outfile) as fh:
                lines = [line.rstrip("\n") for line in fh]
        return lines

    def test_blank_af_field_preserved(self):
        """A blank AF value must survive as its own empty column, not vanish."""
        out = self._run([_row(af="")])
        header, data = out[0], out[1]
        self.assertEqual(len(header.split("\t")), 22)
        fields = data.split("\t")
        self.assertEqual(len(fields), 22, f"expected 22 fields, got {len(fields)}: {fields}")
        self.assertEqual(fields[17], "", "AF column should be blank, not dropped/shifted")
        # everything after AF must still be at its original index
        self.assertEqual(fields[18], "chr22_10730857_G_T")
        self.assertEqual(fields[20], "0.0")
        self.assertEqual(fields[21], "0.0")

    def test_populated_af_field_unaffected(self):
        """Normal (non-blank) rows must be unaffected by the parsing fix."""
        out = self._run([_row(af="0.02")])
        fields = out[1].split("\t")
        self.assertEqual(len(fields), 22)
        self.assertEqual(fields[17], "0.02")
        self.assertEqual(fields[21], "0.0")

    def test_annotation_type_column_still_gets_rewritten(self):
        """annotate_target() writes into fields[14] (Annotation_Type); confirm
        the explicit-split fix didn't disturb that in-place rewrite."""
        out = self._run([_row(af="0.1")])
        fields = out[1].split("\t")
        self.assertEqual(fields[14], "n")  # annotation=None -> always "n"

    def test_mixed_blank_and_populated_rows(self):
        """A blank-AF row and a normal row in the same file must both come out
        with 22 fields (guards against any per-file, not just per-row, state)."""
        out = self._run([_row(af=""), _row(af="0.5")])
        for data in out[1:]:
            self.assertEqual(len(data.split("\t")), 22)


if __name__ == "__main__":
    unittest.main()
