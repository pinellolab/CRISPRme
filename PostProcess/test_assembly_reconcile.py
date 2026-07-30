"""Unit test suite for CRISPRme's diploid-assembly reconciliation library
(``PostProcess/assembly_reconcile.py``).

Run with:

    <env>/bin/python3 -m unittest discover -s PostProcess -p 'test_assembly_reconcile.py' -v

or directly:

    <env>/bin/python3 PostProcess/test_assembly_reconcile.py -v

Each test builds its fixtures on the fly (in-memory DataFrames or files under
``tempfile.TemporaryDirectory()``) so the suite is self-contained and needs
neither real haplotype genomes nor the external ``liftOver`` binary/env --
``run_liftover``'s only external-tool call is stubbed out via
``unittest.mock.patch`` wherever a test needs to exercise the code around it
(``reconcile_haplotypes``). The two duplicate-row bug cases (exact-coordinate
and bulge-shifted near-duplicate) fixed by ``cluster_collapse()`` are the
same cases documented in ``assembly_search_plan.md`` and caught against real
HG01255 data; here they're reproduced with small synthetic fixtures instead.
"""

import os
import sys
import tempfile
import unittest
from unittest.mock import patch

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import assembly_reconcile as ar  # noqa: E402


class TestClusterCollapse(unittest.TestCase):
    def test_distinct_clusters_kept_separate(self):
        df = pd.DataFrame({
            "chrom": ["chr1", "chr1"],
            "strand": ["+", "+"],
            "pos": [100, 200],
            "score": [0.5, 0.9],
        })
        out = ar.cluster_collapse(df, "chrom", "strand", "pos", "score", merge_bp=3)
        self.assertEqual(len(out), 2)
        self.assertCountEqual(out["pos"].tolist(), [100, 200])

    def test_exact_coordinate_duplicates_collapse_to_best_score(self):
        # Mirrors the chr9:67041469 case: many rows, same locus, best CFD wins.
        df = pd.DataFrame({
            "chrom": ["chr9"] * 4,
            "strand": ["-"] * 4,
            "pos": [67041469] * 4,
            "score": [1.000, 0.197, 0.5, 0.8],
        })
        out = ar.cluster_collapse(df, "chrom", "strand", "pos", "score", merge_bp=3)
        self.assertEqual(len(out), 1)
        self.assertAlmostEqual(out["score"].iloc[0], 1.000)

    def test_bulge_shifted_near_duplicates_chain_into_one_cluster(self):
        # Mirrors the chr1:607952/607953 case: 1bp apart, same site, best score wins.
        df = pd.DataFrame({
            "chrom": ["chr1", "chr1"],
            "strand": ["+", "+"],
            "pos": [607952, 607953],
            "score": [0.034, 0.080],
        })
        out = ar.cluster_collapse(df, "chrom", "strand", "pos", "score", merge_bp=3)
        self.assertEqual(len(out), 1)
        self.assertEqual(out["pos"].iloc[0], 607953)

    def test_chained_gaps_merge_even_beyond_single_step_threshold(self):
        # 100 -> 102 -> 104 -> 106: each hop is 2bp (<= merge_bp=3), so all
        # four chain into one cluster even though 100->106 spans 6bp.
        df = pd.DataFrame({
            "chrom": ["chr1"] * 4,
            "strand": ["+"] * 4,
            "pos": [100, 102, 104, 106],
            "score": [0.1, 0.9, 0.2, 0.3],
        })
        out = ar.cluster_collapse(df, "chrom", "strand", "pos", "score", merge_bp=3)
        self.assertEqual(len(out), 1)
        self.assertEqual(out["pos"].iloc[0], 102)

    def test_gap_over_threshold_starts_new_cluster(self):
        df = pd.DataFrame({
            "chrom": ["chr1", "chr1"],
            "strand": ["+", "+"],
            "pos": [100, 105],
            "score": [0.5, 0.5],
        })
        out = ar.cluster_collapse(df, "chrom", "strand", "pos", "score", merge_bp=3)
        self.assertEqual(len(out), 2)

    def test_different_chrom_or_strand_never_merge(self):
        df = pd.DataFrame({
            "chrom": ["chr1", "chr2"],
            "strand": ["+", "+"],
            "pos": [100, 101],
            "score": [0.5, 0.5],
        })
        out = ar.cluster_collapse(df, "chrom", "strand", "pos", "score", merge_bp=3)
        self.assertEqual(len(out), 2)

    def test_nan_score_sorts_last_and_loses_ties(self):
        df = pd.DataFrame({
            "chrom": ["chr1", "chr1"],
            "strand": ["+", "+"],
            "pos": [100, 101],
            "score": [float("nan"), 0.1],
        })
        out = ar.cluster_collapse(df, "chrom", "strand", "pos", "score", merge_bp=3)
        self.assertEqual(len(out), 1)
        self.assertEqual(out["pos"].iloc[0], 101)


class TestFindResultsPrefix(unittest.TestCase):
    def test_finds_unique_prefix(self):
        with tempfile.TemporaryDirectory() as d:
            open(os.path.join(d, "guideX_PAM_genome_mm4_bMax2_integrated_results.tsv"), "w").close()
            prefix = ar.find_results_prefix(d)
            self.assertEqual(prefix, "guideX_PAM_genome_mm4_bMax2")

    def test_raises_when_no_match(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(FileNotFoundError):
                ar.find_results_prefix(d)

    def test_raises_when_multiple_matches(self):
        with tempfile.TemporaryDirectory() as d:
            open(os.path.join(d, "a_integrated_results.tsv"), "w").close()
            open(os.path.join(d, "b_integrated_results.tsv"), "w").close()
            with self.assertRaises(FileNotFoundError):
                ar.find_results_prefix(d)


class TestLoadChromAlias(unittest.TestCase):
    def test_parses_assembly_ucsc_genbank_columns(self):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "sample.chromAlias.txt")
            with open(path, "w") as f:
                f.write("# assembly\tucsc\tgenbank\n")
                f.write("chr1_h1\tchr1\tCM000001.1\n")
                f.write("chr2_h1\tchr2\tCM000002.1\n")
            assembly_to_ucsc, ucsc_to_genbank = ar.load_chrom_alias(path)
            self.assertEqual(assembly_to_ucsc["chr1_h1"], "chr1")
            self.assertEqual(ucsc_to_genbank["chr2"], "CM000002.1")


class TestBuildOfftargetBed(unittest.TestCase):
    def test_writes_bed_and_drops_unmapped_chrom(self):
        preds = pd.DataFrame({
            "Chromosome": ["chr1", "chrUn_unknown"],
            "Start_coordinate_(fewest_mm+b)": [1000, 2000],
            "off_target_id": ["0", "1"],
        })
        ucsc_to_genbank = {"chr1": "CM000001.1"}
        with tempfile.TemporaryDirectory() as d:
            bed_path = os.path.join(d, "out.bed")
            ar.build_offtarget_bed(preds, ucsc_to_genbank, bed_path)
            with open(bed_path) as f:
                lines = f.read().splitlines()
            self.assertEqual(len(lines), 1)
            chrom, start, end, oid = lines[0].split("\t")
            self.assertEqual(chrom, "CM000001.1")
            self.assertEqual(int(start), 999)
            self.assertEqual(int(end), 1000)
            self.assertEqual(oid, "0")


class TestLoadLiftedBed(unittest.TestCase):
    def test_parses_mapped_bed(self):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "mapped.bed")
            with open(path, "w") as f:
                f.write("chr1\t999\t1000\t0\n")
            df = ar.load_lifted_bed(path)
            self.assertEqual(df.loc[0, "hg38_chr"], "chr1")
            self.assertEqual(df.loc[0, "off_target_id"], "0")


class TestLoadUnliftedIds(unittest.TestCase):
    def test_skips_comment_and_blank_lines(self):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "unmapped.bed")
            with open(path, "w") as f:
                f.write("#Deleted in new\n")
                f.write("chrUn\t1\t2\t7\n")
                f.write("\n")
            ids = ar.load_unlifted_ids(path)
            self.assertEqual(ids, {"7"})


class TestReconcileHaplotypes(unittest.TestCase):
    def test_rejects_ploidy_other_than_two(self):
        with self.assertRaises(ValueError):
            ar.reconcile_haplotypes({"paternal": {}}, workdir="/tmp/unused")

    def _make_haplotype_results(self, tmpdir, name, rows):
        """Writes a minimal *_integrated_results.tsv (+ empty alt-alignments
        file) for one haplotype, using the real PRED_COLS schema."""
        results_dir = os.path.join(tmpdir, f"{name}_results")
        os.makedirs(results_dir, exist_ok=True)
        prefix = f"guide_PAM_{name}_mm4_bMax2"
        df = pd.DataFrame(rows, columns=ar.PRED_COLS)
        df.to_csv(
            os.path.join(results_dir, f"{prefix}_integrated_results.tsv"),
            sep="\t", index=False,
        )
        pd.DataFrame(columns=ar.PRED_COLS).to_csv(
            os.path.join(results_dir, f"{prefix}_all_results_with_alternative_alignments.tsv"),
            sep="\t", index=False,
        )
        return results_dir

    def _make_chrom_alias(self, tmpdir, name):
        path = os.path.join(tmpdir, f"{name}.chromAlias.txt")
        with open(path, "w") as f:
            f.write("# assembly\tucsc\tgenbank\n")
            f.write("chr1_asm\tchr1\tCM000001.1\n")
        return path

    def _pred_row(self, chrom, pos, cfd):
        return {
            "Spacer+PAM": "ACGT", "Chromosome": chrom,
            "Start_coordinate_(fewest_mm+b)": pos, "Strand_(fewest_mm+b)": "+",
            "Aligned_spacer+PAM_(fewest_mm+b)": "ACGT",
            "Aligned_protospacer+PAM_REF_(fewest_mm+b)": "ACGT",
            "Aligned_protospacer+PAM_ALT_(fewest_mm+b)": "ACGT",
            "Mismatches_(fewest_mm+b)": 0, "Bulges_(fewest_mm+b)": 0,
            "CFD_score_(fewest_mm+b)": cfd,
        }

    def test_categorizes_both_only_and_non_mappable(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            paternal_dir = self._make_haplotype_results(
                tmpdir, "paternal",
                [self._pred_row("chr1", 1000, 0.9), self._pred_row("chr1", 5000, 0.5)],
            )
            maternal_dir = self._make_haplotype_results(
                tmpdir, "maternal",
                [self._pred_row("chr1", 1000, 0.9), self._pred_row("chr1", 9000, 0.4)],
            )
            haplotypes = {
                "paternal": {
                    "chrom_alias_file": self._make_chrom_alias(tmpdir, "paternal"),
                    "chain_file": "unused.chain",
                    "results_dir": paternal_dir,
                },
                "maternal": {
                    "chrom_alias_file": self._make_chrom_alias(tmpdir, "maternal"),
                    "chain_file": "unused.chain",
                    "results_dir": maternal_dir,
                },
            }

            # Stub run_liftover: shared locus (1000) lifts identically for both
            # haplotypes (-> "both"); each haplotype's second locus lifts to a
            # haplotype-unique hg38 position (-> "<name>_only"); nothing is
            # left un-liftable in this fixture (covered separately below).
            def fake_run_liftover(bed_path, chain_file, mapped_path, unmapped_path, env="liftover_env"):
                bed = pd.read_csv(bed_path, sep="\t", header=None,
                                   names=["chrom", "start", "end", "off_target_id"])
                lift_map = {1000: 50000, 5000: 60000, 9000: 70000}
                with open(mapped_path, "w") as mf, open(unmapped_path, "w") as uf:
                    for _, r in bed.iterrows():
                        hg38_end = lift_map[r["end"]]
                        mf.write(f"chr1\t{hg38_end - 1}\t{hg38_end}\t{r['off_target_id']}\n")
                return mapped_path, unmapped_path

            with patch.object(ar, "run_liftover", side_effect=fake_run_liftover):
                combined, summary = ar.reconcile_haplotypes(haplotypes, workdir=tmpdir, merge_bp=3)

            self.assertEqual(summary["both"], 1)
            self.assertEqual(summary["paternal_only"], 1)
            self.assertEqual(summary["maternal_only"], 1)
            self.assertEqual(summary["paternal_non_mappable"], 0)
            self.assertEqual(summary["maternal_non_mappable"], 0)
            self.assertEqual(len(combined), 3)

    def test_unliftable_locus_counted_as_non_mappable_not_dropped(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            paternal_dir = self._make_haplotype_results(
                tmpdir, "paternal", [self._pred_row("chr1", 1000, 0.9)],
            )
            maternal_dir = self._make_haplotype_results(
                tmpdir, "maternal", [self._pred_row("chr1", 1000, 0.9)],
            )
            haplotypes = {
                "paternal": {
                    "chrom_alias_file": self._make_chrom_alias(tmpdir, "paternal"),
                    "chain_file": "unused.chain",
                    "results_dir": paternal_dir,
                },
                "maternal": {
                    "chrom_alias_file": self._make_chrom_alias(tmpdir, "maternal"),
                    "chain_file": "unused.chain",
                    "results_dir": maternal_dir,
                },
            }

            def fake_run_liftover(bed_path, chain_file, mapped_path, unmapped_path, env="liftover_env"):
                bed = pd.read_csv(bed_path, sep="\t", header=None,
                                   names=["chrom", "start", "end", "off_target_id"])
                with open(mapped_path, "w") as mf, open(unmapped_path, "w") as uf:
                    if "paternal" in bed_path:
                        for _, r in bed.iterrows():
                            mf.write(f"chr1\t49999\t50000\t{r['off_target_id']}\n")
                    else:
                        uf.write("#Deleted in new\n")
                        for _, r in bed.iterrows():
                            uf.write(f"{r['chrom']}\t{r['start']}\t{r['end']}\t{r['off_target_id']}\n")
                return mapped_path, unmapped_path

            with patch.object(ar, "run_liftover", side_effect=fake_run_liftover):
                combined, summary = ar.reconcile_haplotypes(haplotypes, workdir=tmpdir, merge_bp=3)

            self.assertEqual(summary["both"], 0)
            self.assertEqual(summary["paternal_only"], 1)
            self.assertEqual(summary["maternal_only"], 0)
            self.assertEqual(summary["paternal_non_mappable"], 0)
            self.assertEqual(summary["maternal_non_mappable"], 1)


if __name__ == "__main__":
    unittest.main()
