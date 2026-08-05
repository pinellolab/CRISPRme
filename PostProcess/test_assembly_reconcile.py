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


_HG38_ROW_COLS = [
    "off_target_id", "hg38_chr", "hg38_start", "hg38_end",
    "Strand_(fewest_mm+b)", "CFD_score_(fewest_mm+b)",
]


def _hg38_row(off_target_id, chrom, start, strand="+", cfd=0.5, end=None):
    return {
        "off_target_id": off_target_id,
        "hg38_chr": chrom,
        "hg38_start": start,
        "hg38_end": end if end is not None else start + 1,
        "Strand_(fewest_mm+b)": strand,
        "CFD_score_(fewest_mm+b)": cfd,
    }


class TestCombineAcrossHaplotypes(unittest.TestCase):
    def _combine(self, a_rows, b_rows, merge_bp=3):
        # even an empty haplotype's real predictions frame still has the
        # full column schema (it's the result of a merge/filter, not a bare
        # empty list) -- pd.DataFrame([]) would have zero columns, which
        # isn't representative of real usage and hides schema-dependent bugs
        preds = {
            "a": pd.DataFrame(a_rows, columns=_HG38_ROW_COLS),
            "b": pd.DataFrame(b_rows, columns=_HG38_ROW_COLS),
        }
        return ar._combine_across_haplotypes(preds, ["a", "b"], merge_bp)

    def test_exact_match_is_both(self):
        out = self._combine(
            [_hg38_row("a1", "chr1", 1000)],
            [_hg38_row("b1", "chr1", 1000)],
        )
        self.assertEqual(len(out), 1)
        self.assertEqual(out["origin"].iloc[0], "both")

    def test_within_tolerance_drift_is_both(self):
        # the original bug: 2bp bulge-registration drift used to be reported
        # as two separate "_only" hits instead of one "both" hit
        out = self._combine(
            [_hg38_row("a1", "chr1", 1000)],
            [_hg38_row("b1", "chr1", 1002)],
            merge_bp=3,
        )
        self.assertEqual(len(out), 1)
        self.assertEqual(out["origin"].iloc[0], "both")

    def test_beyond_tolerance_stays_separate(self):
        out = self._combine(
            [_hg38_row("a1", "chr1", 1000)],
            [_hg38_row("b1", "chr1", 1010)],
            merge_bp=3,
        )
        self.assertEqual(len(out), 2)
        self.assertCountEqual(out["origin"].tolist(), ["a_only", "b_only"])

    def test_no_bridging_across_unrelated_distant_sites(self):
        # regression test for the chaining/bridging risk: sorted positions
        # 6000 (a1) -> 6002 (b1) -> 6005 (a2) have consecutive gaps of 2 and
        # 3, both <= merge_bp=3 -- a transitive-chaining approach (like
        # cluster_collapse's, if wrongly reused across haplotypes) would
        # bridge all three into ONE cluster and silently discard 2 of the 3
        # real distinct sites. a1 and a2 are 5bp apart directly (> merge_bp),
        # so they're genuinely separate sites, correctly kept apart within
        # their own haplotype's collapse already. Pairwise-only matching
        # must not silently lose a2 -- it has no direct exclusive claim on
        # b1 (a1 is closer, gap 2 vs 3) but must still surface as its own row.
        out = self._combine(
            [
                _hg38_row("a1", "chr1", 6000, cfd=0.1),
                _hg38_row("a2", "chr1", 6005, cfd=0.2),
            ],
            [_hg38_row("b1", "chr1", 6002, cfd=0.9)],
            merge_bp=3,
        )
        self.assertEqual(len(out), 2)  # never 1 -- that would mean a2 got silently merged away
        self.assertIn("a2", out["off_target_id_a"].dropna().tolist())
        both_rows = out[out["origin"] == "both"]
        self.assertEqual(both_rows["off_target_id_a"].iloc[0], "a1")

    def test_shared_nearest_neighbor_ties_break_to_closest_only(self):
        # two "a" sites each within tolerance of the same "b" site (possible
        # since the two "a" sites are >merge_bp apart from each other, so
        # cluster_collapse wouldn't have merged them within haplotype a) --
        # only the closer "a" site should get credited as "both"
        out = self._combine(
            [
                _hg38_row("a1", "chr1", 1000),
                _hg38_row("a2", "chr1", 1005),
            ],
            [_hg38_row("b1", "chr1", 1002)],
            merge_bp=3,
        )
        self.assertEqual(len(out), 2)
        both_rows = out[out["origin"] == "both"]
        self.assertEqual(len(both_rows), 1)
        self.assertEqual(both_rows["off_target_id_a"].iloc[0], "a1")
        a_only_rows = out[out["origin"] == "a_only"]
        self.assertEqual(a_only_rows["off_target_id_a"].iloc[0], "a2")

    def test_different_strand_never_matches(self):
        out = self._combine(
            [_hg38_row("a1", "chr1", 1000, strand="+")],
            [_hg38_row("b1", "chr1", 1000, strand="-")],
            merge_bp=3,
        )
        self.assertEqual(len(out), 2)
        self.assertCountEqual(out["origin"].tolist(), ["a_only", "b_only"])

    def test_different_chrom_never_matches(self):
        out = self._combine(
            [_hg38_row("a1", "chr1", 1000)],
            [_hg38_row("b1", "chr2", 1000)],
            merge_bp=3,
        )
        self.assertEqual(len(out), 2)
        self.assertCountEqual(out["origin"].tolist(), ["a_only", "b_only"])

    def test_empty_haplotype_b_all_a_only(self):
        out = self._combine(
            [_hg38_row("a1", "chr1", 1000)],
            [],
            merge_bp=3,
        )
        self.assertEqual(len(out), 1)
        self.assertEqual(out["origin"].iloc[0], "a_only")

    def test_multiple_strand_groups_does_not_raise_unsorted_keys(self):
        # regression test: merge_asof requires its "on" column (hg38_start)
        # to be sorted GLOBALLY across the whole frame, not just within each
        # `by` group. Sorting by (chrom, strand, start) -- correct for
        # cluster_collapse's chained clustering -- only guarantees hg38_start
        # increases *within* each (chrom, strand) group; here "+" sorts
        # before "-" but a1's position (5000) is higher than a2's (1000), so
        # the global sequence [5000, 1000] is not monotonic and used to raise
        # "left keys must be sorted". Any real result set with hits on both
        # strands (i.e. virtually all of them) hits this.
        out = self._combine(
            [
                _hg38_row("a1", "chr1", 5000, strand="+"),
                _hg38_row("a2", "chr1", 1000, strand="-"),
            ],
            [
                _hg38_row("b1", "chr1", 5001, strand="+"),
                _hg38_row("b2", "chr1", 1001, strand="-"),
            ],
            merge_bp=3,
        )
        self.assertEqual(len(out), 2)
        both_rows = out[out["origin"] == "both"]
        self.assertEqual(len(both_rows), 2)
        self.assertCountEqual(
            both_rows["off_target_id_a"].tolist(), ["a1", "a2"]
        )


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


class TestHaplotypeSearchComplete(unittest.TestCase):
    def test_false_when_directory_missing(self):
        self.assertFalse(ar.haplotype_search_complete("/no/such/dir"))

    def test_false_when_no_integrated_results(self):
        with tempfile.TemporaryDirectory() as d:
            open(os.path.join(d, ar.LOG_ERROR_NO_CHECK_FILENAME), "w").close()
            self.assertFalse(ar.haplotype_search_complete(d))

    def test_false_when_integrated_results_exist_but_no_completion_rename(self):
        # a run that crashed after writing results but before the pipeline's
        # final log_error.txt -> log_error_no_check.txt rename must not be
        # mistaken for a genuinely finished run
        with tempfile.TemporaryDirectory() as d:
            open(os.path.join(d, "guideX_PAM_genome_mm4_bMax2_integrated_results.tsv"), "w").close()
            self.assertFalse(ar.haplotype_search_complete(d))

    def test_false_when_alt_alignments_file_missing(self):
        # regression test: the integrated-results file and the alt-alignments
        # file are written by two separate (if normally back-to-back) steps
        # in the real pipeline -- a narrow interruption between them could
        # leave one without the other. load_crisprme_predictions reads the
        # alt-alignments file unconditionally, so treating this as
        # "complete" would defer straight to an uncaught FileNotFoundError
        # during reconciliation instead of a clean retry.
        with tempfile.TemporaryDirectory() as d:
            open(os.path.join(d, "guideX_PAM_genome_mm4_bMax2_integrated_results.tsv"), "w").close()
            open(os.path.join(d, ar.LOG_ERROR_NO_CHECK_FILENAME), "w").close()
            self.assertFalse(ar.haplotype_search_complete(d))

    def test_true_when_all_three_present(self):
        with tempfile.TemporaryDirectory() as d:
            open(os.path.join(d, "guideX_PAM_genome_mm4_bMax2_integrated_results.tsv"), "w").close()
            open(os.path.join(d, "guideX_PAM_genome_mm4_bMax2_all_results_with_alternative_alignments.tsv"), "w").close()
            open(os.path.join(d, ar.LOG_ERROR_NO_CHECK_FILENAME), "w").close()
            self.assertTrue(ar.haplotype_search_complete(d))


class TestCleanIncompleteHaplotypeOutput(unittest.TestCase):
    def test_removes_stale_incomplete_directory(self):
        with tempfile.TemporaryDirectory() as tmp:
            stale = os.path.join(tmp, "paternal_output")
            os.makedirs(stale)
            open(os.path.join(stale, "Params.txt"), "w").close()
            open(os.path.join(stale, "log_verbose.txt"), "w").close()

            ar.clean_incomplete_haplotype_output(stale)

            self.assertFalse(os.path.isdir(stale))

    def test_missing_directory_is_a_no_op(self):
        # a fresh, first-time run has nothing to clean up -- must not raise
        with tempfile.TemporaryDirectory() as tmp:
            never_existed = os.path.join(tmp, "does_not_exist")
            ar.clean_incomplete_haplotype_output(never_existed)  # should not raise
            self.assertFalse(os.path.isdir(never_existed))

    def test_does_not_touch_sibling_directories(self):
        # regression guard: cleaning up one haplotype's directory must never
        # affect a sibling directory (e.g. the other haplotype's own,
        # already-complete output), since assembly_search() calls this with
        # one haplotype's specific path at a time
        with tempfile.TemporaryDirectory() as tmp:
            stale = os.path.join(tmp, "paternal_output")
            sibling = os.path.join(tmp, "maternal_output")
            os.makedirs(stale)
            os.makedirs(sibling)
            open(os.path.join(sibling, "guideX_integrated_results.tsv"), "w").close()

            ar.clean_incomplete_haplotype_output(stale)

            self.assertFalse(os.path.isdir(stale))
            self.assertTrue(os.path.isdir(sibling))
            self.assertTrue(os.path.isfile(os.path.join(sibling, "guideX_integrated_results.tsv")))

    def test_custom_reason_used_in_log_message(self):
        with tempfile.TemporaryDirectory() as tmp:
            stale = os.path.join(tmp, "paternal_output")
            os.makedirs(stale)
            with patch("builtins.print") as mock_print:
                ar.clean_incomplete_haplotype_output(
                    stale, reason="results built with different search parameters"
                )
            mock_print.assert_called_once()
            self.assertIn("results built with different search parameters", mock_print.call_args[0][0])
            self.assertFalse(os.path.isdir(stale))


def _write_command_line_file(results_dir, genome, guide, pam, mm, bDNA, bRNA, merge_bp):
    os.makedirs(results_dir, exist_ok=True)
    cmd = (
        f"/usr/bin/python3 /path/crisprme.py complete-search --genome {genome} "
        f"--guide {guide} --pam {pam} --mm {mm} --bDNA {bDNA} --bRNA {bRNA} "
        f"--merge {merge_bp} --output some_name --thread 4"
    )
    with open(os.path.join(results_dir, ar.COMMAND_LINE_FILENAME), "w") as f:
        f.write(f"input_command\t{cmd}\n")


class TestHaplotypeParamsMatch(unittest.TestCase):
    def _params(self, **overrides):
        base = dict(
            genome_dir="/genomes/g1", guide_file="/guides/g.txt",
            pam_file="/pams/p.txt", mm=4, bDNA=1, bRNA=1, merge_bp=3,
        )
        base.update(overrides)
        return base

    def test_matches_when_all_params_identical(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = self._params()
            _write_command_line_file(
                tmp, p["genome_dir"], p["guide_file"], p["pam_file"],
                p["mm"], p["bDNA"], p["bRNA"], p["merge_bp"],
            )
            self.assertTrue(ar.haplotype_params_match(tmp, **p))

    def test_no_match_when_genome_changed(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = self._params()
            _write_command_line_file(
                tmp, p["genome_dir"], p["guide_file"], p["pam_file"],
                p["mm"], p["bDNA"], p["bRNA"], p["merge_bp"],
            )
            changed = self._params(genome_dir="/genomes/g2_different")
            self.assertFalse(ar.haplotype_params_match(tmp, **changed))

    def test_no_match_when_guide_changed(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = self._params()
            _write_command_line_file(
                tmp, p["genome_dir"], p["guide_file"], p["pam_file"],
                p["mm"], p["bDNA"], p["bRNA"], p["merge_bp"],
            )
            changed = self._params(guide_file="/guides/different_guide.txt")
            self.assertFalse(ar.haplotype_params_match(tmp, **changed))

    def test_no_match_when_mm_changed(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = self._params()
            _write_command_line_file(
                tmp, p["genome_dir"], p["guide_file"], p["pam_file"],
                p["mm"], p["bDNA"], p["bRNA"], p["merge_bp"],
            )
            changed = self._params(mm=6)
            self.assertFalse(ar.haplotype_params_match(tmp, **changed))

    def test_no_match_when_command_line_file_missing(self):
        # a genuinely incomplete/crashed run that never got far enough to
        # write .command_line.txt -- must not be treated as reusable
        with tempfile.TemporaryDirectory() as tmp:
            self.assertFalse(ar.haplotype_params_match(tmp, **self._params()))

    def test_no_match_when_command_line_file_malformed(self):
        with tempfile.TemporaryDirectory() as tmp:
            with open(os.path.join(tmp, ar.COMMAND_LINE_FILENAME), "w") as f:
                f.write("not the expected format at all\n")
            self.assertFalse(ar.haplotype_params_match(tmp, **self._params()))


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

    def test_returns_dropped_ids_for_unmapped_chrom(self):
        # the dropped-chrom row must be surfaced, not just silently excluded
        # from the BED -- reconcile_haplotypes folds these into non_mappable
        preds = pd.DataFrame({
            "Chromosome": ["chr1", "chrUn_unknown", "chrUn_other"],
            "Start_coordinate_(fewest_mm+b)": [1000, 2000, 3000],
            "off_target_id": ["0", "1", "2"],
        })
        ucsc_to_genbank = {"chr1": "CM000001.1"}
        with tempfile.TemporaryDirectory() as d:
            bed_path = os.path.join(d, "out.bed")
            _, dropped_ids = ar.build_offtarget_bed(preds, ucsc_to_genbank, bed_path)
            self.assertEqual(dropped_ids, {"1", "2"})

    def test_no_dropped_ids_when_all_chroms_known(self):
        preds = pd.DataFrame({
            "Chromosome": ["chr1"],
            "Start_coordinate_(fewest_mm+b)": [1000],
            "off_target_id": ["0"],
        })
        ucsc_to_genbank = {"chr1": "CM000001.1"}
        with tempfile.TemporaryDirectory() as d:
            bed_path = os.path.join(d, "out.bed")
            _, dropped_ids = ar.build_offtarget_bed(preds, ucsc_to_genbank, bed_path)
            self.assertEqual(dropped_ids, set())


class TestCheckChromAliasCoverage(unittest.TestCase):
    def test_passes_when_all_chroms_known(self):
        preds = pd.DataFrame({"Chromosome": ["chr1", "chr1", "chr2"]})
        ar.check_chrom_alias_coverage(preds, {"chr1": "x", "chr2": "y"}, "paternal")

    def test_passes_when_a_minority_of_chroms_are_unknown(self):
        # a handful of unmapped decoy/unplaced contigs is expected biology,
        # not a naming mismatch -- must not raise
        preds = pd.DataFrame({"Chromosome": ["chr1"] * 9 + ["chrUn_decoy"]})
        ar.check_chrom_alias_coverage(preds, {"chr1": "x"}, "paternal")

    def test_raises_when_majority_of_chroms_are_unknown(self):
        # a systematic mismatch (wrong alias file, wrong naming convention)
        # must be a clear, immediate error, not a silent near-total drop
        preds = pd.DataFrame({"Chromosome": ["chrom1", "chrom1", "chrom2"]})
        with self.assertRaises(ValueError) as ctx:
            ar.check_chrom_alias_coverage(preds, {"chr1": "x"}, "paternal")
        self.assertIn("paternal", str(ctx.exception))
        self.assertIn("chromAlias", str(ctx.exception))

    def test_empty_predictions_do_not_raise(self):
        preds = pd.DataFrame({"Chromosome": []})
        ar.check_chrom_alias_coverage(preds, {"chr1": "x"}, "paternal")


class TestCheckLiftoverFailureRate(unittest.TestCase):
    def test_passes_when_most_predictions_lift_successfully(self):
        ar.check_liftover_failure_rate(attempted_count=100, unlifted_count=10, name="paternal")

    def test_raises_when_majority_rejected(self):
        # a wrong/swapped chain file rejects almost everything -- must be a
        # clear, immediate error, not a silent near-total "haplotype-private"
        # result
        with self.assertRaises(ValueError) as ctx:
            ar.check_liftover_failure_rate(attempted_count=100, unlifted_count=60, name="paternal")
        self.assertIn("paternal", str(ctx.exception))
        self.assertIn("chain file", str(ctx.exception))

    def test_zero_attempted_does_not_raise(self):
        # everything already got dropped for missing chromAlias coverage
        # (that's check_chrom_alias_coverage's job to catch) -- nothing left
        # to divide by here
        ar.check_liftover_failure_rate(attempted_count=0, unlifted_count=0, name="paternal")

    def test_exactly_at_threshold_does_not_raise(self):
        # ratio > threshold raises, ratio == threshold does not
        ar.check_liftover_failure_rate(attempted_count=100, unlifted_count=50, name="paternal")


class TestCheckLiftoverAvailable(unittest.TestCase):
    def test_raises_when_not_on_path(self):
        with patch.object(ar.shutil, "which", return_value=None):
            with self.assertRaises(RuntimeError) as ctx:
                ar.check_liftover_available()
        self.assertIn("liftOver", str(ctx.exception))

    def test_passes_when_on_path(self):
        with patch.object(ar.shutil, "which", return_value="/usr/bin/liftOver"):
            ar.check_liftover_available()  # must not raise


class TestRunLiftover(unittest.TestCase):
    def test_calls_liftover_directly_without_mamba_env_wrapper(self):
        # regression guard: this used to hardcode `mamba run -n liftover_env`,
        # an undocumented dependency not part of crisprme's own environment
        with patch.object(ar.subprocess, "run") as mock_run:
            mock_run.return_value.returncode = 0
            ar.run_liftover("in.bed", "chain.file", "mapped.bed", "unmapped.bed")
        cmd = mock_run.call_args[0][0]
        self.assertNotIn("mamba", cmd)
        self.assertNotIn("liftover_env", cmd)
        self.assertTrue(cmd.startswith("liftOver "))


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
        # 3 loci per haplotype, not 1: check_liftover_failure_rate exists
        # specifically to catch a suspiciously high liftOver rejection rate,
        # so a fixture testing "one genuine rejection is handled correctly"
        # needs enough attempted predictions that one rejection is a normal,
        # low ratio -- not the exact failure mode that check guards against.
        with tempfile.TemporaryDirectory() as tmpdir:
            paternal_dir = self._make_haplotype_results(
                tmpdir, "paternal",
                [
                    self._pred_row("chr1", 1000, 0.9),
                    self._pred_row("chr1", 5000, 0.5),
                    self._pred_row("chr1", 9000, 0.3),
                ],
            )
            maternal_dir = self._make_haplotype_results(
                tmpdir, "maternal",
                [
                    self._pred_row("chr1", 1000, 0.9),
                    self._pred_row("chr1", 5000, 0.5),
                    self._pred_row("chr1", 9000, 0.3),
                ],
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

            # both haplotypes' loci at 5000/9000 lift identically (-> "both");
            # only maternal's locus at 1000 fails to lift at all, which must
            # still be counted as non-mappable, not silently dropped.
            lift_map = {5000: 60000, 9000: 70000}

            def fake_run_liftover(bed_path, chain_file, mapped_path, unmapped_path, env="liftover_env"):
                bed = pd.read_csv(bed_path, sep="\t", header=None,
                                   names=["chrom", "start", "end", "off_target_id"])
                with open(mapped_path, "w") as mf, open(unmapped_path, "w") as uf:
                    for _, r in bed.iterrows():
                        if r["end"] == 1000 and "maternal" in bed_path:
                            uf.write("#Deleted in new\n")
                            uf.write(f"{r['chrom']}\t{r['start']}\t{r['end']}\t{r['off_target_id']}\n")
                        elif r["end"] == 1000:
                            mf.write(f"chr1\t49999\t50000\t{r['off_target_id']}\n")
                        else:
                            hg38_end = lift_map[r["end"]]
                            mf.write(f"chr1\t{hg38_end - 1}\t{hg38_end}\t{r['off_target_id']}\n")
                return mapped_path, unmapped_path

            with patch.object(ar, "run_liftover", side_effect=fake_run_liftover):
                combined, summary = ar.reconcile_haplotypes(haplotypes, workdir=tmpdir, merge_bp=3)

            self.assertEqual(summary["both"], 2)
            self.assertEqual(summary["paternal_only"], 1)
            self.assertEqual(summary["maternal_only"], 0)
            self.assertEqual(summary["paternal_non_mappable"], 0)
            self.assertEqual(summary["maternal_non_mappable"], 1)

    def test_chrom_alias_missing_contig_folded_into_non_mappable(self):
        # a prediction whose chromosome has no chromAlias mapping never even
        # reaches liftOver (build_offtarget_bed drops it beforehand) -- must
        # still be counted as non_mappable, not silently vanish from summary
        with tempfile.TemporaryDirectory() as tmpdir:
            paternal_dir = self._make_haplotype_results(
                tmpdir, "paternal", [self._pred_row("chr1", 1000, 0.9)],
            )
            maternal_dir = self._make_haplotype_results(
                tmpdir, "maternal",
                [self._pred_row("chr1", 1000, 0.9), self._pred_row("chrUn_weird", 2000, 0.8)],
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

            def fake_run_liftover(bed_path, chain_file, mapped_path, unmapped_path):
                bed = pd.read_csv(bed_path, sep="\t", header=None,
                                   names=["chrom", "start", "end", "off_target_id"])
                with open(mapped_path, "w") as mf, open(unmapped_path, "w") as uf:
                    for _, r in bed.iterrows():
                        mf.write(f"chr1\t49999\t50000\t{r['off_target_id']}\n")
                return mapped_path, unmapped_path

            with patch.object(ar, "run_liftover", side_effect=fake_run_liftover):
                combined, summary = ar.reconcile_haplotypes(haplotypes, workdir=tmpdir, merge_bp=3)

            self.assertEqual(summary["both"], 1)
            self.assertEqual(summary["paternal_non_mappable"], 0)
            self.assertEqual(summary["maternal_non_mappable"], 1)

    def test_raises_on_systematic_chrom_alias_mismatch(self):
        # end-to-end: a haplotype whose predictions are on a totally
        # different naming convention than its chromAlias file should fail
        # clearly, not silently report near-total non-mappability
        with tempfile.TemporaryDirectory() as tmpdir:
            paternal_dir = self._make_haplotype_results(
                tmpdir, "paternal", [self._pred_row("chr1", 1000, 0.9)],
            )
            maternal_dir = self._make_haplotype_results(
                tmpdir, "maternal",
                [self._pred_row("chrom1", 1000, 0.9), self._pred_row("chrom1", 2000, 0.8)],
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
            def fake_run_liftover(bed_path, chain_file, mapped_path, unmapped_path):
                # only needs to get paternal through cleanly -- maternal's
                # coverage check should raise before its own liftOver call
                bed = pd.read_csv(bed_path, sep="\t", header=None,
                                   names=["chrom", "start", "end", "off_target_id"])
                with open(mapped_path, "w") as mf, open(unmapped_path, "w") as uf:
                    for _, r in bed.iterrows():
                        mf.write(f"chr1\t49999\t50000\t{r['off_target_id']}\n")
                return mapped_path, unmapped_path

            with patch.object(ar, "run_liftover", side_effect=fake_run_liftover):
                with self.assertRaises(ValueError) as ctx:
                    ar.reconcile_haplotypes(haplotypes, workdir=tmpdir, merge_bp=3)
            self.assertIn("maternal", str(ctx.exception))

            # failure must be recorded in log_error.txt, the same filename
            # complete-search/complete-test write to on failure
            log_error_path = os.path.join(tmpdir, ar.LOG_ERROR_FILENAME)
            with open(log_error_path) as f:
                error_log = f.read()
            self.assertIn("maternal", error_log)

    def _reconcile_simple_success(self, tmpdir):
        """Minimal always-succeeds fixture, for testing the log files
        themselves rather than the reconciliation categories."""
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

        def fake_run_liftover(bed_path, chain_file, mapped_path, unmapped_path):
            bed = pd.read_csv(bed_path, sep="\t", header=None,
                               names=["chrom", "start", "end", "off_target_id"])
            with open(mapped_path, "w") as mf, open(unmapped_path, "w") as uf:
                for _, r in bed.iterrows():
                    mf.write(f"chr1\t49999\t50000\t{r['off_target_id']}\n")
            return mapped_path, unmapped_path

        with patch.object(ar, "run_liftover", side_effect=fake_run_liftover):
            return ar.reconcile_haplotypes(haplotypes, workdir=tmpdir, merge_bp=3)

    def test_writes_log_verbose_and_empty_log_error_on_success(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            self._reconcile_simple_success(tmpdir)

            with open(os.path.join(tmpdir, ar.LOG_VERBOSE_FILENAME)) as f:
                verbose_log = f.read()
            self.assertIn("Reconciliation complete", verbose_log)
            self.assertIn("both: 1", verbose_log)

            with open(os.path.join(tmpdir, ar.LOG_ERROR_FILENAME)) as f:
                error_log = f.read()
            self.assertEqual(error_log, "")

    def test_log_files_truncate_fresh_on_each_call_not_appended(self):
        # matches complete_search()'s own open(..., "w") convention -- a
        # retry into the same workdir must not accumulate log content
        # across attempts
        with tempfile.TemporaryDirectory() as tmpdir:
            self._reconcile_simple_success(tmpdir)
            with open(os.path.join(tmpdir, ar.LOG_VERBOSE_FILENAME)) as f:
                first_call_log = f.read()
            self.assertNotEqual(first_call_log, "")

            self._reconcile_simple_success(tmpdir)
            with open(os.path.join(tmpdir, ar.LOG_VERBOSE_FILENAME)) as f:
                second_call_log = f.read()

            self.assertEqual(second_call_log, first_call_log)  # not doubled


if __name__ == "__main__":
    unittest.main()
