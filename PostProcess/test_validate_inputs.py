"""Unit test suite for CRISPRme's pre-flight input validator
(``PostProcess/validate_inputs.py``).

Run with:

    <env>/bin/python3 -m unittest discover -s PostProcess -p 'test_validate_inputs.py' -v

or directly:

    <env>/bin/python3 PostProcess/test_validate_inputs.py -v

Every check function in ``validate_inputs.py`` is a pure function over file
paths / already-parsed data, so each check gets its own ``TestCase`` with a
"clean" case and the specific failure mode(s) it exists to catch (see
``validate_inputs_plan.md`` at the repo root for the rationale behind each
check). Fixtures are built on the fly under ``tempfile.TemporaryDirectory()``
so the suite is fully self-contained and leaves nothing behind.
"""

import contextlib
import gzip
import io
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from PostProcess import validate_inputs as vi  # noqa: E402


# ---------------------------------------------------------------------------
# Fixture helpers
# ---------------------------------------------------------------------------

DEFAULT_VCF_HEADER = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1"


def write_text(path, text):
    with open(path, "w") as fh:
        fh.write(text)


def write_gzip_text(path, text):
    with gzip.open(path, "wt") as fh:
        fh.write(text)


def vcf_record(chrom="chr1", pos=100, vid=".", ref="A", alt="T", qual=".",
               filt="PASS", info="AF=0.1", fmt="GT", samples=("0/1",)):
    """Builds one tab-separated VCF data record."""
    return "\t".join([chrom, str(pos), vid, ref, alt, qual, filt, info, fmt, *samples])


def make_vcf_gz(path, records, header=DEFAULT_VCF_HEADER,
                 extra_header_lines=("##fileformat=VCFv4.2",)):
    """Writes a (plain-gzip, not true BGZF) VCF fixture -- sufficient since
    every check under test only ever calls `gzip.open`, matching how
    `enricher.py` itself opens VCFs."""
    lines = list(extra_header_lines) + [header] + list(records)
    write_gzip_text(path, "\n".join(lines) + "\n")


def make_fasta(path, chrom, sequence, line_width=60):
    lines = [f">{chrom}"]
    for i in range(0, len(sequence), line_width):
        lines.append(sequence[i:i + line_width])
    write_text(path, "\n".join(lines) + "\n")


# ===========================================================================
# ValidationReport / Issue
# ===========================================================================


class TestIssue(unittest.TestCase):
    def test_is_namedtuple_with_severity_and_message(self):
        issue = vi.Issue(vi.ERROR, "boom")
        self.assertEqual(issue.severity, "ERROR")
        self.assertEqual(issue.message, "boom")


class TestValidationReport(unittest.TestCase):
    def test_ok_records_an_ok_line_without_affecting_counts(self):
        report = vi.ValidationReport("Title")
        report.ok("all good")
        self.assertFalse(report.has_errors())
        self.assertIn("  [OK] all good", report._lines)

    def test_add_with_empty_issues_and_no_ok_message_adds_nothing(self):
        report = vi.ValidationReport("Title")
        report.add([])
        self.assertEqual(report._lines, [])
        self.assertFalse(report.has_errors())

    def test_add_with_empty_issues_and_ok_message_adds_ok_line(self):
        report = vi.ValidationReport("Title")
        report.add([], ok_message="clean")
        self.assertEqual(report._lines, ["  [OK] clean"])

    def test_add_with_issues_ignores_ok_message_and_updates_counts(self):
        report = vi.ValidationReport("Title")
        report.add([vi.Issue(vi.ERROR, "bad"), vi.Issue(vi.WARN, "meh")],
                    ok_message="should not appear")
        self.assertIn("  [ERROR] bad", report._lines)
        self.assertIn("  [WARN] meh", report._lines)
        self.assertNotIn("  [OK] should not appear", report._lines)
        self.assertTrue(report.has_errors())

    def test_has_errors_false_when_only_warnings(self):
        report = vi.ValidationReport("Title")
        report.add([vi.Issue(vi.WARN, "meh")])
        self.assertFalse(report.has_errors())

    def test_has_errors_true_when_any_error_present(self):
        report = vi.ValidationReport("Title")
        report.add([vi.Issue(vi.WARN, "meh")])
        report.add([vi.Issue(vi.ERROR, "bad")])
        self.assertTrue(report.has_errors())

    def test_render_passed_clean(self):
        report = vi.ValidationReport("My Title")
        report.ok("fine")
        text = report.render()
        self.assertIn("My Title", text)
        self.assertTrue(text.strip().endswith("Result: PASSED"))

    def test_render_passed_with_one_warning_singular(self):
        report = vi.ValidationReport("T")
        report.add([vi.Issue(vi.WARN, "w1")])
        text = report.render()
        self.assertIn("Result: PASSED with 1 warning", text)
        self.assertNotIn("1 warnings", text)

    def test_render_passed_with_multiple_warnings_plural(self):
        report = vi.ValidationReport("T")
        report.add([vi.Issue(vi.WARN, "w1"), vi.Issue(vi.WARN, "w2")])
        text = report.render()
        self.assertIn("Result: PASSED with 2 warnings", text)

    def test_render_failed_one_error_singular(self):
        report = vi.ValidationReport("T")
        report.add([vi.Issue(vi.ERROR, "e1")])
        text = report.render()
        self.assertIn("Result: FAILED (1 error)", text)
        self.assertIn("aborting before launching complete-search", text)

    def test_render_failed_multiple_errors_plural(self):
        report = vi.ValidationReport("T")
        report.add([vi.Issue(vi.ERROR, "e1"), vi.Issue(vi.ERROR, "e2")])
        text = report.render()
        self.assertIn("Result: FAILED (2 errors)", text)

    def test_render_failed_takes_priority_over_warnings(self):
        report = vi.ValidationReport("T")
        report.add([vi.Issue(vi.ERROR, "e1"), vi.Issue(vi.WARN, "w1")])
        text = report.render()
        self.assertIn("FAILED", text)
        self.assertNotIn("PASSED", text)

    def test_write_clean_report_goes_to_stdout_not_stderr(self):
        report = vi.ValidationReport("T")
        report.ok("fine")
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            report.write()
        self.assertIn("Result: PASSED", out.getvalue())
        self.assertEqual(err.getvalue(), "")

    def test_write_warn_only_report_goes_to_stdout_not_stderr(self):
        # WARN is non-fatal — must not land on stderr, which the rest of the
        # pipeline treats as a fatal signal regardless of content (issue #94).
        report = vi.ValidationReport("T")
        report.add([vi.Issue(vi.WARN, "w1")])
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            report.write()
        self.assertIn("[WARN] w1", out.getvalue())
        self.assertEqual(err.getvalue(), "")

    def test_write_error_report_goes_to_stderr_not_stdout(self):
        report = vi.ValidationReport("T")
        report.add([vi.Issue(vi.ERROR, "e1")])
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            report.write()
        self.assertIn("[ERROR] e1", err.getvalue())
        self.assertEqual(out.getvalue(), "")


# ===========================================================================
# _af_field
# ===========================================================================


class TestAfField(unittest.TestCase):
    def test_exact_af_first(self):
        self.assertEqual(vi._af_field("AF=0.1;DP=10"), (0, "AF"))

    def test_exact_af_second(self):
        self.assertEqual(vi._af_field("DP=10;AF=0.1"), (1, "AF"))

    def test_prefixed_non_exact_key(self):
        self.assertEqual(vi._af_field("AF_afr=0.2;DP=10"), (0, "AF_afr"))

    def test_no_af_field(self):
        self.assertEqual(vi._af_field("DP=10;AC=2"), (None, None))


class TestAfValueParsesAsNumeric(unittest.TestCase):
    def test_simple_numeric_value(self):
        self.assertTrue(vi._af_value_parses_as_numeric("AF=0.1;DP=10", 0))

    def test_multiallelic_comma_separated_values(self):
        self.assertTrue(vi._af_value_parses_as_numeric("AF=0.1,0.05", 0))

    def test_missing_value_dot(self):
        self.assertFalse(vi._af_value_parses_as_numeric("AF=.", 0))

    def test_empty_value(self):
        self.assertFalse(vi._af_value_parses_as_numeric("AF=", 0))

    def test_no_equals_sign(self):
        self.assertFalse(vi._af_value_parses_as_numeric("AF", 0))

    def test_one_bad_value_among_multiallelic_fails(self):
        self.assertFalse(vi._af_value_parses_as_numeric("AF=0.1,.", 0))


# ===========================================================================
# check_genome_fasta
# ===========================================================================


class TestCheckGenomeFasta(unittest.TestCase):
    def test_valid_multi_chromosome_dir(self):
        with tempfile.TemporaryDirectory() as tmp:
            make_fasta(os.path.join(tmp, "chr1.fa"), "chr1", "ACGT" * 10)
            make_fasta(os.path.join(tmp, "chr2.fasta"), "chr2", "TTTT" * 10)
            issues, chroms = vi.check_genome_fasta(tmp)
            self.assertEqual(issues, [])
            self.assertEqual(sorted(chroms), ["chr1", "chr2"])

    def test_missing_directory(self):
        issues, chroms = vi.check_genome_fasta("/no/such/genome/dir/xyz")
        self.assertEqual(len(issues), 1)
        self.assertEqual(issues[0].severity, vi.ERROR)
        self.assertIn("does not exist", issues[0].message)
        self.assertEqual(chroms, [])

    def test_empty_directory(self):
        with tempfile.TemporaryDirectory() as tmp:
            issues, chroms = vi.check_genome_fasta(tmp)
            self.assertEqual(len(issues), 1)
            self.assertEqual(issues[0].severity, vi.ERROR)
            self.assertIn("No .fa/.fasta files found", issues[0].message)
            self.assertEqual(chroms, [])

    def test_fasta_not_starting_with_gt(self):
        with tempfile.TemporaryDirectory() as tmp:
            write_text(os.path.join(tmp, "chr1.fa"), "ACGTACGT\n")
            issues, chroms = vi.check_genome_fasta(tmp)
            self.assertEqual(len(issues), 1)
            self.assertEqual(issues[0].severity, vi.ERROR)
            self.assertIn("does not start with '>'", issues[0].message)
            # chrom name still collected even though content is invalid
            self.assertEqual(chroms, ["chr1"])


# ===========================================================================
# check_vcf_filename_chrom
# ===========================================================================


class TestCheckVcfFilenameChrom(unittest.TestCase):
    def test_filename_with_chrn_token(self):
        issues, chrom_by_file = vi.check_vcf_filename_chrom(["/x/chr1.vcf.gz"])
        self.assertEqual(issues, [])
        self.assertEqual(chrom_by_file["chr1.vcf.gz"], "chr1")

    def test_filename_with_chrn_token_not_first(self):
        issues, chrom_by_file = vi.check_vcf_filename_chrom(
            ["/x/dataset.chr2.vcf.gz"]
        )
        self.assertEqual(issues, [])
        self.assertEqual(chrom_by_file["dataset.chr2.vcf.gz"], "chr2")

    def test_filename_without_chrn_token(self):
        issues, chrom_by_file = vi.check_vcf_filename_chrom(["/x/1.vcf.gz"])
        self.assertEqual(len(issues), 1)
        self.assertEqual(issues[0].severity, vi.ERROR)
        self.assertIn("no chrN token", issues[0].message)
        self.assertEqual(chrom_by_file["1.vcf.gz"], "")


# ===========================================================================
# check_vcf_chrom_matches_genome
# ===========================================================================


class TestCheckVcfChromMatchesGenome(unittest.TestCase):
    def test_matching_chromosomes(self):
        issues = vi.check_vcf_chrom_matches_genome(
            {"chr1.vcf.gz": "chr1"}, ["chr1", "chr2"]
        )
        self.assertEqual(issues, [])

    def test_mismatched_chr_prefix_bug_warns(self):
        # the AoU-style bug: VCF token "1" vs. genome chromosome "chr1" --
        # WARN, not ERROR (changed 2026-07-31): a single mismatched
        # chromosome shouldn't block the whole multi-dataset run
        issues = vi.check_vcf_chrom_matches_genome(
            {"1.vcf.gz": "1"}, ["chr1", "chr2"]
        )
        self.assertEqual(len(issues), 1)
        self.assertEqual(issues[0].severity, vi.WARN)
        self.assertIn("1 of 1", issues[0].message)
        self.assertIn("1.vcf.gz", issues[0].message)
        self.assertIn("reference-only", issues[0].message)

    def test_multiple_mismatches_aggregated_into_one_warning(self):
        # a handful of legitimately-absent chromosomes (e.g. no chrY/chrM in
        # the genome build) shouldn't produce a separate issue per file --
        # one aggregated WARN with a count is the point of the 2026-07-31 change
        issues = vi.check_vcf_chrom_matches_genome(
            {"chr1.vcf.gz": "chr1", "chrY.vcf.gz": "chrY", "chrM.vcf.gz": "chrM"},
            ["chr1"],
        )
        self.assertEqual(len(issues), 1)
        self.assertEqual(issues[0].severity, vi.WARN)
        self.assertIn("2 of 3", issues[0].message)

    def test_empty_chrom_token_is_skipped(self):
        # a missing chromosome token ("") is already reported by
        # check_vcf_filename_chrom; this check must not double-report it.
        issues = vi.check_vcf_chrom_matches_genome({"weird.vcf.gz": ""}, ["chr1"])
        self.assertEqual(issues, [])


# ===========================================================================
# check_vcf_content
# ===========================================================================


class TestCheckVcfContent(unittest.TestCase):
    def test_valid_vcf_no_issues(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "chr1.vcf.gz")
            make_vcf_gz(path, [vcf_record()])
            self.assertEqual(vi.check_vcf_content(path), [])

    def test_non_gzip_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "chr1.vcf.gz")
            write_text(path, "not actually gzipped\n")
            issues = vi.check_vcf_content(path)
            self.assertEqual(len(issues), 1)
            self.assertEqual(issues[0].severity, vi.ERROR)
            self.assertIn("not gzip-compressed", issues[0].message)

    def test_missing_chrom_header(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "chr1.vcf.gz")
            write_gzip_text(path, "##fileformat=VCFv4.2\n")
            issues = vi.check_vcf_content(path)
            self.assertEqual(len(issues), 1)
            self.assertEqual(issues[0].severity, vi.ERROR)
            self.assertIn("no #CHROM header line found", issues[0].message)

    def test_sites_only_header_too_few_columns(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "chr1.vcf.gz")
            sites_only_header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
            record = "chr1\t100\t.\tA\tT\t.\tPASS\tAF=0.1"
            make_vcf_gz(path, [record], header=sites_only_header)
            issues = vi.check_vcf_content(path)
            self.assertEqual(len(issues), 1)
            self.assertEqual(issues[0].severity, vi.ERROR)
            self.assertIn("sites-only", issues[0].message)

    def test_missing_af_field_entirely(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "chr1.vcf.gz")
            make_vcf_gz(path, [vcf_record(info="DP=10")])
            issues = vi.check_vcf_content(path)
            self.assertEqual(len(issues), 1)
            self.assertEqual(issues[0].severity, vi.ERROR)
            self.assertIn("no AF-prefixed entry", issues[0].message)

    def test_af_field_present_but_not_exact_key_errors(self):
        # traced end-to-end: enricher.py's fixed 3-character slice garbles
        # the value for any non-exact key, which crashes pd.to_numeric() in
        # CRISPRme_plots.py at the end of the run -- a real crash risk, not
        # just a differently-sourced (but valid) number
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "chr1.vcf.gz")
            make_vcf_gz(path, [vcf_record(info="AF_afr=0.2;AF=0.1")])
            issues = vi.check_vcf_content(path)
            self.assertEqual(len(issues), 1)
            self.assertEqual(issues[0].severity, vi.ERROR)
            self.assertIn("AF_afr=", issues[0].message)

    def test_af_value_non_numeric_errors(self):
        # exact "AF" key, but a missing/non-numeric value -- same downstream
        # pd.to_numeric() crash, different root cause than the wrong-key case
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "chr1.vcf.gz")
            make_vcf_gz(path, [vcf_record(info="AF=.")])
            issues = vi.check_vcf_content(path)
            self.assertEqual(len(issues), 1)
            self.assertEqual(issues[0].severity, vi.ERROR)
            self.assertIn("non-numeric or missing", issues[0].message)

    def test_af_value_multiallelic_numeric_passes(self):
        # comma-separated multi-allelic AF values are legitimate, not an error
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "chr1.vcf.gz")
            make_vcf_gz(path, [vcf_record(info="AF=0.1,0.05")])
            issues = vi.check_vcf_content(path)
            self.assertEqual(issues, [])

    def test_malformed_short_first_record(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "chr1.vcf.gz")
            make_vcf_gz(path, ["chr1\t100\t."])
            issues = vi.check_vcf_content(path)
            self.assertEqual(len(issues), 1)
            self.assertEqual(issues[0].severity, vi.ERROR)
            self.assertIn("malformed row", issues[0].message)

    def test_empty_vcf_header_only_passes(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "chr1.vcf.gz")
            make_vcf_gz(path, [])
            self.assertEqual(vi.check_vcf_content(path), [])


# ===========================================================================
# check_tbi_files
# ===========================================================================


class TestCheckTbiFiles(unittest.TestCase):
    def test_stray_tbi_present_warns(self):
        with tempfile.TemporaryDirectory() as tmp:
            write_text(os.path.join(tmp, "chr1.vcf.gz.tbi"), "")
            issues = vi.check_tbi_files(tmp)
            self.assertEqual(len(issues), 1)
            self.assertEqual(issues[0].severity, vi.WARN)
            self.assertIn("chr1.vcf.gz.tbi", issues[0].message)

    def test_no_tbi_files_is_clean(self):
        with tempfile.TemporaryDirectory() as tmp:
            write_text(os.path.join(tmp, "chr1.vcf.gz"), "")
            self.assertEqual(vi.check_tbi_files(tmp), [])


# ===========================================================================
# check_pam_file
# ===========================================================================


class TestCheckPamFile(unittest.TestCase):
    def test_valid_pam(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "pam.txt")
            write_text(path, "NGG 3\n")
            self.assertEqual(vi.check_pam_file(path), [])

    def test_missing_offset_token(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "pam.txt")
            write_text(path, "NGG\n")
            issues = vi.check_pam_file(path)
            self.assertEqual(len(issues), 1)
            self.assertEqual(issues[0].severity, vi.ERROR)
            self.assertIn("1 whitespace-separated token", issues[0].message)

    def test_empty_file_zero_tokens(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "pam.txt")
            write_text(path, "")
            issues = vi.check_pam_file(path)
            self.assertEqual(len(issues), 1)
            self.assertEqual(issues[0].severity, vi.ERROR)
            self.assertIn("0 whitespace-separated token", issues[0].message)

    def test_non_integer_offset(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "pam.txt")
            write_text(path, "NGG abc\n")
            issues = vi.check_pam_file(path)
            self.assertEqual(len(issues), 1)
            self.assertEqual(issues[0].severity, vi.ERROR)
            self.assertIn("not a valid integer", issues[0].message)

    def test_invalid_iupac_character(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "pam.txt")
            write_text(path, "NGGZ 3\n")
            issues = vi.check_pam_file(path)
            self.assertEqual(len(issues), 1)
            self.assertEqual(issues[0].severity, vi.ERROR)
            self.assertIn("invalid IUPAC character", issues[0].message)
            self.assertIn("Z", issues[0].message)


# ===========================================================================
# check_guide_file
# ===========================================================================


class TestCheckGuideFile(unittest.TestCase):
    def test_uniform_length_guides(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "guides.txt")
            write_text(path, "ACGTACGTACGTACGTACGT\nTTTTACGTACGTACGTACGT\n")
            self.assertEqual(vi.check_guide_file(path), [])

    def test_mixed_length_guides(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "guides.txt")
            write_text(path, "ACGTACGTACGTACGTACGT\nACGTACGT\n")
            issues = vi.check_guide_file(path)
            self.assertEqual(len(issues), 1)
            self.assertEqual(issues[0].severity, vi.ERROR)
            self.assertIn("inconsistent lengths", issues[0].message)

    def test_empty_guide_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "guides.txt")
            write_text(path, "")
            issues = vi.check_guide_file(path)
            self.assertEqual(len(issues), 1)
            self.assertEqual(issues[0].severity, vi.ERROR)
            self.assertIn("guide file is empty", issues[0].message)


# ===========================================================================
# check_gzip_compressed
# ===========================================================================


class TestCheckGzipCompressed(unittest.TestCase):
    def test_gzip_magic_passes(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "gene_annotation.bed.gz")
            write_gzip_text(path, "chr1\t0\t100\tfeature\n")
            self.assertEqual(vi.check_gzip_compressed(path, "gene annotation file"), [])

    def test_plain_text_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "gene_annotation.bed.gz")
            write_text(path, "chr1\t0\t100\tfeature\n")
            issues = vi.check_gzip_compressed(path, "gene annotation file")
            self.assertEqual(len(issues), 1)
            self.assertEqual(issues[0].severity, vi.ERROR)
            self.assertIn("not gzip-compressed", issues[0].message)
            self.assertIn("gene annotation file", issues[0].message)


# ===========================================================================
# check_samples_in_vcf_header
# ===========================================================================


class TestCheckSamplesInVcfHeader(unittest.TestCase):
    HEADER = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2"

    def test_all_samples_present_is_clean(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "chr1.vcf.gz")
            make_vcf_gz(path, [vcf_record(samples=("0/1", "1/1"))], header=self.HEADER)
            issues = vi.check_samples_in_vcf_header(path, ["S1", "S2", "S3"])
            self.assertEqual(issues, [])

    def test_missing_vcf_sample_errors(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "chr1.vcf.gz")
            make_vcf_gz(path, [vcf_record(samples=("0/1", "1/1"))], header=self.HEADER)
            issues = vi.check_samples_in_vcf_header(path, ["S1"])
            self.assertEqual(len(issues), 1)
            self.assertEqual(issues[0].severity, vi.ERROR)
            self.assertIn("S2", issues[0].message)
            self.assertIn("KeyError", issues[0].message)

    def test_empty_sample_ids_is_skipped_not_flagged(self):
        # #112 regression: when no sample IDs could be resolved from
        # --samplesID, the check must SKIP (empty result), not flag every VCF
        # sample as missing.
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "chr1.vcf.gz")
            make_vcf_gz(path, [vcf_record(samples=("0/1", "1/1"))], header=self.HEADER)
            self.assertEqual(vi.check_samples_in_vcf_header(path, []), [])


# ===========================================================================
# _read_samples_file  (resolves a --samplesID *config* into real sample IDs)
# ===========================================================================


class TestReadSamplesFile(unittest.TestCase):
    def _make_cwd_with_assoc(self, tmp, config_name, assoc_files):
        """Build <tmp>/samplesIDs/<assoc> files + a --samplesID config listing
        their names. Returns the config path. `assoc_files` maps
        filename -> file text."""
        samplesids_dir = os.path.join(tmp, "samplesIDs")
        os.makedirs(samplesids_dir, exist_ok=True)
        for fname, text in assoc_files.items():
            write_text(os.path.join(samplesids_dir, fname), text)
        config_path = os.path.join(tmp, config_name)
        write_text(config_path, "".join(f"{n}\n" for n in assoc_files))
        return config_path

    def test_resolves_association_file_and_reads_sample_ids(self):
        # The 1000G shape: --samplesID lists an association filename; the
        # association file has SAMPLE_ID/POPULATION_ID/... columns.
        with tempfile.TemporaryDirectory() as tmp:
            assoc = "#SAMPLE_ID\tPOPULATION_ID\tSUPERPOPULATION_ID\tSEX\nHG00096\tGBR\tEUR\tmale\nHG00097\tGBR\tEUR\tfemale\n"
            config = self._make_cwd_with_assoc(
                tmp, "samplesIDs.config.txt", {"hg38_1000G.samplesID.txt": assoc}
            )
            self.assertEqual(
                vi._read_samples_file(config, tmp), ["HG00096", "HG00097"]
            )

    def test_unions_multiple_referenced_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            config = self._make_cwd_with_assoc(
                tmp,
                "cfg.txt",
                {
                    "a.txt": "#SAMPLE_ID\tP\tS\tX\nS1\tp\ts\tm\n",
                    "b.txt": "#SAMPLE_ID\tP\tS\tX\nS2\tp\ts\tf\n",
                },
            )
            self.assertEqual(sorted(vi._read_samples_file(config, tmp)), ["S1", "S2"])

    def test_missing_referenced_file_is_skipped(self):
        # If a referenced association file can't be resolved yet, skip it (do
        # NOT treat the filename as a sample ID -- that was the #112 bug).
        with tempfile.TemporaryDirectory() as tmp:
            config = os.path.join(tmp, "cfg.txt")
            write_text(config, "hg38_1000G.samplesID.txt\n")
            self.assertEqual(vi._read_samples_file(config, tmp), [])

    def test_resolves_relative_to_config_dir_when_no_cwd(self):
        with tempfile.TemporaryDirectory() as tmp:
            write_text(
                os.path.join(tmp, "assoc.txt"),
                "#SAMPLE_ID\tP\tS\tX\nS1\tp\ts\tm\n",
            )
            config = os.path.join(tmp, "cfg.txt")
            write_text(config, "assoc.txt\n")
            self.assertEqual(vi._read_samples_file(config), ["S1"])


# ===========================================================================
# _read_vcf_dataset_names / resolve_vcf_dataset_dirs
# ===========================================================================


class TestReadVcfDatasetNames(unittest.TestCase):
    def test_plain_names(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "vcf_config.txt")
            write_text(path, "datasetA\ndatasetB\n")
            self.assertEqual(vi._read_vcf_dataset_names(path), ["datasetA", "datasetB"])

    def test_trailing_slash_is_stripped(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "vcf_config.txt")
            write_text(path, "datasetA/\n")
            self.assertEqual(vi._read_vcf_dataset_names(path), ["datasetA"])


class TestResolveVcfDatasetDirs(unittest.TestCase):
    def test_resolves_under_vcfs_subdir(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "vcf_config.txt")
            write_text(path, "datasetA\ndatasetB\n")
            dirs = vi.resolve_vcf_dataset_dirs(path, "/cwd")
            self.assertEqual(
                dirs, [os.path.join("/cwd", "VCFs", "datasetA"), os.path.join("/cwd", "VCFs", "datasetB")]
            )


# ===========================================================================
# run_lightweight end-to-end
# ===========================================================================


class TestRunLightweightEndToEnd(unittest.TestCase):
    def _build_clean_fixture(self, tmp):
        genomedir = os.path.join(tmp, "genome")
        os.makedirs(genomedir)
        make_fasta(os.path.join(genomedir, "chr1.fa"), "chr1", "ACGT" * 50)
        make_fasta(os.path.join(genomedir, "chr2.fa"), "chr2", "TTTT" * 50)

        vcf_dir = os.path.join(tmp, "VCFs", "datasetA")
        os.makedirs(vcf_dir)
        make_vcf_gz(os.path.join(vcf_dir, "chr1.vcf.gz"), [vcf_record()])

        pamfile = os.path.join(tmp, "pam.txt")
        write_text(pamfile, "NGG 3\n")

        guidefile = os.path.join(tmp, "guides.txt")
        write_text(guidefile, "ACGTACGTACGTACGTACGT\n")

        return genomedir, [vcf_dir], guidefile, pamfile

    def test_fully_clean_fixture_has_no_errors(self):
        with tempfile.TemporaryDirectory() as tmp:
            genomedir, vcf_dirs, guidefile, pamfile = self._build_clean_fixture(tmp)
            report = vi.run_lightweight(genomedir, vcf_dirs, guidefile, pamfile)
            self.assertFalse(report.has_errors(), report.render())

    def test_genome_vcf_chrom_mismatch_warns_not_errors(self):
        # changed 2026-07-31: a single mismatched chromosome is now WARN,
        # not a run-blocking ERROR (see check_vcf_chrom_matches_genome)
        with tempfile.TemporaryDirectory() as tmp:
            genomedir, vcf_dirs, guidefile, pamfile = self._build_clean_fixture(tmp)
            # rename the vcf so its chrom token (chr3) has no matching genome FASTA
            vcf_dir = vcf_dirs[0]
            os.rename(
                os.path.join(vcf_dir, "chr1.vcf.gz"),
                os.path.join(vcf_dir, "chr3.vcf.gz"),
            )
            report = vi.run_lightweight(genomedir, vcf_dirs, guidefile, pamfile)
            self.assertFalse(report.has_errors())
            self.assertIn("reference-only", report.render())


# ===========================================================================
# _normalize_chrom
# ===========================================================================


class TestNormalizeChrom(unittest.TestCase):
    def test_unprefixed_gets_prefix(self):
        self.assertEqual(vi._normalize_chrom("1"), "chr1")

    def test_already_prefixed_unchanged(self):
        self.assertEqual(vi._normalize_chrom("chr1"), "chr1")


# ===========================================================================
# _genome_chrom_length
# ===========================================================================


class TestGenomeChromLength(unittest.TestCase):
    def test_computes_length_from_fasta(self):
        with tempfile.TemporaryDirectory() as tmp:
            make_fasta(os.path.join(tmp, "chr1.fa"), "chr1", "A" * 137)
            self.assertEqual(vi._genome_chrom_length(tmp, "chr1"), 137)

    def test_returns_none_when_no_matching_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.assertIsNone(vi._genome_chrom_length(tmp, "chr9"))


# ===========================================================================
# check_vcf_full_scan
# ===========================================================================


class TestCheckVcfFullScan(unittest.TestCase):
    def _scan(self, tmp, records, expected_chrom="chr1", chrom_length=None, header=DEFAULT_VCF_HEADER):
        path = os.path.join(tmp, "chr1.vcf.gz")
        make_vcf_gz(path, records, header=header)
        return vi.check_vcf_full_scan(path, expected_chrom, chrom_length)

    def test_empty_vcf_header_only_passes(self):
        with tempfile.TemporaryDirectory() as tmp:
            issues = self._scan(tmp, [])
            self.assertEqual(issues, [])

    def test_clean_file_no_issues(self):
        with tempfile.TemporaryDirectory() as tmp:
            records = [
                vcf_record(pos=100, ref="A", alt="T"),
                vcf_record(pos=200, ref="A", alt="T"),
            ]
            issues = self._scan(tmp, records)
            self.assertEqual(issues, [])

    def test_long_structural_variant_warns(self):
        with tempfile.TemporaryDirectory() as tmp:
            long_ref = "A" * 60
            records = [
                vcf_record(pos=100, ref="A", alt="T"),
                vcf_record(pos=200, ref=long_ref, alt="T"),
            ]
            issues = self._scan(tmp, records)
            self.assertEqual(len(issues), 1)
            self.assertEqual(issues[0].severity, vi.WARN)
            self.assertIn(">50bp", issues[0].message)

    def test_rogue_chromosome_errors(self):
        with tempfile.TemporaryDirectory() as tmp:
            records = [
                vcf_record(chrom="chr1", pos=100, ref="A", alt="T"),
                vcf_record(chrom="chr2", pos=150, ref="A", alt="T"),
            ]
            issues = self._scan(tmp, records, expected_chrom="chr1")
            self.assertEqual(len(issues), 1)
            self.assertEqual(issues[0].severity, vi.ERROR)
            self.assertIn("CHROM value other than", issues[0].message)
            self.assertIn("chr2", issues[0].message)

    def test_af_position_drift_warns(self):
        with tempfile.TemporaryDirectory() as tmp:
            records = [
                vcf_record(pos=100, info="AF=0.1;DP=10"),
                vcf_record(pos=200, info="DP=10;AF=0.2"),
            ]
            issues = self._scan(tmp, records)
            self.assertEqual(len(issues), 1)
            self.assertEqual(issues[0].severity, vi.WARN)
            self.assertIn("AF field position varies", issues[0].message)

    def test_zero_pass_rate_errors(self):
        # unlike the partial <10% case below, 0% is a total, guaranteed
        # void of the whole dataset's contribution to enrichment
        with tempfile.TemporaryDirectory() as tmp:
            records = [
                vcf_record(pos=100, filt="FAIL"),
                vcf_record(pos=200, filt="LowQual"),
            ]
            issues = self._scan(tmp, records)
            self.assertEqual(len(issues), 1)
            self.assertEqual(issues[0].severity, vi.ERROR)
            self.assertIn("0 of 2 variants", issues[0].message)

    def test_low_pass_rate_below_threshold_warns(self):
        with tempfile.TemporaryDirectory() as tmp:
            records = [vcf_record(pos=100 + i, filt="PASS" if i == 0 else "FAIL")
                       for i in range(20)]
            issues = self._scan(tmp, records)
            self.assertEqual(len(issues), 1)
            self.assertEqual(issues[0].severity, vi.WARN)
            self.assertIn("only 1/20", issues[0].message)
            self.assertIn("5.0%", issues[0].message)

    def test_healthy_pass_rate_no_warning(self):
        with tempfile.TemporaryDirectory() as tmp:
            records = [vcf_record(pos=100 + i, filt="PASS" if i < 9 else "FAIL")
                       for i in range(10)]
            issues = self._scan(tmp, records)
            self.assertEqual(issues, [])

    def test_dot_filter_not_treated_as_passing(self):
        # the officially-released crispritz enricher.py only ever hardcodes
        # 'PASS' (verified against GitHub pinellolab/CRISPRitz master); '.'
        # is accepted only by an unsubmitted local patch in one developer's
        # own conda env for HPRC testing, not any shipped crispritz release,
        # so the validator must not treat '.' as passing today
        with tempfile.TemporaryDirectory() as tmp:
            records = [
                vcf_record(pos=100, filt="."),
                vcf_record(pos=200, filt="."),
            ]
            issues = self._scan(tmp, records)
            self.assertEqual(len(issues), 1)
            self.assertEqual(issues[0].severity, vi.ERROR)
            self.assertIn("0 of 2 variants", issues[0].message)
            self.assertIn("FILTER == 'PASS'", issues[0].message)

    def test_pos_out_of_bounds_errors(self):
        with tempfile.TemporaryDirectory() as tmp:
            records = [
                vcf_record(pos=500, ref="A", alt="T"),
                vcf_record(pos=2000, ref="A", alt="T"),
            ]
            issues = self._scan(tmp, records, chrom_length=1000)
            self.assertEqual(len(issues), 1)
            self.assertEqual(issues[0].severity, vi.ERROR)
            self.assertIn("outside", issues[0].message)
            self.assertIn("IndexError", issues[0].message)

    def test_pos_within_bounds_no_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            records = [vcf_record(pos=500, ref="A", alt="T")]
            issues = self._scan(tmp, records, chrom_length=1000)
            self.assertEqual(issues, [])

    def test_high_allele_count_multiallelic_warns(self):
        with tempfile.TemporaryDirectory() as tmp:
            alleles = ["C", "G", "T", "C", "G", "T", "C", "G", "T", "C"]  # 10 alleles
            alt = ",".join(alleles)
            af = ",".join(["0.01"] * len(alleles))  # enough AF values to isolate this check
            records = [vcf_record(pos=1000, ref="A", alt=alt, info=f"AF={af}")]
            issues = self._scan(tmp, records)
            self.assertEqual(len(issues), 1)
            self.assertEqual(issues[0].severity, vi.WARN)
            self.assertIn(">=10 ALT alleles", issues[0].message)

    def test_below_threshold_allele_count_no_warning(self):
        with tempfile.TemporaryDirectory() as tmp:
            alleles = ["C", "G", "T", "C", "G", "T", "C", "G", "T"]  # 9 alleles
            alt = ",".join(alleles)
            af = ",".join(["0.01"] * len(alleles))  # enough AF values to isolate this check
            records = [vcf_record(pos=1000, ref="A", alt=alt, info=f"AF={af}")]
            issues = self._scan(tmp, records)
            self.assertEqual(issues, [])

    def test_af_index_overflow_errors(self):
        # 2 single-char SNP ALT alleles (positions 1, 2) but only 1 AF value
        # -- enricher.py indexes af[2-1] = af[1], an uncaught IndexError
        with tempfile.TemporaryDirectory() as tmp:
            records = [vcf_record(pos=1000, ref="A", alt="C,G", info="AF=0.1")]
            issues = self._scan(tmp, records)
            self.assertEqual(len(issues), 1)
            self.assertEqual(issues[0].severity, vi.ERROR)
            self.assertIn("fewer comma-separated AF values", issues[0].message)

    def test_af_index_sufficient_values_passes(self):
        with tempfile.TemporaryDirectory() as tmp:
            records = [vcf_record(pos=1000, ref="A", alt="C,G", info="AF=0.1,0.05")]
            issues = self._scan(tmp, records)
            self.assertEqual(issues, [])

    def test_af_index_indel_allele_does_not_consume_af_slot(self):
        # second ALT allele is a multi-char indel (not a "SNP" in
        # enricher.py's sense), so it never gets indexed into the AF list --
        # only the single-char "C" at position 1 does, and 1 AF value covers it
        with tempfile.TemporaryDirectory() as tmp:
            records = [vcf_record(pos=1000, ref="A", alt="C,ATCG", info="AF=0.1")]
            issues = self._scan(tmp, records)
            self.assertEqual(issues, [])

    def test_af_index_biallelic_single_alt_not_checked(self):
        # no comma in ALT at all -- enricher.py's biallelic branch uses the
        # whole AF string directly, never splits/indexes it, so a "short" AF
        # here isn't actually at risk
        with tempfile.TemporaryDirectory() as tmp:
            records = [vcf_record(pos=1000, ref="A", alt="T", info="AF=0.1")]
            issues = self._scan(tmp, records)
            self.assertEqual(issues, [])

    def test_breakend_notation_errors(self):
        with tempfile.TemporaryDirectory() as tmp:
            records = [vcf_record(pos=1000, ref="G", alt="G]chr2:500]")]
            issues = self._scan(tmp, records)
            self.assertEqual(len(issues), 1)
            self.assertEqual(issues[0].severity, vi.ERROR)
            self.assertIn("breakend/SV notation", issues[0].message)

    def test_symbolic_del_alt_does_not_trigger_breakend(self):
        with tempfile.TemporaryDirectory() as tmp:
            records = [vcf_record(pos=1000, ref="A", alt="<DEL>")]
            issues = self._scan(tmp, records)
            self.assertEqual(len(issues), 1)
            self.assertNotIn("breakend", issues[0].message)

    def test_symbolic_alt_warns(self):
        # different failure mode than breakend: silently excluded from the
        # search entirely, not corrupted -- WARN, not ERROR
        with tempfile.TemporaryDirectory() as tmp:
            records = [vcf_record(pos=1000, ref="A", alt="<DEL>")]
            issues = self._scan(tmp, records)
            self.assertEqual(len(issues), 1)
            self.assertEqual(issues[0].severity, vi.WARN)
            self.assertIn("symbolic ALT notation", issues[0].message)
            self.assertIn("never searched", issues[0].message)

    def test_duplicate_chr_pos_warns(self):
        with tempfile.TemporaryDirectory() as tmp:
            records = [
                vcf_record(pos=100, ref="A", alt="T"),
                vcf_record(pos=100, ref="A", alt="C"),
            ]
            issues = self._scan(tmp, records)
            self.assertEqual(len(issues), 1)
            self.assertEqual(issues[0].severity, vi.WARN)
            self.assertIn("duplicate CHR+POS", issues[0].message)

    def test_no_duplicates_no_warning(self):
        with tempfile.TemporaryDirectory() as tmp:
            records = [
                vcf_record(pos=100, ref="A", alt="T"),
                vcf_record(pos=200, ref="A", alt="C"),
            ]
            issues = self._scan(tmp, records)
            self.assertEqual(issues, [])

    def test_mixed_phasing_warns(self):
        with tempfile.TemporaryDirectory() as tmp:
            records = [
                vcf_record(pos=100, samples=("0/1",)),
                vcf_record(pos=200, samples=("0|1",)),
            ]
            issues = self._scan(tmp, records)
            self.assertEqual(len(issues), 1)
            self.assertEqual(issues[0].severity, vi.WARN)
            self.assertIn("both phased", issues[0].message)

    def test_only_unphased_no_warning(self):
        with tempfile.TemporaryDirectory() as tmp:
            records = [
                vcf_record(pos=100, samples=("0/1",)),
                vcf_record(pos=200, samples=("1/1",)),
            ]
            issues = self._scan(tmp, records)
            self.assertEqual(issues, [])

    def test_only_phased_no_warning(self):
        with tempfile.TemporaryDirectory() as tmp:
            records = [
                vcf_record(pos=100, samples=("0|1",)),
                vcf_record(pos=200, samples=("1|1",)),
            ]
            issues = self._scan(tmp, records)
            self.assertEqual(issues, [])

    def test_indel_near_chromosome_start_warns(self):
        with tempfile.TemporaryDirectory() as tmp:
            records = [vcf_record(pos=10, ref="A", alt="AT")]
            issues = self._scan(tmp, records)
            self.assertEqual(len(issues), 1)
            self.assertEqual(issues[0].severity, vi.WARN)
            self.assertIn("chromosome start", issues[0].message)

    def test_indel_far_from_start_no_warning(self):
        with tempfile.TemporaryDirectory() as tmp:
            records = [vcf_record(pos=1000, ref="A", alt="AT")]
            issues = self._scan(tmp, records)
            self.assertEqual(issues, [])


# ===========================================================================
# run_full end-to-end
# ===========================================================================


class TestRunFullEndToEnd(unittest.TestCase):
    def test_clean_fixture_has_no_errors(self):
        with tempfile.TemporaryDirectory() as tmp:
            genomedir = os.path.join(tmp, "genome")
            os.makedirs(genomedir)
            make_fasta(os.path.join(genomedir, "chr1.fa"), "chr1", "A" * 5000)

            vcf_dir = os.path.join(tmp, "VCFs", "datasetA")
            os.makedirs(vcf_dir)
            records = [
                vcf_record(chrom="chr1", pos=1000, ref="A", alt="T"),
                vcf_record(chrom="chr1", pos=2000, ref="A", alt="C"),
            ]
            make_vcf_gz(os.path.join(vcf_dir, "chr1.vcf.gz"), records)

            report = vi.run_full(genomedir, [vcf_dir])
            self.assertFalse(report.has_errors(), report.render())

    def test_defective_fixture_has_errors(self):
        with tempfile.TemporaryDirectory() as tmp:
            genomedir = os.path.join(tmp, "genome")
            os.makedirs(genomedir)
            make_fasta(os.path.join(genomedir, "chr1.fa"), "chr1", "A" * 5000)

            vcf_dir = os.path.join(tmp, "VCFs", "datasetA")
            os.makedirs(vcf_dir)
            records = [
                # rogue chromosome (ERROR)
                vcf_record(chrom="chr9", pos=1000, ref="A", alt="T"),
                # breakend notation (ERROR)
                vcf_record(chrom="chr1", pos=2000, ref="G", alt="G]chr2:500]"),
                # duplicate CHR+POS (WARN)
                vcf_record(chrom="chr1", pos=3000, ref="A", alt="T"),
                vcf_record(chrom="chr1", pos=3000, ref="A", alt="C"),
                # indel near chromosome start (WARN)
                vcf_record(chrom="chr1", pos=5, ref="A", alt="AT"),
            ]
            make_vcf_gz(os.path.join(vcf_dir, "chr1.vcf.gz"), records)

            report = vi.run_full(genomedir, [vcf_dir])
            self.assertTrue(report.has_errors())
            rendered = report.render()
            self.assertIn("breakend/SV notation", rendered)
            self.assertIn("CHROM value other than", rendered)


if __name__ == "__main__":
    unittest.main()
