#!/usr/bin/env python
"""Filter raw CRISPRitz targets to those whose PAM region can satisfy the
requested PAM (IUPAC-aware).

CRISPRitz PAM-specific indexes prune to PAM-matching sites at index time, so
their targets already match. A degenerate/pamless (NNN) index does not, so a
search reusing it for a specific PAM (e.g. NGG) would otherwise report targets
with a non-matching PAM (e.g. TGC). This filter removes those, making results
correct regardless of the index used, and is a no-op for PAM-specific and direct
(-r) searches.

The match is IUPAC-aware (keep a target if *some* allele of its PAM region
satisfies the requested PAM), so variants that CREATE the PAM on the enriched
genome are preserved for the downstream PAM-creation analysis.

Usage: pam_filter.py <targets.txt> <pam_file>   (edits the targets file in place)
"""

import os
import sys

_IUPAC = {
    "A": set("A"), "C": set("C"), "G": set("G"), "T": set("T"),
    "R": set("AG"), "Y": set("CT"), "S": set("GC"), "W": set("AT"),
    "K": set("GT"), "M": set("AC"), "B": set("CGT"), "D": set("AGT"),
    "H": set("ACT"), "V": set("ACG"), "N": set("ACGT"),
}


def _pam_can_match(target_pam: str, requested_pam: str) -> bool:
    """True if the target PAM region can satisfy the requested PAM at every
    position (IUPAC intersection non-empty)."""
    if len(target_pam) != len(requested_pam):
        return False
    for tc, rc in zip(target_pam.upper(), requested_pam.upper()):
        if not (_IUPAC.get(tc, set()) & _IUPAC.get(rc, set())):
            return False
    return True


def _requested_pam(pam_file: str):
    """Return (motif, three_prime) for the requested PAM, mirroring how
    index-genome / build-index-only extract it from the PAM file."""
    with open(pam_file) as handle:
        first = handle.readline().split()
    seq, pos = first[0], int(first[1])
    n = abs(pos)
    if pos < 0:  # 5' PAM (e.g. Cas12a): motif at the start
        return seq[:n], False
    return seq[-n:], True  # 3' PAM (e.g. SpCas9): motif at the end


def main() -> None:
    targets_file, pam_file = sys.argv[1], sys.argv[2]
    motif, three_prime = _requested_pam(pam_file)
    n = len(motif)
    tmp = targets_file + ".pamfilt"
    kept = dropped = 0
    with open(targets_file) as fin, open(tmp, "w") as fout:
        for row in fin:
            cols = row.rstrip("\n").split("\t")
            if len(cols) < 3:  # header or malformed -> pass through untouched
                fout.write(row)
                continue
            target = cols[2]  # aligned protospacer+PAM (CRISPRitz column 3)
            tgt_pam = target[-n:] if three_prime else target[:n]
            if _pam_can_match(tgt_pam, motif):
                fout.write(row)
                kept += 1
            else:
                dropped += 1
    os.replace(tmp, targets_file)
    # write progress to STDOUT, never STDERR: the caller treats a non-empty
    # stderr log as a fatal post-analysis failure.
    sys.stdout.write(
        f"pam_filter: {os.path.basename(targets_file)} kept {kept}, "
        f"dropped {dropped} non-{motif} target(s)\n"
    )


if __name__ == "__main__":
    main()
