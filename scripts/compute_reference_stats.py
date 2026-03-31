#!/usr/bin/env python3
"""
compute_reference_stats.py
--------------------------
Parses cdr3TestCases from test.js, computes 6 physicochemical metrics for
every VH and VL sequence, then writes the resulting APPROVED_MAB_STATS
constant directly into index.html between sentinel comments.

Usage (from repository root):
    python scripts/compute_reference_stats.py
"""

import re
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TEST_JS = ROOT / "test.js"
INDEX_HTML = ROOT / "index.html"

BEGIN_SENTINEL = "// BEGIN APPROVED_MAB_STATS"
END_SENTINEL = "// END APPROVED_MAB_STATS"

# ---------------------------------------------------------------------------
# Same constants as index.html
# ---------------------------------------------------------------------------

AA_RESIDUE_MW = {
    "A": 71.0788,  "R": 156.1875, "N": 114.1038, "D": 115.0886, "C": 103.1388,
    "E": 129.1155, "Q": 128.1307, "G": 57.0519,  "H": 137.1411, "I": 113.1594,
    "L": 113.1594, "K": 128.1741, "M": 131.1926, "F": 147.1766, "P": 97.1167,
    "S": 87.0782,  "T": 101.1051, "W": 186.2132, "Y": 163.1760, "V": 99.1326,
}
WATER_MW = 18.0153

PK = {
    "nTerm": 8.00, "cTerm": 3.10,
    "D": 3.65, "E": 4.25, "C": 8.18, "Y": 10.07,
    "H": 6.00, "K": 10.53, "R": 12.48,
}

KD_HYDROPHOBICITY = {
    "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5, "M": 1.9, "A": 1.8,
    "G": -0.4, "T": -0.7, "S": -0.8, "W": -0.9, "Y": -1.3, "P": -1.6, "H": -3.2,
    "D": -3.5, "E": -3.5, "N": -3.5, "Q": -3.5, "K": -3.9, "R": -4.5,
}


def calc_mw(seq):
    mw = WATER_MW + sum(AA_RESIDUE_MW.get(aa, 0) for aa in seq)
    cys_count = seq.count("C")
    mw -= (cys_count // 2) * 2.016
    return mw


def calc_net_charge(seq, pH):
    charge = 0.0
    charge += 10 ** PK["nTerm"] / (10 ** PK["nTerm"] + 10 ** pH)
    charge -= 10 ** pH / (10 ** PK["cTerm"] + 10 ** pH)
    cys_count = 0
    for aa in seq:
        if aa == "K":
            charge += 10 ** PK["K"] / (10 ** PK["K"] + 10 ** pH)
        elif aa == "R":
            charge += 10 ** PK["R"] / (10 ** PK["R"] + 10 ** pH)
        elif aa == "H":
            charge += 10 ** PK["H"] / (10 ** PK["H"] + 10 ** pH)
        elif aa == "D":
            charge -= 10 ** pH / (10 ** PK["D"] + 10 ** pH)
        elif aa == "E":
            charge -= 10 ** pH / (10 ** PK["E"] + 10 ** pH)
        elif aa == "C":
            cys_count += 1
        elif aa == "Y":
            charge -= 10 ** pH / (10 ** PK["Y"] + 10 ** pH)
            
    free_cys = cys_count % 2
    charge -= free_cys * (10 ** pH / (10 ** PK["C"] + 10 ** pH))
    return charge


def calc_pi(seq):
    lo, hi = 0.0, 14.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if calc_net_charge(seq, mid) > 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def calc_gravy(seq):
    return sum(KD_HYDROPHOBICITY.get(aa, 0) for aa in seq) / len(seq)


def compute_props(seq):
    return {
        "length": len(seq),
        "mw": calc_mw(seq),
        "pI": calc_pi(seq),
        "charge74": calc_net_charge(seq, 7.4),
        "charge55": calc_net_charge(seq, 5.5),
        "gravy": calc_gravy(seq),
    }


# ---------------------------------------------------------------------------
# Parse test.js
# ---------------------------------------------------------------------------

def parse_test_cases(path):
    text = path.read_text()
    # Find the cdr3TestCases array block
    start = text.index("const cdr3TestCases = [")
    end = text.index("];", start) + 2
    block = text[start:end]

    # Extract each object: pull id, chain, seq
    pattern = re.compile(
        r'\{\s*id:\s*"([^"]+)".*?chain:\s*"([^"]+)".*?seq:\s*"([^"]+)"',
        re.DOTALL,
    )
    cases = []
    for m in pattern.finditer(block):
        cases.append({"id": m.group(1), "chain": m.group(2), "seq": m.group(3)})
    return cases


# ---------------------------------------------------------------------------
# Distribution stats helpers
# ---------------------------------------------------------------------------

def percentile(sorted_vals, p):
    """Linear interpolation percentile (same as numpy default)."""
    n = len(sorted_vals)
    if n == 0:
        return 0
    idx = (p / 100) * (n - 1)
    lo = int(math.floor(idx))
    hi = int(math.ceil(idx))
    if lo == hi:
        return sorted_vals[lo]
    return sorted_vals[lo] + (idx - lo) * (sorted_vals[hi] - sorted_vals[lo])


def distribution_stats(values):
    sv = sorted(values)
    return {
        "q1": round(percentile(sv, 25), 4),
        "q2": round(percentile(sv, 50), 4),
        "q3": round(percentile(sv, 75), 4),
        "values": [round(v, 4) for v in sv],
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    cases = parse_test_cases(TEST_JS)
    print(f"Parsed {len(cases)} sequences from test.js", file=sys.stderr)

    by_chain = {"VH": [], "VL": []}
    for c in cases:
        chain = c["chain"]
        if chain in by_chain:
            by_chain[chain].append(c["seq"])

    for chain, seqs in by_chain.items():
        print(f"  {chain}: {len(seqs)} sequences", file=sys.stderr)

    metric_keys = ["length", "mw", "pI", "charge74", "charge55", "gravy"]

    stats = {}
    for chain, seqs in by_chain.items():
        all_props = [compute_props(s) for s in seqs]
        stats[chain] = {"n": len(seqs)}
        for key in metric_keys:
            vals = [p[key] for p in all_props]
            stats[chain][key] = distribution_stats(vals)

    # Serialize compactly: values arrays on one line per metric
    # Build the JS constant string
    lines = ["const APPROVED_MAB_STATS = " + json.dumps(stats, separators=(", ", ": ")) + ";"]
    js_block = "\n".join(lines)

    # Write into index.html between sentinels
    html = INDEX_HTML.read_text()
    begin_idx = html.find(BEGIN_SENTINEL)
    end_idx = html.find(END_SENTINEL)

    if begin_idx == -1 or end_idx == -1:
        print("ERROR: sentinel comments not found in index.html", file=sys.stderr)
        sys.exit(1)

    after_begin = html.index("\n", begin_idx) + 1
    replacement = js_block + "\n"
    new_html = html[:after_begin] + replacement + html[end_idx:]

    INDEX_HTML.write_text(new_html)
    print(f"Written APPROVED_MAB_STATS to {INDEX_HTML}", file=sys.stderr)


if __name__ == "__main__":
    main()
