#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# ///
"""Parse germline FASTA & IgBLAST TSV files and inject V_GENE_DB and J_GENE_DB into index.html.

Usage (from repository root):
    python scripts/populate_germlines.py
"""

import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
GERMLINES_DIR = ROOT / "germlines"
HTML_FILE = ROOT / "index.html"

SPECIES = ["human", "mouse"]
V_GENES = ["ighv", "igkv", "iglv"]
J_GENES = ["ighj", "igkj", "iglj"]


def parse_fasta(path: Path) -> dict[str, str]:
    """Return {gene_name: aa_sequence} from a FASTA file, skipping partial entries."""
    entries: dict[str, str] = {}
    current_name: str | None = None
    current_seq: list[str] = []

    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if current_name is not None:
                entries[current_name] = "".join(current_seq)
            if "partial" in line.lower():
                current_name = None
                current_seq = []
                continue
            fields = line[1:].split("|")
            current_name = fields[1] if len(fields) >= 2 else None
            current_seq = []
        elif current_name is not None:
            current_seq.append(line.strip())

    if current_name is not None:
        entries[current_name] = "".join(current_seq)

    return entries


def parse_j_fasta(path: Path) -> dict[str, str]:
    """Return {gene_name: aa_sequence} from a J gene FASTA (amino acid sequences).

    Only functional entries (F or (F)) are included; ORFs and pseudogenes are skipped.
    Entries marked as partial in the header are also skipped.
    """
    entries: dict[str, str] = {}
    current_name: str | None = None
    current_seq: list[str] = []

    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if current_name is not None:
                entries[current_name] = "".join(current_seq)

            current_name = None
            current_seq = []

            if "partial" in line.lower():
                continue

            fields = line[1:].split("|")
            if len(fields) < 2:
                continue

            # Field 3 (0-based): functionality — include only F and (F)
            func = fields[3].strip() if len(fields) >= 4 else ""
            if func not in ("F", "(F)"):
                continue

            current_name = fields[1]

        elif current_name is not None:
            current_seq.append(line.strip())

    if current_name is not None:
        entries[current_name] = "".join(current_seq)

    return entries


def parse_igblast_tsv(path: Path) -> dict[str, dict]:
    """Extract CDR1/CDR2 amino acid positions from IgBLAST AIRR TSV output.

    Returns {gene_name: {"cdr1": [aa_start, aa_end], "cdr2": [aa_start, aa_end]}}
    where positions are 0-based into the amino acid sequence, end is exclusive.
    """
    cdr_data: dict[str, dict] = {}
    lines = path.read_text().splitlines()
    if not lines:
        return cdr_data

    header = lines[0].split("\t")
    col = {name: i for i, name in enumerate(header)}

    for line in lines[1:]:
        fields = line.split("\t")
        gene = fields[col["sequence_id"]]

        try:
            c1s = int(fields[col["cdr1_start"]])
            c1e = int(fields[col["cdr1_end"]])
            c2s = int(fields[col["cdr2_start"]])
            c2e = int(fields[col["cdr2_end"]])
        except (ValueError, KeyError, IndexError):
            continue

        cdr_data[gene] = {
            "cdr1": [(c1s - 1) // 3, (c1e - 1) // 3 + 1],
            "cdr2": [(c2s - 1) // 3, (c2e - 1) // 3 + 1],
        }

    return cdr_data


def build_v_js_object(db: dict, cdr_db: dict) -> str:
    """Build the JS source for the V_GENE_DB constant."""
    lines = ["const V_GENE_DB = {"]
    for species in SPECIES:
        lines.append(f"  {species}: {{")
        for gene in V_GENES:
            gene_upper = gene.upper()
            entries = db.get(species, {}).get(gene_upper, {})
            lines.append(f"    {gene_upper}: {{")
            for name, seq in entries.items():
                cdr = cdr_db.get(species, {}).get(name)
                if cdr:
                    obj = {"s": seq, "cdr1": cdr["cdr1"], "cdr2": cdr["cdr2"]}
                else:
                    obj = {"s": seq}
                lines.append(f"      {json.dumps(name)}: {json.dumps(obj)},")
            lines.append("    },")
        lines.append("  },")
    lines.append("};")
    return "\n".join(lines)


def build_j_js_object(db: dict) -> str:
    """Build the JS source for the J_GENE_DB constant."""
    lines = ["const J_GENE_DB = {"]
    for species in SPECIES:
        lines.append(f"  {species}: {{")
        for gene in J_GENES:
            gene_upper = gene.upper()
            entries = db.get(species, {}).get(gene_upper, {})
            lines.append(f"    {gene_upper}: {{")
            for name, seq in entries.items():
                lines.append(f"      {json.dumps(name)}: {json.dumps({'s': seq})},")
            lines.append("    },")
        lines.append("  },")
    lines.append("};")
    return "\n".join(lines)


def main() -> None:
    # ── V genes ──────────────────────────────────────────────────────
    v_db: dict[str, dict[str, dict[str, str]]] = {}
    v_total = 0
    print("V genes:")
    for species in SPECIES:
        v_db[species] = {}
        for gene in V_GENES:
            fasta_path = GERMLINES_DIR / f"{species}-{gene}.fasta"
            entries = parse_fasta(fasta_path)
            v_db[species][gene.upper()] = entries
            v_total += len(entries)
            print(f"  {fasta_path.name}: {len(entries)} sequences")

    print(f"  Total: {v_total} sequences (partials excluded)")

    cdr_db: dict[str, dict[str, dict]] = {}
    for species in SPECIES:
        tsv_path = GERMLINES_DIR / f"igblast_{species}_V.tsv"
        if tsv_path.exists():
            cdr_db[species] = parse_igblast_tsv(tsv_path)
            print(f"  {tsv_path.name}: CDR annotations for {len(cdr_db[species])} genes")
        else:
            cdr_db[species] = {}

    # ── J genes ──────────────────────────────────────────────────────
    j_db: dict[str, dict[str, dict[str, str]]] = {}
    j_total = 0
    print("J genes:")
    for species in SPECIES:
        j_db[species] = {}
        for gene in J_GENES:
            fasta_path = GERMLINES_DIR / f"{species}-{gene}.fasta"
            entries = parse_j_fasta(fasta_path)
            j_db[species][gene.upper()] = entries
            j_total += len(entries)
            print(f"  {fasta_path.name}: {len(entries)} sequences")

    print(f"  Total: {j_total} sequences (ORFs/pseudogenes/partials excluded)")

    # ── Inject into index.html ────────────────────────────────────────
    html = HTML_FILE.read_text()

    v_pattern = re.compile(r"const V_GENE_DB = \{.*?\n\};", re.DOTALL)
    if not v_pattern.search(html):
        raise SystemExit("Could not locate V_GENE_DB block in index.html")
    html = v_pattern.sub(build_v_js_object(v_db, cdr_db), html)

    j_pattern = re.compile(r"const J_GENE_DB = \{.*?\n\};", re.DOTALL)
    if not j_pattern.search(html):
        raise SystemExit("Could not locate J_GENE_DB block in index.html")
    html = j_pattern.sub(build_j_js_object(j_db), html)

    HTML_FILE.write_text(html)
    print(f"  Updated {HTML_FILE}")


if __name__ == "__main__":
    main()
