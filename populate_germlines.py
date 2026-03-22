#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# ///
"""Parse germline FASTA & IgBLAST TSV files and inject V_GENE_DB into index.html."""

import json
import re
from pathlib import Path

GERMLINES_DIR = Path(__file__).parent / "germlines"
HTML_FILE = Path(__file__).parent / "index.html"

SPECIES = ["human", "mouse"]
GENES = ["ighv", "igkv", "iglv"]


def parse_fasta(path: Path) -> dict[str, str]:
    """Return {gene_name: sequence} from a FASTA file, skipping partial entries."""
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


def build_js_object(db: dict, cdr_db: dict) -> str:
    """Build the JS source for the V_GENE_DB constant."""
    lines = ["const V_GENE_DB = {"]
    for species in SPECIES:
        lines.append(f"  {species}: {{")
        for gene in GENES:
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


def main() -> None:
    db: dict[str, dict[str, dict[str, str]]] = {}
    total = 0
    for species in SPECIES:
        db[species] = {}
        for gene in GENES:
            fasta_path = GERMLINES_DIR / f"{species}-{gene}.fasta"
            entries = parse_fasta(fasta_path)
            db[species][gene.upper()] = entries
            total += len(entries)
            print(f"  {fasta_path.name}: {len(entries)} sequences")

    print(f"  Total: {total} sequences (partials excluded)")

    cdr_db: dict[str, dict[str, dict]] = {}
    for species in SPECIES:
        tsv_path = GERMLINES_DIR / f"igblast_{species}_V.tsv"
        if tsv_path.exists():
            cdr_db[species] = parse_igblast_tsv(tsv_path)
            print(f"  {tsv_path.name}: CDR annotations for {len(cdr_db[species])} genes")
        else:
            cdr_db[species] = {}

    html = HTML_FILE.read_text()

    pattern = re.compile(
        r"const V_GENE_DB = \{.*?\n\};",
        re.DOTALL,
    )
    if not pattern.search(html):
        raise SystemExit("Could not locate V_GENE_DB block in index.html")

    new_block = build_js_object(db, cdr_db)
    html = pattern.sub(new_block, html)
    HTML_FILE.write_text(html)
    print(f"  Updated {HTML_FILE}")


if __name__ == "__main__":
    main()
