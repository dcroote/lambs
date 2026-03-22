# Original prompt

```
I would like to build a analysis dashboard for monoclonal antibody variable region sequences.

# User Inputs
- One VH amino acid sequence and one VL amino acid sequence

# Requirements
- MANDATORY: EVERYTHING MUST BE IN THE SINGLE index.html file. @index.html
- ALL ANALYSIS SHOULD OCCUR USING EMBEDDED <script> tags.
- All text should be monospace font.
- There is an input box for VH sequence and input box for VL sequence

# Features

## Sequence Liability Highlighting:
-Unpaired cysteine
- Methionine oxidation
- N-linked glycosylation
- Deamidation
- Isomerization
- Other common sequence liabilities
It should be very clear visually where each of these occurs in the sequence

## Species identification:
- Identify species (limited to mouse and human) by closest V gene using sequence similarity. The closest v gene will be the "germline." This should be done for both the VH and VL separately against the appropriate gene set
- For now, create a placeholder for where the V gene sequences should go, and I'll input them from IMGT

## Sequence and germline alignment:
- Alignment of sequence to identified germine V gene. In between these aligned sequences should be a line that has an X for every mismatch and a space otherwise.

## Amino acid coloring
- option (toggle) to color amino acid text by property

```
