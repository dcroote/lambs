---
name: mAb Sequence Dashboard
overview: Build a single-file HTML dashboard that analyzes monoclonal antibody VH/VL amino acid sequences, highlighting sequence liabilities, identifying closest germline V gene (human/mouse), displaying sequence-germline alignment with mismatch markers, and offering toggleable amino acid property coloring.
todos:
  - id: html-structure
    content: "Build HTML skeleton: header, input section (VH/VL textareas, Analyze button, AA-coloring toggle), results container, legend"
    status: completed
  - id: css-styling
    content: "Embedded CSS: light theme, monospace font, annotation track styling, AA property colors, layout grid, responsive design"
    status: completed
  - id: liability-detection
    content: "JavaScript: liability scanner functions for all 7 motif types, return per-position boolean arrays for each liability"
    status: completed
  - id: liability-rendering
    content: "JavaScript: render sequence line plus one annotation track row per liability type beneath it, with position numbering"
    status: completed
  - id: vgene-placeholder
    content: "JavaScript: V_GENE_DB object with structured placeholders for human/mouse IGHV, IGKV, IGLV with clear comments for user to paste IMGT data"
    status: completed
  - id: nw-alignment
    content: "JavaScript: Needleman-Wunsch global alignment function (match +1, mismatch -1, gap -2), returns aligned pair and percent identity"
    status: completed
  - id: species-identification
    content: "JavaScript: iterate V gene DB, align input to each reference, select best match, report gene name + species + identity%"
    status: completed
  - id: alignment-display
    content: "JavaScript: render three-line alignment block (query / match-line / germline) with mismatch highlighting"
    status: completed
  - id: aa-coloring-toggle
    content: "JavaScript: toggle event listener that applies/removes AA property CSS classes across all rendered residues"
    status: completed
  - id: integration-testing
    content: Test end-to-end with a sample VH/VL pair to verify all panels render correctly
    status: completed
isProject: false
---

# mAb Variable Region Sequence Analysis Dashboard

Everything lives in `[index.html](index.html)` -- HTML structure, CSS styles, and all JavaScript logic in embedded `<script>` tags. No external dependencies.

## Layout / UI Structure

- **Header**: Dashboard title
- **Input Section**: Two `<textarea>` inputs (VH and VL), an "Analyze" button, and a toggle checkbox for amino acid property coloring
- **Results Section** (appears after analysis):
  - **Liability Analysis** panel -- one sub-panel per chain (VH, VL)
  - **Species / Germline Identification** panel -- closest V gene match for each chain
  - **Germline Alignment** panel -- three-line alignment display per chain
- **Legend**: Color key for liabilities and (when toggled) amino acid property groups
- All text uses `font-family: monospace`

## 1. Sequence Liability Analysis (Annotation Tracks)

The sequence is displayed as plain monospace text on the top line. Below it, each liability type gets its own **annotation track** -- a dedicated row of the same character width as the sequence. At positions where the liability applies, a marker character is printed; elsewhere the row is blank. This avoids any color overlap and keeps the sequence itself uncluttered.

### Liability Detection Rules

- **Unpaired Cysteine** -- odd total count of C in the chain (mark every C)
- **Methionine Oxidation** -- any M
- **N-linked Glycosylation** -- N-X-S/T where X != P (mark all residues in the motif)
- **Deamidation** -- NG, NS, NT, NN, NH (mark the N)
- **Isomerization** -- DG, DS, DT, DH (mark the D)
- **Acid-labile / Cleavage** -- DP motif (mark the D)
- **Tryptophan Oxidation** -- any W

### Display Format

Each track row is the same width as the sequence but mostly blank spaces. At liability positions, the **actual amino acid characters** are shown, plus **3 flanking residues on each side** for positional context. The liability residues themselves are rendered in **red text**; the flanking context residues are rendered in a muted gray. This makes it trivial to see exactly where in the sequence each liability falls, even on lower tracks.

Tracks with no hits are hidden entirely to save space.

```
Position: 1         11        21        31
Sequence: EVQLVESGGGLVQPGGSLRLSCAASGFTFSNYWMNWVR
Unp. Cys:                   RLSCAAS
Met Ox:                               YWMNWVR
Deamid:                           GFTFSNYWM
Trp Ox:                              FSNYWMNWVR
```

(In the actual rendering, the red characters are the flagged liability residues; the gray characters are the 3-residue flanking context. All characters are uppercase and maintain their exact column position matching the sequence above.)

Each track row has a short label on the left. A compact **legend** beneath explains the track labels and the red/gray convention.

## 2. Species Identification via V Gene Matching

### Germline Database Placeholder

A JavaScript object at the top of the `<script>` block will hold V gene reference sequences. Structure:

```javascript
const V_GENE_DB = {
  human: {
    IGHV: { "IGHV1-1*01": "PLACEHOLDER..." /* ~50+ entries */ },
    IGKV: { "IGKV1-5*01": "PLACEHOLDER..." /* entries */ },
    IGLV: { "IGLV1-40*01": "PLACEHOLDER..." /* entries */ },
  },
  mouse: {
    IGHV: { "IGHV1-1*01": "PLACEHOLDER..." /* entries */ },
    IGKV: { "IGKV1-5*01": "PLACEHOLDER..." /* entries */ },
    IGLV: { "IGLV1-40*01": "PLACEHOLDER..." /* entries */ },
  },
};
```

Each value is a placeholder string (e.g., `"EVQLVESGGGLVQPGG..."`) that the user will replace with real IMGT amino acid sequences. A clear comment block will mark where to paste data. I will seed 3-5 example entries per gene family with realistic-length placeholder strings so the structure is unambiguous.

### Matching Algorithm

For a given input sequence (VH or VL):

1. Iterate over all V gene entries in the appropriate families (IGHV for VH; IGKV + IGLV for VL) across both species.
2. Compute **percent identity** using a simple Needleman-Wunsch global alignment (BLOSUM-like scoring not needed -- simple match/mismatch scoring: +1 match, -1 mismatch, -2 gap).
3. Select the gene with the highest percent identity. Report: gene name, species, and percent identity.

### Needleman-Wunsch Implementation

A compact ~40-line NW function embedded in the script:

- Scoring: match = +1, mismatch = -1, gap = -2
- Returns: aligned sequence pair and identity percentage
- This same alignment is reused for the alignment display (Feature 3)

## 3. Sequence-to-Germline Alignment Display

For each chain (VH, VL), display a three-line block in monospace:

```
Query:    EVQLVESGGGLVQPGGSLRLSCAASGFTFS...
                   X
Germline: EVQLVESGGG-VQPGGSLRLSCAASGFTFS...
```

- Line 1: Input sequence (with gaps from alignment)
- Line 2: Match line -- a space for identity, `X` for mismatch, a space for gap positions (clean, only mismatches are marked)
- Line 3: Germline sequence (with gaps from alignment)

Mismatched positions on the query/germline lines will also be rendered with a subtle red text color to draw the eye.

## 4. Amino Acid Property Coloring (Toggle)

A checkbox toggle labeled "Color by amino acid property." When enabled, each residue `<span>` gets a text/background color based on its physicochemical group:

- **Hydrophobic** (A, V, I, L, M, F, W): dark blue text
- **Polar uncharged** (S, T, N, Q): green text
- **Positively charged** (K, R, H): red text
- **Negatively charged** (D, E): magenta text
- **Aromatic** (F, W, Y): gold/amber text (F, W overlap with hydrophobic -- aromatic takes precedence)
- **Special** (C, G, P): orange text

This coloring applies to both the liability display and the alignment display. When the toggle is off, text is default monospace black/white.

## 5. Styling and UX

- **Light-themed UI** (white/light gray background, dark text)
- Monospace font everywhere (`Courier New, monospace`)
- Responsive layout using CSS grid/flexbox
- Sequences displayed with position numbering every 10 residues
- Hover tooltips on annotation track markers showing the full liability name and motif
- Clear section dividers between VH and VL results
