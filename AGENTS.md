# LAMBS — guide for AI agents

This file serves **two audiences**. Jump to the section that matches the task.

| If you are… | Read |
|-------------|------|
| **Using LAMBS** to analyze antibody sequences (VH/VL, basket, liabilities, germline) | [Using LAMBS (analysis agents)](#using-lambs-analysis-agents) |
| **Changing LAMBS** (bugs, features, tests, germline DB, UI in `index.html`) | [Developing LAMBS (contributor agents)](#developing-lambs-contributor-agents) |

**Product rule (both audiences):** end users and analysis agents rely on **one** offline [`index.html`](index.html). Do not introduce a second runtime artifact, npm install for users, or a parallel analysis implementation.

---

## Using LAMBS (analysis agents)

Your job is to run analysis and consume **structured results**, not to edit the repo or reimplement algorithms.

### How to run analysis

1. Open [`index.html`](index.html) in a browser (file URL or static server).
2. Prefer **`window.LAMBS`** over scraping the DOM.

| Method | Returns |
|--------|---------|
| `LAMBS.analyzeMabReportFromRaw(vhRaw, vlRaw)` | `{ ok, report? \| errors? }` — preferred |
| `LAMBS.analyzeSingle(vhRaw, vlRaw)` | Report object, or `{ ok: false, errors }` on validation failure |
| `LAMBS.analyzeMabReport(vhSeq, vlSeq)` | Report from **cleaned** AA strings |
| `LAMBS.buildChainReport(varSeq, 'VH' \| 'VL')` | Variable-region chain report |
| `LAMBS.buildClusterReport(basket, threshold, filterNames)` | Cluster report |
| `LAMBS.getHumanIgGEuNumbers(isotype)` | Eu index array for `IGHG1`–`IGHG4`, else `null` (IMGT Hu_IGHGnber) |
| `LAMBS.mapEuOntoConstantAlignment(aligned2, species, isotype)` | Per-column Eu numbers for a CH alignment, or `null` |
| `LAMBS.lambsReportToJSON(report)` | JSON string (no `_constMatch` / alignment blobs) |
| `LAMBS.lastReport` | Last single-tab report after **Analyze** |
| `LAMBS.lastClusterReport` | Last cluster report after **Cluster & Analyze** |

Playwright / browser MCP:

```javascript
await page.goto('file:///path/to/index.html');
const result = await page.evaluate((vh, vl) => LAMBS.analyzeMabReportFromRaw(vh, vl), vhRaw, vlRaw);
```

Human-visible shortcuts (same file): **Copy JSON** after **Analyze** (`#lambs-last-report`) or after **Cluster & Analyze** (`#lambs-last-cluster-report`), each with `dataset.json`.

### Input contract

- **FASTA or raw amino acid** — `LAMBS.cleanSequence` / `LAMBS.validateSequence`.
- **Basket CSV** — handled in the UI; column logic is `csvChainColumnIndices` inside `index.html` if you need to mirror import rules.
- **Do not** parse liability highlight colors from HTML when `report.liabilities` or `LAMBS` report builders are available.

### JSON report schema (`reportVersion: 1`)

#### Single mAb (`kind: "single-mab"`)

```json
{
  "reportVersion": 1,
  "kind": "single-mab",
  "chains": {
    "vh": {
      "chainType": "VH",
      "fullSequence": "...",
      "variable": {
        "chainType": "VH",
        "sequence": "...",
        "properties": { "length", "mw", "pI", "charge74", "charge55", "gravy", "percentiles?", "quartiles?" },
        "germline": {
          "vGene": { "gene", "species", "identity", "matches", "total", "cdr1?", "cdr2?" },
          "jGene": { "gene", "species", "identity", "matches", "total" } | null,
          "cdr3": { "start", "end", "length", "sequence" } | null
        },
        "liabilities": [{ "type": "Met Ox", "positions": [0, 1] }]
      },
      "constantRegion": { "isotype", "species", "regionKind", "identity", "matches", "total", "splitPos" } | null
    },
    "vl": { "...": "same shape" }
  }
}
```

#### Cluster (`kind: "cluster"`)

```json
{
  "reportVersion": 1,
  "kind": "cluster",
  "thresholdPercent": 70,
  "liabilityFilters": [],
  "inputCount": 2,
  "filteredCount": 2,
  "clusterCount": 1,
  "clusters": [{
    "index": 1,
    "size": 2,
    "representativeName": "A",
 "clonotype": {
  "chainType": "VH",
  "vGene": { "gene", "species" },
  "jGene": { "gene", "species" },
  "cdr3": { "length", "sequence" }
 } | null,
    "members": [{ "name", "vh": { "...panel..." }, "vl": { "...panel..." } }],
    "vh": {
      "chainType": "VH",
      "memberNames": [],
      "variableSequences": [],
      "consensus": "...",
      "consensusAnalysis": { "...buildChainReport..." },
      "conservation": { "perColumnIdentity": [] } | null
    },
    "vl": { "...": "same" },
    "constantRegions": [{ "name", "ch", "cl", "vhJGene", "vlJGene" }] | null
  }]
}
```

Exports omit alignment strings (`aligned1` / `aligned2`). `lambsReportToJSON` also strips internal `_constMatch`.

### What analysis agents should avoid

- Running **`node test.js`** expecting JSON output — that command is the **test suite** (see contributor section).
- Reading all of `index.html` into context (huge IMGT blocks) when `LAMBS.*` or grep for `window.LAMBS` / report schema suffices.
- Duplicating `detectLiabilities`, germline search, or clustering in new code instead of calling the shipped logic.

---

## Developing LAMBS (contributor agents)

Your job is to change **`index.html`** (and repo tooling), keep the single-file product, and keep tests green.

### Verify changes

```bash
node test.js
```

[`test.js`](test.js) extracts the first `<script>` from `index.html`, cuts at `// TAB SWITCHING`, and exports pure functions for Node. That mirrors how CI validates logic; it is **not** a supported analysis API for end users or analysis agents.

### Editing `index.html` efficiently

- **Do not** split germline databases into separate URLs for normal releases.
- **Do not** load the full file when exploring: grep section comments (`// SEQUENCE LIABILITY DETECTION`, `// MAIN ANALYSIS`, `// STRUCTURED REPORTS`, `window.LAMBS`).
- UI and JSON reports should share one code path (`buildChainReport` → render helpers) so agents and humans stay aligned.
- After UI changes, regenerate README screenshots per [`README.md`](README.md) Developer notes.

### Exported test surface

`test.js` re-exports many functions (see its export list ~lines 52–100): `cleanSequence`, `validateSequence`, `detectLiabilities`, `findClosestGermline`, `clusterSequences`, `buildChainReport`, `analyzeMabReport`, `buildClusterReport`, etc. Use these in **tests**, not as a second product entry point.

### When adding agent-facing behavior

- Add or extend APIs on **`window.LAMBS`** inside `index.html`.
- Document JSON field changes here under **Using LAMBS** (`reportVersion` bumps if breaking).
- Prefer **Copy JSON** / `lambsReportToJSON` over new DOM scraping hooks.
