#!/usr/bin/env node
//
// Extracts the <script> block from index.html and runs unit tests
// against the pure-logic functions (no DOM required).
//
// Usage: node test.js
//

const fs = require("fs");
const path = require("path");

// --- Extract pure-logic JS from index.html (everything before DOM code) ---
const html = fs.readFileSync(path.join(__dirname, "index.html"), "utf-8");
const scriptMatch = html.match(/<script>([\s\S]*?)<\/script>/);
if (!scriptMatch) {
  console.error("No <script> block found");
  process.exit(1);
}

// Cut at the MAIN ANALYSIS section — everything before it is pure logic
const CUTOFF =
  "// ================================================================\n// MAIN ANALYSIS";
const pureJS = scriptMatch[1].split(CUTOFF)[0];

const wrappedJS = `
(function(exports) {
  ${pureJS}
  exports.cleanSequence = cleanSequence;
  exports.buildRuler = buildRuler;
  exports.detectLiabilities = detectLiabilities;
  exports.alignVGene = alignVGene;
  exports.findClosestGermline = findClosestGermline;
  exports.AA_PROPERTY = AA_PROPERTY;
  exports.V_GENE_DB = V_GENE_DB;
})(this);
`;
const ctx = {};
new Function(wrappedJS).call(ctx);

const {
  cleanSequence,
  buildRuler,
  detectLiabilities,
  alignVGene,
  findClosestGermline,
  AA_PROPERTY,
  V_GENE_DB,
} = ctx;

// --- Test harness ---
let passed = 0,
  failed = 0;
function assert(condition, msg) {
  if (condition) {
    passed++;
  } else {
    failed++;
    console.error(`  FAIL: ${msg}`);
  }
}
function section(name) {
  console.log(`\n--- ${name} ---`);
}

// ============================================================
// cleanSequence
// ============================================================
section("cleanSequence");

assert(cleanSequence("evqlves") === "EVQLVES", "converts to uppercase");
assert(cleanSequence("EVQ LVES") === "EVQLVES", "strips spaces");
assert(cleanSequence("EVQ123LVES") === "EVQLVES", "strips digits");
assert(
  cleanSequence(">header\nEVQL\nVES") === "EVQLVES",
  "strips FASTA header",
);
assert(cleanSequence("  \n  ") === "", "empty after stripping");

// ============================================================
// buildRuler
// ============================================================
section("buildRuler");

const r25 = buildRuler(25).join("");
assert(r25[0] === "1", "position 1 at column 0");
assert(r25[8] === "1" && r25[9] === "0", "position 10 at columns 8-9");
assert(r25[18] === "2" && r25[19] === "0", "position 20 at columns 18-19");
assert(r25.length === 25, "ruler length matches input");

const r5 = buildRuler(5).join("");
assert(r5[0] === "1", "short ruler starts with 1");
assert(r5.length === 5, "short ruler correct length");

const r110 = buildRuler(110).join("");
assert(
  r110[97] === "1" && r110[98] === "0" && r110[99] === "0",
  "position 100 right-aligned to column 99",
);

// ============================================================
// detectLiabilities
// ============================================================
section("detectLiabilities");

// Unpaired cysteine (odd count)
assert(
  detectLiabilities("ACDA")["Unp. Cys"]?.length === 1,
  "single C = unpaired",
);
assert(!detectLiabilities("ACCA")["Unp. Cys"], "two C = paired, no flag");

// Methionine oxidation
assert(detectLiabilities("AMGA")["Met Ox"]?.includes(1), "M detected at pos 1");
assert(!detectLiabilities("AGGA")["Met Ox"], "no M = no met ox");

// N-linked glycosylation: N-X-S/T where X!=P
assert(
  detectLiabilities("ANST")["Glycosyl"]?.length === 3,
  "NAS motif flags 3 positions",
);
assert(
  !detectLiabilities("ANPA")["Glycosyl"],
  "NPS should not flag (but NPA has X=N, check...)",
);
// More specific: NPT should NOT flag because X=P
assert(!detectLiabilities("NPTA")["Glycosyl"], "NPT not flagged (X=P)");
// NAS should flag
const glyc = detectLiabilities("XNASY");
assert(glyc["Glycosyl"]?.includes(1), "glycosyl flags N");
assert(glyc["Glycosyl"]?.includes(2), "glycosyl flags A (middle)");
assert(glyc["Glycosyl"]?.includes(3), "glycosyl flags S (third)");

// Deamidation: NG, NS, NT, NN, NH — both residues flagged
const deamNG = detectLiabilities("ANGA")["Deamid"];
assert(
  deamNG?.includes(1) && deamNG?.includes(2),
  "NG deamidation flags both N and G",
);
const deamNS = detectLiabilities("ANSA")["Deamid"];
assert(
  deamNS?.includes(1) && deamNS?.includes(2),
  "NS deamidation flags both N and S",
);
const deamNT = detectLiabilities("ANTA")["Deamid"];
assert(
  deamNT?.includes(1) && deamNT?.includes(2),
  "NT deamidation flags both N and T",
);
const deamNN = detectLiabilities("ANNA")["Deamid"];
assert(
  deamNN?.includes(1) && deamNN?.includes(2),
  "NN deamidation flags both Ns",
);
const deamNH = detectLiabilities("ANHA")["Deamid"];
assert(
  deamNH?.includes(1) && deamNH?.includes(2),
  "NH deamidation flags both N and H",
);
assert(!detectLiabilities("ANDA")["Deamid"], "ND is NOT deamidation");

// Isomerization: DG, DS — both residues flagged
const isoDG = detectLiabilities("ADGA")["Isomer"];
assert(
  isoDG?.includes(1) && isoDG?.includes(2),
  "DG isomerization flags both D and G",
);
const isoDS = detectLiabilities("ADSA")["Isomer"];
assert(
  isoDS?.includes(1) && isoDS?.includes(2),
  "DS isomerization flags both D and S",
);
assert(!detectLiabilities("ADPA")["Isomer"], "DP is NOT isomerization");

// Acid-labile cleavage: DP — both residues flagged
const cleavDP = detectLiabilities("ADPA")["Cleavage"];
assert(
  cleavDP?.includes(1) && cleavDP?.includes(2),
  "DP cleavage flags both D and P",
);
assert(!detectLiabilities("ADGA")["Cleavage"], "DG is NOT cleavage");

// N-terminal Q
assert(
  detectLiabilities("QAAA")["N-term Q"]?.[0] === 0,
  "leading Q flagged at index 0",
);
assert(!detectLiabilities("AQAA")["N-term Q"], "Q not at N-term = no flag");
assert(!detectLiabilities("")["N-term Q"], "empty sequence = no N-term Q");

// Tryptophan oxidation
// assert(detectLiabilities("AWGA")["Trp Ox"]?.includes(1), "W detected");
// assert(!detectLiabilities("AGGA")["Trp Ox"], "no W = no trp ox");

// ============================================================
// Full VH/VL liability check
// ============================================================
section("detectLiabilities (full sequences)");

const VH =
  "EVQLVESGGGLVQPGGSLRLSCAASGFTFSDYWMHWVRQAPGKGLEWVSRINSDGSSTSYADSVKGRFTISRDNAKNTLYLQMNSLRAEDTAVYYCTTDGWGFDYWGQGTLVTVSS";
const vhL = detectLiabilities(VH);
assert(!vhL["Unp. Cys"], "VH has even cysteine count (2)");
// assert(vhL["Met Ox"]?.length === 2, "VH has 2 methionines");
// assert(vhL["Trp Ox"]?.length === 5, "VH has 5 tryptophans");
assert(vhL["Deamid"]?.length > 0, "VH has deamidation sites");
assert(vhL["Isomer"]?.length > 0, "VH has isomerization sites");
assert(!vhL["Glycosyl"], "VH has no N-linked glycosylation");

const VL =
  "DIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSYSTPLTFGGGTKVEIK";
const vlL = detectLiabilities(VL);
assert(vlL["Met Ox"]?.length === 1, "VL has 1 methionine");
// assert(vlL["Trp Ox"]?.length === 1, "VL has 1 tryptophan");

// ============================================================
// alignVGene (semi-global: germline fully consumed, query free suffix)
// ============================================================
section("alignVGene");

// Identical sequences — same as global alignment
const sg1 = alignVGene("ABCDEF", "ABCDEF");
assert(sg1.identity === 1.0, "identical sequences = 100% identity");
assert(sg1.aligned1 === "ABCDEF", "identical align1");
assert(sg1.aligned2 === "ABCDEF", "identical align2");
assert(sg1.matches === 6 && sg1.total === 6, "identical counts");

// Single mismatch
const sg2 = alignVGene("ABCDEF", "ABXDEF");
assert(
  sg2.aligned1.replace(/-/g, "") === "ABCDEF",
  "mismatch: query chars preserved",
);
assert(
  sg2.aligned2.replace(/-/g, "") === "ABXDEF",
  "mismatch: germline chars preserved",
);
assert(sg2.matches === 5, "mismatch: 5 matches");

// Query longer than germline (D/J suffix should be free)
const sg3 = alignVGene("ABCDEFGHIJ", "ABCDEF");
assert(sg3.identity === 1.0, "prefix match = 100% V-gene identity");
assert(
  sg3.aligned1.replace(/-/g, "") === "ABCDEFGHIJ",
  "long query: all query chars present",
);
const germChars = sg3.aligned2.replace(/-/g, "");
assert(germChars === "ABCDEF", "long query: all germline chars present");
const trailingGaps = sg3.aligned2.slice(-4);
assert(
  trailingGaps === "----",
  "long query: trailing gaps in germline for D/J suffix",
);

// Germline with insertion relative to query
const sg4 = alignVGene("ABDEF", "ABCDEF");
assert(
  sg4.aligned1.replace(/-/g, "") === "ABDEF",
  "insertion: query chars preserved",
);
assert(
  sg4.aligned2.replace(/-/g, "") === "ABCDEF",
  "insertion: germline chars preserved",
);
const gapCount = (sg4.aligned1.match(/-/g) || []).length;
assert(gapCount >= 1, "insertion: at least 1 gap in query");

// Empty sequences
const sg5 = alignVGene("", "");
assert(sg5.identity === 0, "empty sequences = 0 identity");
assert(sg5.aligned1 === "" && sg5.aligned2 === "", "empty aligned strings");

// One empty (query only, no germline)
const sg6 = alignVGene("ABC", "");
assert(sg6.aligned2 === "---", "empty germline: gaps in germline row");

// ============================================================
// findClosestGermline (with a test sequence injected)
// ============================================================
section("findClosestGermline (with test data)");

// Temporarily inject a test germline
V_GENE_DB.human.IGHV["TEST-GENE*01"] = VH.substring(0, 98);
const germTest = findClosestGermline(VH, "VH");
assert(germTest !== null, "finds match with test data");
assert(germTest.gene === "TEST-GENE*01", "correct gene name");
assert(germTest.species === "human", "correct species");
assert(germTest.identity > 0.8, "high identity for near-match");
delete V_GENE_DB.human.IGHV["TEST-GENE*01"];

// ============================================================
// AA_PROPERTY map
// ============================================================
section("AA_PROPERTY");

assert(AA_PROPERTY["A"] === "hydrophobic", "A is hydrophobic");
assert(
  AA_PROPERTY["F"] === "aromatic",
  "F is aromatic (overrides hydrophobic)",
);
assert(AA_PROPERTY["W"] === "aromatic", "W is aromatic");
assert(AA_PROPERTY["S"] === "polar", "S is polar");
assert(AA_PROPERTY["K"] === "positive", "K is positive");
assert(AA_PROPERTY["D"] === "negative", "D is negative");
assert(AA_PROPERTY["C"] === "special", "C is special");
assert(AA_PROPERTY["G"] === "special", "G is special");
assert(AA_PROPERTY["P"] === "special", "P is special");

// ============================================================
// Summary
// ============================================================
console.log(`\n========================================`);
console.log(`  ${passed} passed, ${failed} failed`);
console.log(`========================================`);
process.exit(failed > 0 ? 1 : 0);
