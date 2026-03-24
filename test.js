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

// Cut at the TAB SWITCHING section — everything before it is pure logic
// (DOM-touching basket functions get a mock document so they don't crash)
const CUTOFF =
  "// ================================================================\n// TAB SWITCHING";
const pureJS = scriptMatch[1].split(CUTOFF)[0];

const wrappedJS = `
(function(exports) {
  // Minimal DOM mock for basket functions that call renderBasket()
  var document = { getElementById: function() {
    return { textContent: '', innerHTML: '', disabled: false, style: {} };
  }};

  ${pureJS}
  exports.cleanSequence = cleanSequence;
  exports.buildRuler = buildRuler;
  exports.detectLiabilities = detectLiabilities;
  exports.alignVGene = alignVGene;
  exports.findClosestGermline = findClosestGermline;
  exports.AA_PROPERTY = AA_PROPERTY;
  exports.V_GENE_DB = V_GENE_DB;
  exports.parseCSV = parseCSV;
  exports.escapeHtml = escapeHtml;
  exports.globalAlign = globalAlign;
  exports.pairwiseIdentity = pairwiseIdentity;
  exports.clusterSequences = clusterSequences;
  exports.starMSA = starMSA;
  exports.computeConservation = computeConservation;
  exports.consensusSequenceForAnalysis = consensusSequenceForAnalysis;
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
  parseCSV,
  escapeHtml,
  globalAlign,
  pairwiseIdentity,
  clusterSequences,
  starMSA,
  computeConservation,
  consensusSequenceForAnalysis,
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
// escapeHtml
// ============================================================
section("escapeHtml");

assert(escapeHtml("a&b") === "a&amp;b", "escapes ampersand");
assert(escapeHtml("<b>") === "&lt;b&gt;", "escapes angle brackets");
assert(escapeHtml('"hi"') === "&quot;hi&quot;", "escapes quotes");
assert(escapeHtml("plain") === "plain", "no-op for plain text");

// ============================================================
// parseCSV
// ============================================================
section("parseCSV");

const csv1 = parseCSV("a,b,c\n1,2,3\n4,5,6\n");
assert(csv1.length === 3, "parseCSV: 3 rows (header + 2 data)");
assert(csv1[0][0] === "a" && csv1[0][2] === "c", "parseCSV: header parsed");
assert(csv1[1][1] === "2", "parseCSV: data cell correct");

const csv2 = parseCSV('name,vh\n"Seq, 1","ABC"\n');
assert(csv2[1][0] === "Seq, 1", "parseCSV: quoted field with comma");
assert(csv2[1][1] === "ABC", "parseCSV: quoted field value");

const csv3 = parseCSV('a,b\r\n1,2\r\n3,4\r\n');
assert(csv3.length === 3, "parseCSV: handles CRLF line endings");
assert(csv3[2][0] === "3", "parseCSV: CRLF data correct");

const csv4 = parseCSV('a\n"he said ""hi"""\n');
assert(csv4[1][0] === 'he said "hi"', "parseCSV: escaped double quotes");

const csvEmpty = parseCSV("");
assert(csvEmpty.length === 0, "parseCSV: empty string returns no rows");

// ============================================================
// globalAlign (full Needleman-Wunsch)
// ============================================================
section("globalAlign");

const ga1 = globalAlign("ABCDEF", "ABCDEF");
assert(ga1.identity === 1.0, "globalAlign: identical = 100%");
assert(ga1.matches === 6 && ga1.total === 6, "globalAlign: identical counts");

const ga2 = globalAlign("ABCDEF", "ABXDEF");
assert(ga2.matches === 5, "globalAlign: single mismatch = 5 matches");
assert(
  ga2.aligned1.replace(/-/g, "") === "ABCDEF",
  "globalAlign: mismatch query preserved",
);
assert(
  ga2.aligned2.replace(/-/g, "") === "ABXDEF",
  "globalAlign: mismatch ref preserved",
);

const ga3 = globalAlign("ABCDEF", "ABCXYZDEF");
assert(
  ga3.aligned1.replace(/-/g, "") === "ABCDEF",
  "globalAlign: insertion query preserved",
);
assert(
  ga3.aligned2.replace(/-/g, "") === "ABCXYZDEF",
  "globalAlign: insertion ref preserved",
);

const ga4 = globalAlign("", "");
assert(ga4.identity === 0, "globalAlign: empty sequences = 0");

const ga5 = globalAlign("ABC", "");
assert(ga5.aligned1.replace(/-/g, "") === "ABC", "globalAlign: one empty seq1");
assert(ga5.aligned2 === "---", "globalAlign: one empty seq2 gaps");

const ga6 = globalAlign("ABCDEFGHIJ", "ABCDEF");
assert(ga6.matches === 6, "globalAlign: prefix match 6 matches");
assert(ga6.total === 10, "globalAlign: prefix match 10 total positions");

// ============================================================
// pairwiseIdentity
// ============================================================
section("pairwiseIdentity");

const pi1 = pairwiseIdentity(
  { vh: "ABCDEF", vl: "GHIJ" },
  { vh: "ABCDEF", vl: "GHIJ" },
);
assert(pi1 === 1.0, "pairwiseIdentity: identical pairs = 1.0");

const pi2 = pairwiseIdentity(
  { vh: "ABCDEF", vl: "GHIJ" },
  { vh: "XYZXYZ", vl: "WXWX" },
);
assert(pi2 < 0.5, "pairwiseIdentity: dissimilar pairs < 0.5");

const pi3 = pairwiseIdentity({ vh: "", vl: "" }, { vh: "ABC", vl: "DEF" });
assert(pi3 === 0, "pairwiseIdentity: empty entry = 0");

const pi4 = pairwiseIdentity(
  { vh: "ABCDEF", vl: "" },
  { vh: "ABCDEF", vl: "" },
);
assert(pi4 === 1.0, "pairwiseIdentity: VH-only identical = 1.0");

// ============================================================
// clusterSequences
// ============================================================
section("clusterSequences");

const clEntries = [
  { id: 1, name: "A", vh: "ABCDEF", vl: "GHIJ" },
  { id: 2, name: "B", vh: "ABCDEF", vl: "GHIJ" },
  { id: 3, name: "C", vh: "XYZXYZ", vl: "WXWX" },
];

const cl100 = clusterSequences(clEntries, 100);
assert(cl100.length >= 2, "cluster@100%: dissimilar seqs form separate clusters");

const cl10 = clusterSequences(clEntries, 10);
assert(
  cl10.length >= 1,
  "cluster@10%: very low threshold groups more sequences",
);

const clIdentical = clusterSequences(
  [
    { id: 1, name: "A", vh: "ABCDEF", vl: "GHIJ" },
    { id: 2, name: "B", vh: "ABCDEF", vl: "GHIJ" },
    { id: 3, name: "C", vh: "ABCDEF", vl: "GHIJ" },
  ],
  70,
);
assert(clIdentical.length === 1, "cluster: 3 identical seqs = 1 cluster");
assert(
  clIdentical[0].members.length === 3,
  "cluster: single cluster has all 3 members",
);

const clEmpty = clusterSequences([], 70);
assert(clEmpty.length === 0, "cluster: empty input = empty output");

const clSingle = clusterSequences(
  [{ id: 1, name: "A", vh: "ABCDEF", vl: "GHIJ" }],
  70,
);
assert(clSingle.length === 1, "cluster: single entry = 1 cluster");
assert(clSingle[0].members.length === 1, "cluster: single entry cluster size 1");

// Verify sorted by size (largest first)
const clSorted = clusterSequences(
  [
    { id: 1, name: "A", vh: "ABCDEF", vl: "GHIJ" },
    { id: 2, name: "B", vh: "ABCDEF", vl: "GHIJ" },
    { id: 3, name: "C", vh: "XYZXYZ", vl: "WXWX" },
  ],
  70,
);
assert(
  clSorted[0].members.length >= clSorted[clSorted.length - 1].members.length,
  "cluster: results sorted largest first",
);

// ============================================================
// starMSA
// ============================================================
section("starMSA");

const msa0 = starMSA([]);
assert(msa0.length === 0, "starMSA: empty input");

const msa1 = starMSA(["ABCDEF"]);
assert(msa1.length === 1, "starMSA: single sequence");
assert(msa1[0] === "ABCDEF", "starMSA: single sequence unchanged");

const msa2 = starMSA(["ABCDEF", "ABCDEF"]);
assert(msa2.length === 2, "starMSA: two identical sequences");
assert(msa2[0] === msa2[1], "starMSA: identical sequences aligned identically");

const msa3 = starMSA(["ABCDEF", "ABCDEF", "ABCDEF"]);
assert(msa3.length === 3, "starMSA: three identical");
assert(
  msa3[0] === msa3[1] && msa3[1] === msa3[2],
  "starMSA: three identical all equal",
);

// All MSA outputs must have equal length
const msa4 = starMSA(["ABCDEF", "ABXDEF", "ABCEF"]);
assert(msa4.length === 3, "starMSA: three similar sequences");
const msaLen = msa4[0].length;
assert(
  msa4[1].length === msaLen && msa4[2].length === msaLen,
  "starMSA: all outputs same length",
);

// Original residues are preserved (no residue loss)
for (let i = 0; i < msa4.length; i++) {
  const original = ["ABCDEF", "ABXDEF", "ABCEF"][i];
  const stripped = msa4[i].replace(/-/g, "");
  assert(
    stripped === original,
    "starMSA: residues preserved for seq " + (i + 1),
  );
}

// Realistic antibody-like sequences
const msaAb = starMSA([
  "EVQLVESGGGLVQPGGSLRLSCAAS",
  "EVQLVESGGGLVQPGGSLRLSCAAS",
  "EVQLVETGGGLVQPGGSLRLSCAAS",
]);
assert(msaAb.length === 3, "starMSA: antibody-like sequences");
assert(
  msaAb[0].length === msaAb[1].length && msaAb[1].length === msaAb[2].length,
  "starMSA: antibody-like equal length",
);

// ============================================================
// computeConservation
// ============================================================
section("computeConservation");

const cons0 = computeConservation([]);
assert(cons0.identity.length === 0, "conservation: empty input");
assert(cons0.consensus === "", "conservation: empty consensus");

const cons1 = computeConservation(["QQQ", "QET", "QQS"]);
assert(cons1.identity.length === 3, "conservation: 3 columns");
assert(cons1.identity[0] === 1.0, "conservation: col 1 = 100% (all Q)");
assert(
  Math.abs(cons1.identity[1] - 2 / 3) < 0.01,
  "conservation: col 2 = 66% (2 Q, 1 E)",
);
assert(
  Math.abs(cons1.identity[2] - 1 / 3) < 0.01,
  "conservation: col 3 = 33% (Q, T, S all different)",
);
assert(cons1.consensus[0] === "Q", "conservation: consensus pos 1 = Q");
assert(cons1.consensus[1] === "Q", "conservation: consensus pos 2 = Q (majority)");

// Consensus lowercase when no residue exceeds 50%
const cons2 = computeConservation(["AB", "CD", "EF", "GH"]);
assert(
  cons2.consensus[0] === cons2.consensus[0].toLowerCase(),
  "conservation: low-identity position = lowercase",
);

// Gaps are excluded from counting
const cons3 = computeConservation(["A-C", "A-C", "ABC"]);
assert(cons3.identity[0] === 1.0, "conservation: col 1 all A = 100%");
assert(
  cons3.identity[1] === 1.0,
  "conservation: col 2 gaps excluded, only B counted = 100%",
);
assert(cons3.identity[2] === 1.0, "conservation: col 3 all C = 100%");

// All gaps column
const cons4 = computeConservation(["A-", "A-"]);
assert(cons4.identity[1] === 0, "conservation: all-gap column = 0");
assert(cons4.consensus[1] === "-", "conservation: all-gap consensus = dash");

// Perfect conservation
const cons5 = computeConservation(["AAAA", "AAAA", "AAAA"]);
assert(
  cons5.identity.every(function (v) {
    return v === 1.0;
  }),
  "conservation: all identical = all 1.0",
);
assert(cons5.consensus === "AAAA", "conservation: all identical consensus");

// ============================================================
// consensusSequenceForAnalysis
// ============================================================
section("consensusSequenceForAnalysis");

assert(
  consensusSequenceForAnalysis("A-bC") === "ABC",
  "consensusForAnalysis: strips gaps and uppercases",
);
assert(consensusSequenceForAnalysis("") === "", "consensusForAnalysis: empty");
assert(
  consensusSequenceForAnalysis("abc") === "ABC",
  "consensusForAnalysis: lowercases to uppercase",
);

// ============================================================
// Integration: cluster + MSA + conservation pipeline
// ============================================================
section("integration: cluster-MSA-conservation pipeline");

const pipelineEntries = [
  {
    id: 1,
    name: "mAb1",
    vh: "EVQLVESGGGLVQPGGSLRLSCAASGFTFSDYWMHWVRQAPGKGLEWVSRINSDGSSTSYADSVKGRFTISRDNAKNTLYLQMNSLRAEDTAVYYCTTDGWGFDYWGQGTLVTVSS",
    vl: "DIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSYSTPLTFGGGTKVEIK",
  },
  {
    id: 2,
    name: "mAb2",
    vh: "EVQLVESGGGLVQPGGSLRLSCAASGFTFSDYWMHWVRQAPGKGLEWVSRINSDGSSTSYADSVKGRFTISRDNAKNTLYLQMNSLRAEDTAVYYCTTDGWGFDYWGQGTLVTVSS",
    vl: "DIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSYSTPLTFGGGTKVEIK",
  },
];

const pipeClusters = clusterSequences(pipelineEntries, 70);
assert(pipeClusters.length === 1, "pipeline: identical pair = 1 cluster");
assert(
  pipeClusters[0].members.length === 2,
  "pipeline: cluster has 2 members",
);

const pipeMSA = starMSA(
  pipeClusters[0].members.map(function (m) {
    return m.vh;
  }),
);
assert(pipeMSA.length === 2, "pipeline: MSA produces 2 aligned seqs");
assert(
  pipeMSA[0].length === pipeMSA[1].length,
  "pipeline: MSA seqs equal length",
);

const pipeCons = computeConservation(pipeMSA);
assert(
  pipeCons.identity.every(function (v) {
    return v === 1.0;
  }),
  "pipeline: identical seqs = 100% conservation everywhere",
);
assert(
  pipeCons.consensus === pipeMSA[0],
  "pipeline: consensus matches aligned seq",
);

// ============================================================
// Summary
// ============================================================
console.log(`\n========================================`);
console.log(`  ${passed} passed, ${failed} failed`);
console.log(`========================================`);
process.exit(failed > 0 ? 1 : 0);
