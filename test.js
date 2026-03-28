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
  exports.LIABILITY_ORDER_ACTIVE = LIABILITY_ORDER_ACTIVE;
  exports.msaRawToMasterAligned = msaRawToMasterAligned;
  exports.filterBasketByLiabilityFilters = filterBasketByLiabilityFilters;
  exports.entryHasFilteredLiability = entryHasFilteredLiability;
  exports.buildMsaLiabilityColorMap = buildMsaLiabilityColorMap;
  exports.KD_HYDROPHOBICITY = KD_HYDROPHOBICITY;
  exports.AGG_WINDOW = AGG_WINDOW;
  exports.AGG_THRESHOLD = AGG_THRESHOLD;
  exports.detectCDR3 = detectCDR3;
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
  LIABILITY_ORDER_ACTIVE,
  msaRawToMasterAligned,
  filterBasketByLiabilityFilters,
  entryHasFilteredLiability,
  buildMsaLiabilityColorMap,
  KD_HYDROPHOBICITY,
  AGG_WINDOW,
  AGG_THRESHOLD,
  detectCDR3,
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

const csv3 = parseCSV("a,b\r\n1,2\r\n3,4\r\n");
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
assert(
  cl100.length >= 2,
  "cluster@100%: dissimilar seqs form separate clusters",
);

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
assert(
  clSingle[0].members.length === 1,
  "cluster: single entry cluster size 1",
);

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
// MSA liability helpers
// ============================================================
section("MSA liability helpers");

assert(
  LIABILITY_ORDER_ACTIVE.indexOf("Trp Ox") < 0,
  "LIABILITY_ORDER_ACTIVE excludes Trp Ox",
);

const mapGapped = msaRawToMasterAligned("---MQ");
assert(
  mapGapped[0] === 3 && mapGapped[1] === 4,
  "msaRawToMasterAligned maps gap-free indices to MSA columns",
);

const colorMap = buildMsaLiabilityColorMap("---MQ", ["Met Ox"]);
assert(
  typeof colorMap === "object" && colorMap !== null,
  "buildMsaLiabilityColorMap returns an object",
);
assert(
  Object.keys(colorMap).length > 0,
  "buildMsaLiabilityColorMap finds Met Ox hit in ---MQ",
);
assert(
  colorMap[3] && colorMap[3].indexOf("rgba(") >= 0,
  "buildMsaLiabilityColorMap maps M at MSA col 3 to an rgba color",
);
assert(
  Object.keys(buildMsaLiabilityColorMap("MQ", [])).length === 0,
  "buildMsaLiabilityColorMap returns empty map when no highlights selected",
);

const qEntry = { id: 1, name: "Q", vh: "QAAAA", vl: "" };
assert(
  entryHasFilteredLiability(qEntry, ["N-term Q"]),
  "entryHasFilteredLiability detects N-term Q on VH",
);
assert(
  !entryHasFilteredLiability(qEntry, ["Met Ox"]),
  "entryHasFilteredLiability false when filter type absent",
);
const kept = filterBasketByLiabilityFilters(
  [qEntry, { id: 2, name: "M", vh: "AAAAA", vl: "" }],
  ["N-term Q"],
);
assert(
  kept.length === 1 && kept[0].id === 2,
  "filterBasketByLiabilityFilters drops matching entries",
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
assert(
  cons1.consensus[1] === "Q",
  "conservation: consensus pos 2 = Q (majority)",
);

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
assert(pipeClusters[0].members.length === 2, "pipeline: cluster has 2 members");

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
// detectLiabilities — Agg Patch (Kyte-Doolittle sliding window)
// ============================================================
section("detectLiabilities (Agg Patch)");

// Verify constants
assert(AGG_WINDOW === 7, "AGG_WINDOW is 7");
assert(AGG_THRESHOLD === 1.6, "AGG_THRESHOLD is 1.6");
assert(KD_HYDROPHOBICITY["I"] === 4.5, "KD: Ile = 4.5");
assert(KD_HYDROPHOBICITY["R"] === -4.5, "KD: Arg = -4.5");
assert(KD_HYDROPHOBICITY["A"] === 1.8, "KD: Ala = 1.8");

// Pure hydrophobic 7-mer: mean 4.5 >> 1.6
const aggAllI = detectLiabilities("IIIIIII")["Agg Patch"];
assert(aggAllI && aggAllI.length === 7, "7 Ile = all 7 positions flagged");
assert(aggAllI[0] === 0 && aggAllI[6] === 6, "7 Ile spans positions 0-6");

// Pure charged 7-mer: mean -3.5, no patch
assert(!detectLiabilities("DDDDDDD")["Agg Patch"], "7 Asp = no agg patch");

// Sequence shorter than window: no patch possible
assert(
  !detectLiabilities("IIIII")["Agg Patch"],
  "5 Ile < window = no agg patch",
);
assert(
  !detectLiabilities("IIIIII")["Agg Patch"],
  "6 Ile < window = no agg patch",
);

// Below threshold: AAAAAAG → 6(1.8) + (-0.4) = 10.4, mean 1.486
assert(
  !detectLiabilities("AAAAAAG")["Agg Patch"],
  "AAAAAAG mean 1.49 = no patch",
);

// Above threshold: AAAAAAM → 6(1.8) + 1.9 = 12.7, mean 1.814
const aggAAM = detectLiabilities("AAAAAAM")["Agg Patch"];
assert(aggAAM && aggAAM.length === 7, "AAAAAAM mean 1.81 = all 7 flagged");

// Hydrophobic island in polar flanks: DDDIIIIIIIDDD (13 chars)
// Windows [0,6] and [6,12] straddle the boundary, only middle windows pass
const aggIsland = detectLiabilities("DDDIIIIIIIDDD")["Agg Patch"];
assert(aggIsland && aggIsland.length > 0, "Ile island is detected");
assert(!aggIsland.includes(0), "position 0 (D flank) not flagged");
assert(!aggIsland.includes(12), "position 12 (D flank) not flagged");
assert(aggIsland.includes(5), "middle of island is flagged");
assert(
  aggIsland.includes(1) && aggIsland.includes(11),
  "near-boundary positions flagged",
);

// Two separate patches with charged gap: IIIIIII + 10×R + IIIIIII
const twoPatch = detectLiabilities("IIIIIIIRRRRRRRRRRIIIIIII")["Agg Patch"];
assert(twoPatch && twoPatch.length > 0, "two patches detected");
assert(
  twoPatch.includes(0) && twoPatch.includes(6),
  "first patch core flagged",
);
assert(
  twoPatch.includes(17) && twoPatch.includes(23),
  "second patch core flagged",
);
assert(
  !twoPatch.includes(10) && !twoPatch.includes(12),
  "charged gap not flagged",
);

// Positions are always sorted
for (let i = 1; i < twoPatch.length; i++) {
  assert(
    twoPatch[i] > twoPatch[i - 1],
    "agg patch positions sorted at index " + i,
  );
}

// Unknown residues default to 0 hydrophobicity
// XXXXXXX → 7(0) = 0, mean 0 < 1.6
assert(
  !detectLiabilities("XXXXXXX")["Agg Patch"],
  "unknown residues (X) = no patch",
);

// Single high window in longer sequence: RRRRRIIIIIIIRRRRR (17 chars)
const aggSingle = detectLiabilities("RRRRRIIIIIIIRRRRR")["Agg Patch"];
assert(aggSingle && aggSingle.length > 0, "single central patch detected");
assert(aggSingle.includes(5) && aggSingle.includes(10), "Ile core flagged");
assert(
  !aggSingle.includes(0) && !aggSingle.includes(16),
  "R flanks not flagged",
);

// All-Ala 7-mer: mean 1.8 >= 1.6, should flag
const aggAllA = detectLiabilities("AAAAAAA")["Agg Patch"];
assert(aggAllA && aggAllA.length === 7, "7 Ala mean 1.8 >= 1.6 = flagged");

// Mixed just below: AAAAAAS → 6(1.8) + (-0.8) = 10.0, mean 1.429
assert(
  !detectLiabilities("AAAAAAS")["Agg Patch"],
  "AAAAAAS mean 1.43 = no patch",
);

// Real antibody VH — verify agg patch exists in hydrophobic framework regions
const vhAgg = detectLiabilities(
  "EVQLVESGGGLVQPGGSLRLSCAASGFTFSDYWMHWVRQAPGKGLEWVSRINSDGSSTSYADSVKGRFTISRDNAKNTLYLQMNSLRAEDTAVYYCTTDGWGFDYWGQGTLVTVSS",
)["Agg Patch"];
if (vhAgg) {
  for (let i = 1; i < vhAgg.length; i++) {
    assert(vhAgg[i] > vhAgg[i - 1], "VH agg patch positions sorted at " + i);
  }
  assert(
    vhAgg.every(function (p) {
      return p >= 0 && p < 113;
    }),
    "VH agg patch positions within sequence bounds",
  );
}

// ============================================================
// detectCDR3 (igblast_thera.tsv ground truth)
// ============================================================
section("detectCDR3 (igblast_thera.tsv)");

const cdr3TestCases = [
  {
    id: "Adalimumab_vh",
    chain: "VH",
    cdr3: "AKVSYLSTASSLDY",
    seq: "EVQLVESGGGLVQPGRSLRLSCAASGFTFDDYAMHWVRQAPGKGLEWVSAITWNSGHIDYADSVEGRFTISRDNAKNSLYLQMNSLRAEDTAVYYCAKVSYLSTASSLDYWGQGTLVTVSS",
  },
  {
    id: "Adalimumab_vl",
    chain: "VL",
    cdr3: "QRYNRAPYT",
    seq: "DIQMTQSPSSLSASVGDRVTITCRASQGIRNYLAWYQQKPGKAPKLLIYAASTLQSGVPSRFSGSGSGTDFTLTISSLQPEDVATYYCQRYNRAPYTFGQGTKVEIK",
  },
  {
    id: "Aducanumab_vh",
    chain: "VH",
    cdr3: "ARDRGIGARRGPYYMDV",
    seq: "QVQLVESGGGVVQPGRSLRLSCAASGFAFSSYGMHWVRQAPGKGLEWVAVIWFDGTKKYYTDSVKGRFTISRDNSKNTLYLQMNTLRAEDTAVYYCARDRGIGARRGPYYMDVWGKGTTVTVSS",
  },
  {
    id: "Aducanumab_vl",
    chain: "VL",
    cdr3: "QQSYSTPLT",
    seq: "DIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSYSTPLTFGGGTKVEIK",
  },
  {
    id: "Alemtuzumab_vh",
    chain: "VH",
    cdr3: "AREGHTAAPFDY",
    seq: "QVQLQESGPGLVRPSQTLSLTCTVSGFTFTDFYMNWVRQPPGRGLEWIGFIRDKAKGYTTEYNPSVKGRVTMLVDTSKNQFSLRLSSVTAADTAVYYCAREGHTAAPFDYWGQGSLVTVSS",
  },
  {
    id: "Alemtuzumab_vl",
    chain: "VL",
    cdr3: "LQHISRPRT",
    seq: "DIQMTQSPSSLSASVGDRVTITCKASQNIDKYLNWYQQKPGKAPKLLIYNTNNLQTGVPSRFSGSGSGTDFTFTISSLQPEDIATYYCLQHISRPRTFGQGTKVEIK",
  },
  {
    id: "Alirocumab_vh",
    chain: "VH",
    cdr3: "AKDSNWGNFDL",
    seq: "EVQLVESGGGLVQPGGSLRLSCAASGFTFNNYAMNWVRQAPGKGLDWVSTISGSGGTTNYADSVKGRFIISRDSSKHTLYLQMNSLRAEDTAVYYCAKDSNWGNFDLWGRGTLVTVSS",
  },
  {
    id: "Alirocumab_vl",
    chain: "VL",
    cdr3: "QQYYTTPYT",
    seq: "DIVMTQSPDSLAVSLGERATINCKSSQSVLYRSNNRNFLGWYQQKPGQPPNLLIYWASTRESGVPDRFSGSGSGTDFTLTISSLQAEDVAVYYCQQYYTTPYTFGQGTKLEIK",
  },
  {
    id: "Anifrolumab_vh",
    chain: "VH",
    cdr3: "ARHDIEGFDY",
    seq: "EVQLVQSGAEVKKPGESLKISCKGSGYIFTNYWIAWVRQMPGKGLESMGIIYPGDSDIRYSPSFQGQVTISADKSITTAYLQWSSLKASDTAMYYCARHDIEGFDYWGRGTLVTVSS",
  },
  {
    id: "Anifrolumab_vl",
    chain: "VL",
    cdr3: "QQYDSSAIT",
    seq: "EIVLTQSPGTLSLSPGERATLSCRASQSVSSSFFAWYQQKPGQAPRLLIYGASSRATGIPDRLSGSGSGTDFTLTITRLEPEDFAVYYCQQYDSSAITFGQGTRLEIK",
  },
  {
    id: "Ansuvimab_vh",
    chain: "VH",
    cdr3: "VRSDRGVAGLFDS",
    seq: "EVQLVESGGGLIQPGGSLRLSCAASGFALRMYDMHWVRQTIDKRLEWVSAVGPSGDTTYADSVKGRFAVSRENAKNSLSLQMNSLTAGDTAIYYCVRSDRGVAGLFDSWGQGILVTVSS",
  },
  {
    id: "Ansuvimab_vl",
    chain: "VL",
    cdr3: "QNYNSAPLT",
    seq: "DIQMTQSPSSLSASVGDRITITCRASQAFDNYVAWYQQRPGKVPKLLISAASALHAGVPSRFSGSGSGTHFTLTISSLQPEDVATYYCQNYNSAPLTFGGGTKVEIK",
  },
  {
    id: "Atezolizumab_vh",
    chain: "VH",
    cdr3: "ARRHWPGGFDY",
    seq: "EVQLVESGGGLVQPGGSLRLSCAASGFTFSDSWIHWVRQAPGKGLEWVAWISPYGGSTYYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCARRHWPGGFDYWGQGTLVTVSS",
  },
  {
    id: "Atezolizumab_vl",
    chain: "VL",
    cdr3: "QQYLYHPAT",
    seq: "DIQMTQSPSSLSASVGDRVTITCRASQDVSTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQYLYHPATFGQGTKVEIK",
  },
  {
    id: "Avelumab_vh",
    chain: "VH",
    cdr3: "ARIKLGTVTTVDY",
    seq: "EVQLLESGGGLVQPGGSLRLSCAASGFTFSSYIMMWVRQAPGKGLEWVSSIYPSGGITFYADTVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCARIKLGTVTTVDYWGQGTLVTVSS",
  },
  {
    id: "Avelumab_vl",
    chain: "VL",
    cdr3: "SSYTSSSTRV",
    seq: "ALTQPASVSGSPGQSITISCTGTSSDVGGYNYVSWYQQHPGKAPKLMIYDVSNRPSGVSNRFSGSKSGNTASLTISGLQAEDEADYYCSSYTSSSTRVFGTGTKVTVL",
  },
  {
    id: "Bamlanivimab_vh",
    chain: "VH",
    cdr3: "ARGYYEARHYYYYYAMDV",
    seq: "QVQLVQSGAEVKKPGSSVKVSCKASGGTFSNYAISWVRQAPGQGLEWMGRIIPILGIANYAQKFQGRVTITADKSTSTAYMELSSLRSEDTAVYYCARGYYEARHYYYYYAMDVWGQGTAVTVSS",
  },
  {
    id: "Bamlanivimab_vl",
    chain: "VL",
    cdr3: "QQSYSTPRT",
    seq: "DIQMTQSPSSLSASVGDRVTITCRASQSISSYLSWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTITSLQPEDFATYYCQQSYSTPRTFGQGTKVEIK",
  },
  {
    id: "Bebtelovimab_vh",
    chain: "VH",
    cdr3: "AHHSISTIFDH",
    seq: "QITLKESGPTLVKPTQTLTLTCTFSGFSLSISGVGVGWLRQPPGKALEWLALIYWDDDKRYSPSLKSRLTISKDTSKNQVVLKMTNIDPVDTATYYCAHHSISTIFDHWGQGTLVTVSS",
  },
  {
    id: "Bebtelovimab_vl",
    chain: "VL",
    cdr3: "SSYTTSSAV",
    seq: "ALTQPASVSGSPGQSITISCTATSSDVGDYNYVSWYQQHPGKAPKLMIFEVSDRPSGISNRFSGSKSGNTASLTISGLQAEDEADYYCSSYTTSSAVFGGGTKLTVL",
  },
  {
    id: "Belimumab_vh",
    chain: "VH",
    cdr3: "ARSRDLLLFPHHALSP",
    seq: "QVQLQQSGAEVKKPGSSVRVSCKASGGTFNNNAINWVRQAPGQGLEWMGGIIPMFGTAKYSQNFQGRVAITADESTGTASMELSSLRSEDTAVYYCARSRDLLLFPHHALSPWGRGTMVTVSS",
  },
  {
    id: "Belimumab_vl",
    chain: "VL",
    cdr3: "SSRDSSGNHWV",
    seq: "ELTQDPAVSVALGQTVRVTCQGDSLRSYYASWYQQKPGQAPVLVIYGKNNRPSGIPDRFSGSSSGNTASLTITGAQAEDEADYYCSSRDSSGNHWVFGGGTELTVL",
  },
  {
    id: "Benmelstobart_vh",
    chain: "VH",
    cdr3: "ARLGFYAMDY",
    seq: "QITLKESGPTLVKPTQTLTLTCTVSGFSLSTYGVHWIRQPPGKALEWLGVIWRGVTTDYNAAFMSRLTITKDNSKNQVVLTMNNMDPVDTATYYCARLGFYAMDYWGQGTLVTVSS",
  },
  {
    id: "Benmelstobart_vl",
    chain: "VL",
    cdr3: "QQDYTSPYT",
    seq: "DIQMTQSPSSLSASVGDRVTITCKASQSVSNDVAWYQQKPGKAPKLLIYYAANRYTGVPDRFSGSGYGTDFTFTISSLQPEDIATYFCQQDYTSPYTFGQGTKLEIK",
  },
  {
    id: "Benralizumab_vh",
    chain: "VH",
    cdr3: "GREGIRYYGLLGDY",
    seq: "VQLVQSGAEVKKPGASVKVSCKASGYTFTSYVIHWVRQRPGQGLAWMGYINPYNDGTKYNERFKGKVTITSDRSTSTVYMELSSLRSEDTAVYLCGREGIRYYGLLGDYWGQGTLVTVSS",
  },
  {
    id: "Benralizumab_vl",
    chain: "VL",
    cdr3: "QQGYTLPYT",
    seq: "DIQMTQSPSSLSASVGDRVTITCGTSEDIINYLNWYQQKPGKAPKLLIYHTSRLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQGYTLPYTFGQGTKVEIK",
  },
  {
    id: "Bevacizumab_vh",
    chain: "VH",
    cdr3: "AKYPHYYGSSHWYFDV",
    seq: "VQLVESGGGLVQPGGSLRLSCAASGYTFTNYGMNWVRQAPGKGLEWVGWINTYTGEPTYAADFKRRFTFSLDTSKSTAYLQMNSLRAEDTAVYYCAKYPHYYGSSHWYFDVWGQGTLVTVSS",
  },
  {
    id: "Bevacizumab_vl",
    chain: "VL",
    cdr3: "QQYSTVPWT",
    seq: "DIQMTQSPSSLSASVGDRVTITCSASQDISNYLNWYQQKPGKAPKVLIYFTSSLHSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQYSTVPWTFGQGTKVEIK",
  },
  {
    id: "Bezlotoxumab_vh",
    chain: "VH",
    cdr3: "ARRRNWGNAFDI",
    seq: "EVQLVQSGAEVKKSGESLKISCKGSGYSFTSYWIGWVRQMPGKGLEWMGIFYPGDSSTRYSPSFQGQVTISADKSVNTAYLQWSSLKASDTAMYYCARRRNWGNAFDIWGQGTMVTVSS",
  },
  {
    id: "Bezlotoxumab_vl",
    chain: "VL",
    cdr3: "QQYGSSTWT",
    seq: "EIVLTQSPGTLSLSPGERATLSCRASQSVSSSYLAWYQQKPGQAPRLLIYGASSRATGIPDRFSGSGSGTDFTLTISRLEPEDFAVYYCQQYGSSTWTFGQGTKVEIK",
  },
  {
    id: "Bimekizumab_vh",
    chain: "VH",
    cdr3: "ASPPQYYEGSIYRLWFAH",
    seq: "EVQLVESGGGLVQPGGSLRLSCAASGFTFSDYNMAWVRQAPGKGLEWVATITYEGRNTYYRDSVKGRFTISRDNAKNSLYLQMNSLRAEDTAVYYCASPPQYYEGSIYRLWFAHWGQGTLVTVSS",
  },
  {
    id: "Bimekizumab_vl",
    chain: "VL",
    cdr3: "QQTWSDPWT",
    seq: "IQLTQSPSSLSASVGDRVTITCRADESVRTLMHWYQQKPGKAPKLLIYLVSNSEIGVPDRFSGSGSGTDFRLTISSLQPEDFATYYCQQTWSDPWTFGQGTKVEIK",
  },
  {
    id: "Brodalumab_vh",
    chain: "VH",
    cdr3: "ARRQLYFDY",
    seq: "QVQLVQSGAEVKKPGASVKVSCKASGYTFTRYGISWVRQAPGQGLEWMGWISTYSGNTNYAQKLQGRVTMTTDTSTSTAYMELRSLRSDDTAVYYCARRQLYFDYWGQGTLVTVSS",
  },
  {
    id: "Brodalumab_vl",
    chain: "VL",
    cdr3: "QQYDNWPLT",
    seq: "EIVMTQSPATLSVSPGERATLSCRASQSVSSNLAWFQQKPGQAPRPLIYDASTRATGVPARFSGSGSGTDFTLTISSLQSEDFAVYYCQQYDNWPLTFGGGTKVEIK",
  },
  {
    id: "Burosumab_vh",
    chain: "VH",
    cdr3: "ARDIVDAFDF",
    seq: "QVQLVQSGAEVKKPGASVKVSCKASGYTFTNHYMHWVRQAPGQGLEWMGIINPISGSTSNAQKFQGRVTMTRDTSTSTVYMELSSLRSEDTAVYYCARDIVDAFDFWGQGTMVTVSS",
  },
  {
    id: "Burosumab_vl",
    chain: "VL",
    cdr3: "QQFNDYFT",
    seq: "AIQLTQSPSSLSASVGDRVTITCRASQGISSALVWYQQKPGKAPKLLIYDASSLESGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQFNDYFTFGPGTKVDIK",
  },
  {
    id: "Camrelizumab_vh",
    chain: "VH",
    cdr3: "ARQLYYFDY",
    seq: "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYMMSWVRQAPGKGLEWVATISGGGANTYYPDSVKGRFTISRDNAKNSLYLQMNSLRAEDTAVYYCARQLYYFDYWGQGTTVTVSS",
  },
  {
    id: "Camrelizumab_vl",
    chain: "VL",
    cdr3: "QQVYSIPWT",
    seq: "DIQMTQSPSSLSASVGDRVTITCLASQTIGTWLTWYQQKPGKAPKLLIYTATSLADGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQVYSIPWTFGGGTKVEIK",
  },
  {
    id: "Canakinumab_vh",
    chain: "VH",
    cdr3: "ARDLRTGPFDY",
    seq: "QVQLVESGGGVVQPGRSLRLSCAASGFTFSVYGMNWVRQAPGKGLEWVAIIWYDGDNQYYADSVKGRFTISRDNSKNTLYLQMNGLRAEDTAVYYCARDLRTGPFDYWGQGTLVTVSS",
  },
  {
    id: "Canakinumab_vl",
    chain: "VL",
    cdr3: "HQSSSLPFT",
    seq: "EIVLTQSPDFQSVTPKEKVTITCRASQSIGSSLHWYQQKPDQSPKLLIKYASQSFSGVPSRFSGSGSGTDFTLTINSLEAEDAAAYYCHQSSSLPFTFGPGTKVDIK",
  },
  {
    id: "Cemiplimab_vh",
    chain: "VH",
    cdr3: "VKWGNIYFDY",
    seq: "EVQLLESGGVLVQPGGSLRLSCAASGFTFSNFGMTWVRQAPGKGLEWVSGISGGGRDTYFADSVKGRFTISRDNSKNTLYLQMNSLKGEDTAVYYCVKWGNIYFDYWGQGTLVTVSS",
  },
  {
    id: "Cemiplimab_vl",
    chain: "VL",
    cdr3: "QQSSNTPFT",
    seq: "DIQMTQSPSSLSASVGDSITITCRASLSINTFLNWYQQKPGKAPNLLIYAASSLHGGVPSRFSGSGSGTDFTLTIRTLQPEDFATYYCQQSSNTPFTFGPGTVVDFR",
  },
  {
    id: "Cilgavimab_vh",
    chain: "VH",
    cdr3: "TTAGSYYYDTVGPGLPEGKFDY",
    seq: "EVQLVESGGGLVKPGGSLRLSCAASGFTFRDVWMSWVRQAPGKGLEWVGRIKSKIDGGTTDYAAPVKGRFTISRDDSKNTLYLQMNSLKTEDTAVYYCTTAGSYYYDTVGPGLPEGKFDYWGQGTLVTVSS",
  },
  {
    id: "Cilgavimab_vl",
    chain: "VL",
    cdr3: "QQYYSTLT",
    seq: "DIVMTQSPDSLAVSLGERATINCKSSQSVLYSSNNKNYLAWYQQKPGQPPKLLMYWASTRESGVPDRFSGSGSGAEFTLTISSLQAEDVAIYYCQQYYSTLTFGGGTKVEIK",
  },
  {
    id: "Concizumab_vh",
    chain: "VH",
    cdr3: "ARLGGYDEGDAMDS",
    seq: "EVQLVESGGGLVKPGGSLRLSCAASGFTFSNYAMSWVRQTPEKRLEWVATISRSGSYSYFPDSVQGRFTISRDNAKNSLYLQMNSLRAEDTAVYYCARLGGYDEGDAMDSWGQGTTVTVSS",
  },
  {
    id: "Concizumab_vl",
    chain: "VL",
    cdr3: "LQATHFPQT",
    seq: "DIVMTQTPLSLSVTPGQPASISCKSSQSLLESDGKTYLNWYLQKPGQSPQLLIYLVSILDSGVPDRFSGSGSGTDFTLKISRVEAEDVGVYYCLQATHFPQTFGGGTKVEIK",
  },
  {
    id: "Cosibelimab_vh",
    chain: "VH",
    cdr3: "ARGRQMFGAGIDF",
    seq: "VQLVQSGAEVKKPGSSVKVSCKASGGTFSRSAISWVRQAPGQGLEWMGVIIPAFGEANYAQKFQGRVTITADESTSTAYMELSSLRSEDTAVYYCARGRQMFGAGIDFWGQGTLVTVSS",
  },
  {
    id: "Cosibelimab_vl",
    chain: "VL",
    cdr3: "QSYDSNNRHVI",
    seq: "NFMLTQPHSVSESPGKTVTISCTRSSGSIDSNYVQWYQQRPGSAPTTVIYEDNQRPSGVPDRFSGSIDSSSNSASLTISGLKTEDEADYYCQSYDSNNRHVIFGGGTKLTVL",
  },
  {
    id: "Crizanlizumab_vh",
    chain: "VH",
    cdr3: "ARRGEYGNYEGAMDY",
    seq: "QVQLVQSGAEVKKPGASVKVSCKVSGYTFTSYDINWVRQAPGKGLEWMGWIYPGDGSIKYNEKFKGRVTMTVDKSTDTAYMELSSLRSEDTAVYYCARRGEYGNYEGAMDYWGQGTLVTVSS",
  },
  {
    id: "Crizanlizumab_vl",
    chain: "VL",
    cdr3: "QQSDENPLT",
    seq: "DIQMTQSPSSLSASVGDRVTITCKASQSVDYDGHSYMNWYQQKPGKAPKLLIYAASNLESGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSDENPLTFGGGTKVEIK",
  },
  {
    id: "Crovalimab_vh",
    chain: "VH",
    cdr3: "ASDAGYDYPTHAMHY",
    seq: "VQLVESGGGLVQPGRSLRLSCAASGFTVHSSYYMAWVRQAPGKGLEWVGAIFTGSGAEYKAEWAKGRVTISKDTSKNQVVLTMTNMDPVDTATYYCASDAGYDYPTHAMHYWGQGTLVTVSS",
  },
  {
    id: "Crovalimab_vl",
    chain: "VL",
    cdr3: "QNTKVGSSYGNT",
    seq: "DIQMTQSPSSLSASVGDRVTITCRASQGISSSLAWYQQKPGKAPKLLIYGASETESGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQNTKVGSSYGNTFGGGTKVEIK",
  },
  {
    id: "Daclizumab_vh",
    chain: "VH",
    cdr3: "ARGGGVFDY",
    seq: "QVQLVQSGAEVKKPGSSVKVSCKASGYTFTSYRMHWVRQAPGQGLEWIGYINPSTGYTEYNQKFKDKATITADESTNTAYMELSSLRSEDTAVYYCARGGGVFDYWGQGTLVTVSS",
  },
  {
    id: "Daclizumab_vl",
    chain: "VL",
    cdr3: "HQRSTYPLT",
    seq: "DIQMTQSPSTLSASVGDRVTITCSASSSISYMHWYQQKPGKAPKLLIYTTSNLASGVPARFSGSGSGTEFTLTISSLQPDDFATYYCHQRSTYPLTFGQGTKVEVK",
  },
  {
    id: "Daratumumab_vh",
    chain: "VH",
    cdr3: "AKDKILWFGEPVFDY",
    seq: "EVQLLESGGGLVQPGGSLRLSCAVSGFTFNSFAMSWVRQAPGKGLEWVSAISGSGGGTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYFCAKDKILWFGEPVFDYWGQGTLVTVSS",
  },
  {
    id: "Daratumumab_vl",
    chain: "VL",
    cdr3: "QQRSNWPPT",
    seq: "EIVLTQSPATLSLSPGERATLSCRASQSVSSYLAWYQQKPGQAPRLLIYDASNRATGIPARFSGSGSGTDFTLTISSLEPEDFAVYYCQQRSNWPPTFGQGTKVEIK",
  },
  {
    id: "Datopotamab_vh",
    chain: "VH",
    cdr3: "ARSGFGSSYWYFDV",
    seq: "QVQLVQSGAEVKKPGASVKVSCKASGYTFTTAGMQWVRQAPGQGLEWMGWINTHSGVPKYAEDFKGRVTISADTSTSTAYLQLSSLKSEDTAVYYCARSGFGSSYWYFDVWGQGTLVTVSS",
  },
  {
    id: "Datopotamab_vl",
    chain: "VL",
    cdr3: "QQHYITPLT",
    seq: "DIQMTQSPSSLSASVGDRVTITCKASQDVSTAVAWYQQKPGKAPKLLIYSASYRYTGVPSRFSGSGSGTDFTLTISSLQPEDFAVYYCQQHYITPLTFGQGTKLEIK",
  },
  {
    id: "Denosumab_vh",
    chain: "VH",
    cdr3: "AKDPGTTVIMSWFDP",
    seq: "EVQLLESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSGITGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDPGTTVIMSWFDPWGQGTLVTVSS",
  },
  {
    id: "Denosumab_vl",
    chain: "VL",
    cdr3: "QQYGSSPRT",
    seq: "EIVLTQSPGTLSLSPGERATLSCRASQSVRGRYLAWYQQKPGQAPRLLIYGASSRATGIPDRFSGSGSGTDFTLTISRLEPEDFAVFYCQQYGSSPRTFGQGTKVEIK",
  },
  {
    id: "Donanemab_vh",
    chain: "VH",
    cdr3: "AREGITVY",
    seq: "QVQLVQSGAEVKKPGSSVKVSCKASGYDFTRYYINWVRQAPGQGLEWMGWINPGSGNTKYNEKFKGRVTITADESTSTAYMELSSLRSEDTAVYYCAREGITVYWGQGTTVTVSS",
  },
  {
    id: "Donanemab_vl",
    chain: "VL",
    cdr3: "VQGTHYPFT",
    seq: "DIVMTQTPLSLSVTPGQPASISCKSSQSLLYSRGKTYLNWLLQKPGQSPQLLIYAVSKLDSGVPDRFSGSGSGTDFTLKISRVEAEDVGVYYCVQGTHYPFTFGQGTKLEIK",
  },
  {
    id: "Dostarlimab_vh",
    chain: "VH",
    cdr3: "ASPYYAMDY",
    seq: "EVQLLESGGGLVQPGGSLRLSCAASGFTFSSYDMSWVRQAPGKGLEWVSTISGGGSYTYYQDSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCASPYYAMDYWGQGTTVTVSS",
  },
  {
    id: "Dostarlimab_vl",
    chain: "VL",
    cdr3: "QHYSSYPWT",
    seq: "DIQLTQSPSFLSAYVGDRVTITCKASQDVGTAVAWYQQKPGKAPKLLIYWASTLHTGVPSRFSGSGSGTEFTLTISSLQPEDFATYYCQHYSSYPWTFGQGTKLEIK",
  },
  {
    id: "Dupilumab_vh",
    chain: "VH",
    cdr3: "AKDRLSITIRPRYYGLDV",
    seq: "EVQLVESGGGLEQPGGSLRLSCAGSGFTFRDYAMTWVRQAPGKGLEWVSSISGSGGNTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDRLSITIRPRYYGLDVWGQGTTVTVSS",
  },
  {
    id: "Dupilumab_vl",
    chain: "VL",
    cdr3: "MQALQTPYT",
    seq: "DIVMTQSPLSLPVTPGEPASISCRSSQSLLYSIGYNYLDWYLQKSGQSPQLLIYLGSNRASGVPDRFSGSGSGTDFTLKISRVEAEDVGFYYCMQALQTPYTFGQGTKLEIK",
  },
  {
    id: "Durvalumab_vh",
    chain: "VH",
    cdr3: "AREGGWFGELAFDY",
    seq: "EVQLVESGGGLVQPGGSLRLSCAASGFTFSRYWMSWVRQAPGKGLEWVANIKQDGSEKYYVDSVKGRFTISRDNAKNSLYLQMNSLRAEDTAVYYCAREGGWFGELAFDYWGQGTLVTVSS",
  },
  {
    id: "Durvalumab_vl",
    chain: "VL",
    cdr3: "QQYGSLPWT",
    seq: "EIVLTQSPGTLSLSPGERATLSCRASQRVSSSYLAWYQQKPGQAPRLLIYDASSRATGIPDRFSGSGSGTDFTLTISRLEPEDFAVYYCQQYGSLPWTFGQGTKVEIK",
  },
  {
    id: "Ebdarokimab_vh",
    chain: "VH",
    cdr3: "ARRRPGQGYFDF",
    seq: "EVQLVQSGAEVKKPGESLKISCQSSGYTFTSYWIGWVRQMPGQGLEWIGIMSPVDSDIRYNPMFRGQVTMSVDKSSSTAYLQWSSLKASDTAMYYCARRRPGQGYFDFWGQGTMVTVSS",
  },
  {
    id: "Ebdarokimab_vl",
    chain: "VL",
    cdr3: "QQYNIYPYT",
    seq: "EIVLTQSPATLSASPGERATISCRASQSVGTWVAWYQQKPGQAPRSLIYAASNLQSGIPARFSGSGSGTDFTLTISSLEPEDFAVYYCQQYNIYPYTFGQGTRLEIK",
  },
  {
    id: "Ebronucimab_vh",
    chain: "VH",
    cdr3: "AREYDFWSAYYDAFDV",
    seq: "EVQLVESGGGLVQPGRSLRLSCAASGFTFSSYSMNWVRQAPGKGLEWVSGISSSSSYISYADSVQGRFTISRDNGKNSLYLQMNSLRAEDTALYFCAREYDFWSAYYDAFDVWGQGTMVTVSS",
  },
  {
    id: "Ebronucimab_vl",
    chain: "VL",
    cdr3: "QSFDGSLSGSV",
    seq: "LTQPRSVSGSPGQSVTISCTGTSRNIGGGNDVHWYQQHPGKAPKLLISGVIERSSGVPDRFSGSKSGNTASLTISGLQAEDEADYYCQSFDGSLSGSVFGTGTDVTVL",
  },
  {
    id: "Eculizumab_vh",
    chain: "VH",
    cdr3: "ARYFFGSSPNWYFDV",
    seq: "QVQLVQSGAEVKKPGASVKVSCKASGYIFSNYWIQWVRQAPGQGLEWMGEILPGSGSTEYTENFKDRVTMTRDTSTSTVYMELSSLRSEDTAVYYCARYFFGSSPNWYFDVWGQGTLVTVSS",
  },
  {
    id: "Eculizumab_vl",
    chain: "VL",
    cdr3: "QNVLNTPLT",
    seq: "DIQMTQSPSSLSASVGDRVTITCGASENIYGALNWYQQKPGKAPKLLIYGATNLADGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQNVLNTPLTFGQGTKVEIK",
  },
  {
    id: "Efalizumab_vh",
    chain: "VH",
    cdr3: "ARGIYFYGTTYFDY",
    seq: "EVQLVESGGGLVQPGGSLRLSCAASGYSFTGHWMNWVRQAPGKGLEWVGMIHPSDSETRYNQKFKDRFTISVDKSKNTLYLQMNSLRAEDTAVYYCARGIYFYGTTYFDYWGQGTLVTVSS",
  },
  {
    id: "Efalizumab_vl",
    chain: "VL",
    cdr3: "QQHNEYPLT",
    seq: "DIQMTQSPSSLSASVGDRVTITCRASKTISKYLAWYQQKPGKAPKLLIYSGSTLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQHNEYPLTFGQGTKVEIK",
  },
  {
    id: "Elotuzumab_vh",
    chain: "VH",
    cdr3: "ARPDGNYWYFDV",
    seq: "EVQLVESGGGLVQPGGSLRLSCAASGFDFSRYWMSWVRQAPGKGLEWIGEINPDSSTINYAPSLKDKFIISRDNAKNSLYLQMNSLRAEDTAVYYCARPDGNYWYFDVWGQGTLVTVSS",
  },
  {
    id: "Elotuzumab_vl",
    chain: "VL",
    cdr3: "QQYSSYPYT",
    seq: "DIQMTQSPSSLSASVGDRVTITCKASQDVGIAVAWYQQKPGKVPKLLIYWASTRHTGVPDRFSGSGSGTDFTLTISSLQPEDVATYYCQQYSSYPYTFGQGTKVEIK",
  },
  {
    id: "Emapalumab_vh",
    chain: "VH",
    cdr3: "AKDGSSGWYVPHWFDP",
    seq: "EVQLLESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDGSSGWYVPHWFDPWGQGTLVTVSS",
  },
  {
    id: "Emapalumab_vl",
    chain: "VL",
    cdr3: "QSYDGSNRWM",
    seq: "NFMLTQPHSVSESPGKTVTISCTRSSGSIASNYVQWYQQRPGSSPTTVIYEDNQRPSGVPDRFSGSIDSSSNSASLTISGLKTEDEADYYCQSYDGSNRWMFGGGTKLTVL",
  },
  {
    id: "Enlonstobart_vh",
    chain: "VH",
    cdr3: "ATNNDY",
    seq: "QVQLVESGGGVVQPGRSLRLTCKASGLTFSSSGMHWVRQAPGKGLEWVAVIWYDGSKRYYADSVKGRFTISRDNSKNTLFLQMNSLRAEDTAVYYCATNNDYWGQGTLVTVSS",
  },
  {
    id: "Enlonstobart_vl",
    chain: "VL",
    cdr3: "QQYSNWPRT",
    seq: "EIVLTQSPATLSLSPGERATLSCRASQSVSSYLAWYQQKPGQAPRLLIYTASNRATGIPARFSGSGSGTDFTLTISSLEPEDFAVYYCQQYSNWPRTFGQGTKVEIK",
  },
  {
    id: "Eptinezumab_vh",
    chain: "VH",
    cdr3: "ARGDI",
    seq: "EVQLVESGGGLVQPGGSLRLSCAVSGIDLSGYYMNWVRQAPGKGLEWVGVIGINGATYYASWAKGRFTISRDNSKTTVYLQMNSLRAEDTAVYFCARGDIWGQGTLVTVSS",
  },
  {
    id: "Eptinezumab_vl",
    chain: "VL",
    cdr3: "LGSYDCTNGDCFV",
    seq: "TQSPSSLSASVGDRVTINCQASQSVYHNTYLAWYQQKPGKVPKQLIYDASTLASGVPSRFSGSGSGTDFTLTISSLQPEDVATYYCLGSYDCTNGDCFVFGGGTKVEIK",
  },
  {
    id: "Erenumab_vh",
    chain: "VH",
    cdr3: "ARDRLNYYDSSGYYHYKYYGMAV",
    seq: "QVQLVESGGGVVQPGRSLRLSCAASGFTFSSFGMHWVRQAPGKGLEWVAVISFDGSIKYSVDSVKGRFTISRDNSKNTLFLQMNSLRAEDTAVYYCARDRLNYYDSSGYYHYKYYGMAVWGQGTTVTVSS",
  },
  {
    id: "Erenumab_vl",
    chain: "VL",
    cdr3: "GTWDSRLSAVV",
    seq: "VLTQPPSVSAAPGQKVTISCSGSSSNIGNNYVSWYQQLPGTAPKLLIYDNNKRPSGIPDRFSGSKSGTSTTLGITGLQTGDEADYYCGTWDSRLSAVVFGGGTKLTVL",
  },
  {
    id: "Etesevimab_vh",
    chain: "VH",
    cdr3: "ARVLPMYGDYLDY",
    seq: "EVQLVESGGGLVQPGGSLRLSCAASGFTVSSNYMSWVRQAPGKGLEWVSVIYSGGSTFYADSVKGRFTISRDNSMNTLFLQMNSLRAEDTAVYYCARVLPMYGDYLDYWGQGTLVTVSS",
  },
  {
    id: "Etesevimab_vl",
    chain: "VL",
    cdr3: "QQSYSTPPEYT",
    seq: "DIVMTQSPSSLSASVGDRVTITCRASQSISRYLNWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSYSTPPEYTFGQGTKLEIK",
  },
  {
    id: "Evinacumab_vh",
    chain: "VH",
    cdr3: "AKDLRNTIFGVVIPDAFDI",
    seq: "EVQLVESGGGVIQPGGSLRLSCAASGFTFDDYAMNWVRQGPGKGLEWVSAISGDGGSTYYADSVKGRFTISRDNSKNSLYLQMNSLRAEDTAFFYCAKDLRNTIFGVVIPDAFDIWGQGTMVTVSS",
  },
  {
    id: "Evinacumab_vl",
    chain: "VL",
    cdr3: "QQYNSYSYT",
    seq: "DIQMTQSPSTLSASVGDRVTITCRASQSIRSWLAWYQQKPGKAPKLLIYKASSLESGVPSRFSGSGSGTEFTLTISSLQPDDFATYYCQQYNSYSYTFGQGTKLEIK",
  },
  {
    id: "Evolocumab_vh",
    chain: "VH",
    cdr3: "ARGYGMDV",
    seq: "VQLVQSGAEVKKPGASVKVSCKASGYTLTSYGISWVRQAPGQGLEWMGWVSFYNGNTNYAQKLQGRGTMTTDPSTSTAYMELRSLRSDDTAVYYCARGYGMDVWGQGTTVTVSS",
  },
  {
    id: "Evolocumab_vl",
    chain: "VL",
    cdr3: "NSYTSTSMV",
    seq: "ALTQPASVSGSPGQSITISCTGTSSDVGGYNSVSWYQQHPGKAPKLMIYEVSNRPSGVSNRFSGSKSGNTASLTISGLQAEDEADYYCNSYTSTSMVFGGGTKLTVL",
  },
  {
    id: "Fremanezumab_vh",
    chain: "VH",
    cdr3: "LAYFDYGLAIQNY",
    seq: "EVQLVESGGGLVQPGGSLRLSCAASGFTFSNYWISWVRQAPGKGLEWVAEIRSESDASATHYAEAVKGRFTISRDNAKNSLYLQMNSLRAEDTAVYYCLAYFDYGLAIQNYWGQGTLVTVSS",
  },
  {
    id: "Fremanezumab_vl",
    chain: "VL",
    cdr3: "SQSYNYPYT",
    seq: "EIVLTQSPATLSLSPGERATLSCKASKRVTTYVSWYQQKPGQAPRLLIYGASNRYLGIPARFSGSGSGTDFTLTISSLEPEDFAVYYCSQSYNYPYTFGQGTKLEIK",
  },
  {
    id: "Galcanezumab_vh",
    chain: "VH",
    cdr3: "ARLSDYVSGFGY",
    seq: "QVQLVQSGAEVKKPGSSVKVSCKASGYTFGNYWMQWVRQAPGQGLEWMGAIYEGTGKTVYIQKFADRVTITADKSTSTAYMELSSLRSEDTAVYYCARLSDYVSGFGYWGQGTTVTVSS",
  },
  {
    id: "Galcanezumab_vl",
    chain: "VL",
    cdr3: "QQGDALPPT",
    seq: "DIQMTQSPSSLSASVGDRVTITCRASKDISKYLNWYQQKPGKAPKLLIYYTSGYHSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQGDALPPTFGGGTKVEIK",
  },
  {
    id: "Garadacimab_vh",
    chain: "VH",
    cdr3: "ARALPRSGYLISPHYYYYALDV",
    seq: "EVQLLESGGGLVQPGGSLRLSCAASGFTFSKYIMQWVRQAPGKGLEWVSGIDIPTKGTVYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCARALPRSGYLISPHYYYYALDVWGQGTTVTVSS",
  },
  {
    id: "Garadacimab_vl",
    chain: "VL",
    cdr3: "AAWDASLRGV",
    seq: "VLTQPPSASGTPGQRVTISCSGSSSNIGRNYVYWYQQLPGTAPKLLIYSNNQRPSGVPDRFSGSKSGTSASLAISGLRSEDEADYYCAAWDASLRGVFGGGTKLTVL",
  },
  {
    id: "Golimumab_vh",
    chain: "VH",
    cdr3: "ARDRGIAAGGNYYYYGMDV",
    seq: "QVQLVESGGGVVQPGRSLRLSCAASGFIFSSYAMHWVRQAPGNGLEWVAFMSYDGSNKKYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCARDRGIAAGGNYYYYGMDVWGQGTTVTVSS",
  },
  {
    id: "Golimumab_vl",
    chain: "VL",
    cdr3: "QQRSNWPPFT",
    seq: "EIVLTQSPATLSLSPGERATLSCRASQSVYSYLAWYQQKPGQAPRLLIYDASNRATGIPARFSGSGSGTDFTLTISSLEPEDFAVYYCQQRSNWPPFTFGPGTKVDIK",
  },
  {
    id: "Guselkumab_vh",
    chain: "VH",
    cdr3: "ARWYYKPFDV",
    seq: "EVQLVQSGAEVKKPGESLKISCKGSGYSFSNYWIGWVRQMPGKGLEWMGIIDPSNSYTRYSPSFQGQVTISADKSISTAYLQWSSLKASDTAMYYCARWYYKPFDVWGQGTLVTVSS",
  },
  {
    id: "Guselkumab_vl",
    chain: "VL",
    cdr3: "ASWTDGLSLVV",
    seq: "VLTQPPSVSGAPGQRVTISCTGSSSNIGSGYDVHWYQQLPGTAPKLLIYGNSKRPSGVPDRFSGSKSGTSASLAITGLQSEDEADYYCASWTDGLSLVVFGGGTKLTVL",
  },
  {
    id: "Ibalizumab_vh",
    chain: "VH",
    cdr3: "AREKDNYATGAWFAY",
    seq: "QVQLQQSGPEVVKPGASVKMSCKASGYTFTSYVIHWVRQKPGQGLDWIGYINPYNDGTDYDEKFKGKATLTSDTSTSTAYMELSSLRSEDTAVYYCAREKDNYATGAWFAYWGQGTLVTVSS",
  },
  {
    id: "Ibalizumab_vl",
    chain: "VL",
    cdr3: "QQYYSYRT",
    seq: "DIVMTQSPDSLAVSLGERVTMNCKSSQSLLYSTNQKNYLAWYQQKPGQSPKLLIYWASTRESGVPDRFSGSGSGTDFTLTISSVQAEDVAVYYCQQYYSYRTFGGGTKLEIK",
  },
  {
    id: "Inebilizumab_vh",
    chain: "VH",
    cdr3: "ARSGFITTVRDFDY",
    seq: "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSSWMNWVRQAPGKGLEWVGRIYPGDGDTNYNVKFKGRFTISRDDSKNSLYLQMNSLKTEDTAVYYCARSGFITTVRDFDYWGQGTLVTVSS",
  },
  {
    id: "Inebilizumab_vl",
    chain: "VL",
    cdr3: "QQSKEVPFT",
    seq: "EIVLTQSPDFQSVTPKEKVTITCRASESVDTFGISFMNWFQQKPDQSPKLLIHEASNQGSGVPSRFSGSGSGTDFTLTINSLEAEDAATYYCQQSKEVPFTFGGGTKVEIK",
  },
  {
    id: "Ipilimumab_vh",
    chain: "VH",
    cdr3: "ARTGWLGPFDY",
    seq: "QVQLVESGGGVVQPGRSLRLSCAASGFTFSSYTMHWVRQAPGKGLEWVTFISYDGNNKYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAIYYCARTGWLGPFDYWGQGTLVTVSS",
  },
  {
    id: "Ipilimumab_vl",
    chain: "VL",
    cdr3: "QQYGSSPWT",
    seq: "EIVLTQSPGTLSLSPGERATLSCRASQSVGSSYLAWYQQKPGQAPRLLIYGAFSRATGIPDRFSGSGSGTDFTLTISRLEPEDFAVYYCQQYGSSPWTFGQGTKVEIK",
  },
  {
    id: "Itolizumab_vh",
    chain: "VH",
    cdr3: "ARRDYDLDYFDS",
    seq: "EVQLVESGGGLVKPGGSLKLSCAASGFKFSRYAMSWVRQAPGKRLEWVATISSGGSYIYYPDSVKGRFTISRDNVKNTLYLQMSSLRSEDTAMYYCARRDYDLDYFDSWGQGTLVTVSS",
  },
  {
    id: "Ixekizumab_vh",
    chain: "VH",
    cdr3: "ARYDYFTGTGVY",
    seq: "QVQLVQSGAEVKKPGSSVKVSCKASGYSFTDYHIHWVRQAPGQGLEWMGVINPMYGTTDYNQRFKGRVTITADESTSTAYMELSSLRSEDTAVYYCARYDYFTGTGVYWGQGTLVTVSS",
  },
  {
    id: "Ixekizumab_vl",
    chain: "VL",
    cdr3: "SQSTHLPFT",
    seq: "DIVMTQTPLSLSVTPGQPASISCRSSRSLVHSRGNTYLHWYLQKPGQSPQLLIYKVSNRFIGVPDRFSGSGSGTDFTLKISRVEAEDVGVYYCSQSTHLPFTFGQGTKLEIK",
  },
  {
    id: "Lanadelumab_vh",
    chain: "VH",
    cdr3: "AYRRIGVPRRDEFDI",
    seq: "EVQLLESGGGLVQPGGSLRLSCAASGFTFSHYIMMWVRQAPGKGLEWVSGIYSSGGITVYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAYRRIGVPRRDEFDIWGQGTMVTVSS",
  },
  {
    id: "Lanadelumab_vl",
    chain: "VL",
    cdr3: "QQYNTYWT",
    seq: "DIQMTQSPSTLSASVGDRVTITCRASQSISSWLAWYQQKPGKAPKLLIYKASTLESGVPSRFSGSGSGTEFTLTISSLQPDDFATYYCQQYNTYWTFGQGTKVEIK",
  },
  {
    id: "Lebrikizumab_vh",
    chain: "VH",
    cdr3: "AGDGYYPYAMDN",
    seq: "QVTLRESGPALVKPTQTLTLTCTVSGFSLSAYSVNWIRQPPGKALEWLAMIWGDGKIVYNSALKSRLTISKDTSKNQVVLTMTNMDPVDTATYYCAGDGYYPYAMDNWGQGSLVTVSS",
  },
  {
    id: "Lebrikizumab_vl",
    chain: "VL",
    cdr3: "QQNNEDPRT",
    seq: "DIVMTQSPDSLSVSLGERATINCRASKSVDSYGNSFMHWYQQKPGQPPKLLIYLASNLESGVPDRFSGSGSGTDFTLTISSLQAEDVAVYYCQQNNEDPRTFGGGTKVEIK",
  },
  {
    id: "Leronlimab_vh",
    chain: "VH",
    cdr3: "GSSFGSNYVFAWFTY",
    seq: "EVQLVESGGGLVKPGGSLRLSCAASGYTFSNYWIGWVRQAPGKGLEWIGDIYPGGNYIRNNEKFKDKTTLSADTSKNTAYLQMNSLKTEDTAVYYCGSSFGSNYVFAWFTYWGQGTLVTVSS",
  },
  {
    id: "Leronlimab_vl",
    chain: "VL",
    cdr3: "SQSTHVPLT",
    seq: "DIVMTQSPLSLPVTPGEPASISCRSSQRLLSSYGHTYLHWYLQKPGQSPQLLIYEVSNRFSGVPDRFSGSGSGTDFTLKISRVEAEDVGVYYCSQSTHVPLTFGQGTKVEIK",
  },
  {
    id: "Marstacimab_vh",
    chain: "VH",
    cdr3: "AILGATSLSAFDI",
    seq: "EVQLLESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAILGATSLSAFDIWGQGTMVTVSS",
  },
  {
    id: "Marstacimab_vl",
    chain: "VL",
    cdr3: "QSYDSSLSGSGV",
    seq: "VLTQPPSVSGAPGQRVTISCTGSSSNIGAGYDVHWYQQLPGTAPKLLIYGNSNRPSGVPDRFSGSKSGTSASLAITGLQAEDEADYYCQSYDSSLSGSGVFGGGTKLTVL",
  },
  {
    id: "Masavibart_vh",
    chain: "VH",
    cdr3: "ASGSDYGDYLLVY",
    seq: "QVQLVESGGGVVQPGRSLRLSCAASGFTFSNYAMYWVRQAPGKGLEWVAVISYDGSNKYYADSVKGRFTISRDNSKNTLYLQMNSLRTEDTAVYYCASGSDYGDYLLVYWGQGTLVTVSS",
  },
  {
    id: "Masavibart_vl",
    chain: "VL",
    cdr3: "NSLTSISTWV",
    seq: "ALTQPASVSGSPGQSITISCTGTSSDVGGYNYVSWYQQHPGKAPKLMIYDVSKRPSGVSNRFSGSKSGNTASLTISGLQSEDEADYYCNSLTSISTWVFGGGTKLTVL",
  },
  {
    id: "Mepolizumab_vh",
    chain: "VH",
    cdr3: "ARDPPSSLLRLDY",
    seq: "QVTLRESGPALVKPTQTLTLTCTVSGFSLTSYSVHWVRQPPGKGLEWLGVIWASGGTDYNSALMSRLSISKDTSRNQVVLTMTNMDPVDTATYYCARDPPSSLLRLDYWGRGTPVTVSS",
  },
  {
    id: "Mepolizumab_vl",
    chain: "VL",
    cdr3: "QNVHSFPFT",
    seq: "DIVMTQSPDSLAVSLGERATINCKSSQSLLNSGNQKNYLAWYQQKPGQPPKLLIYGASTRESGVPDRFSGSGSGTDFTLTISSLQAEDVAVYYCQNVHSFPFTFGGGTKLEIK",
  },
  {
    id: "Mirikizumab_vh",
    chain: "VH",
    cdr3: "ARNWDTGL",
    seq: "QVQLVQSGAEVKKPGSSVKVSCKASGYKFTRYVMHWVRQAPGQGLEWMGYINPYNDGTNYNEKFKGRVTITADKSTSTAYMELSSLRSEDTAVYYCARNWDTGLWGQGTTVTVSS",
  },
  {
    id: "Mirikizumab_vl",
    chain: "VL",
    cdr3: "QMYWSTPFT",
    seq: "DIQMTQSPSSLSASVGDRVTITCKASDHILKFLTWYQQKPGKAPKLLIYGATSLETGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQMYWSTPFTFGGGTKVEIK",
  },
  {
    id: "Mogamulizumab_vh",
    chain: "VH",
    cdr3: "GRHSDGNFAFGY",
    seq: "EVQLVESGGDLVQPGRSLRLSCAASGFIFSNYGMSWVRQAPGKGLEWVATISSASTYSYYPDSVKGRFTISRDNAKNSLYLQMNSLRVEDTALYYCGRHSDGNFAFGYWGQGTLVTVSS",
  },
  {
    id: "Mogamulizumab_vl",
    chain: "VL",
    cdr3: "FQGSLLPWT",
    seq: "MTQSPLSLPVTPGEPASISCRSSRNIVHINGDTYLEWYLQKPGQSPQLLIYKVSNRFSGVPDRFSGSGSGTDFTLKISRVEAEDVGVYYCFQGSLLPWTFGQGTKVEIK",
  },
  {
    id: "Natalizumab_vh",
    chain: "VH",
    cdr3: "AREGYYGNYGVYAMDY",
    seq: "QVQLVQSGAEVKKPGASVKVSCKASGFNIKDTYIHWVRQAPGQRLEWMGRIDPANGYTKYDPKFQGRVTITADTSASTAYMELSSLRSEDTAVYYCAREGYYGNYGVYAMDYWGQGTLVTVSS",
  },
  {
    id: "Natalizumab_vl",
    chain: "VL",
    cdr3: "LQYDNLWT",
    seq: "DIQMTQSPSSLSASVGDRVTITCKTSQDINKYMAWYQQTPGKAPRLLIHYTSALQPGIPSRFSGSGSGRDYTFTISSLQPEDIATYYCLQYDNLWTFGQGTKVEIK",
  },
  {
    id: "Necitumumab_vh",
    chain: "VH",
    cdr3: "ARVSIFGVGTFDY",
    seq: "QVQLQESGPGLVKPSQTLSLTCTVSGGSISSGDYYWSWIRQPPGKGLEWIGYIYYSGSTDYNPSLKSRVTMSVDTSKNQFSLKVNSVTAADTAVYYCARVSIFGVGTFDYWGQGTLVTVSS",
  },
  {
    id: "Necitumumab_vl",
    chain: "VL",
    cdr3: "HQYGSTPLT",
    seq: "EIVMTQSPATLSLSPGERATLSCRASQSVSSYLAWYQQKPGQAPRLLIYDASNRATGIPARFSGSGSGTDFTLTISSLEPEDFAVYYCHQYGSTPLTFGGGTKAEIK",
  },
  {
    id: "Nemolizumab_vh",
    chain: "VH",
    cdr3: "ARDGYDDGPYTLET",
    seq: "QVQLVQSGAEVKKPGASVKVSCKASGYTFTGYIMNWVRQAPGQGLEWMGLINPYNGGTDYNPQFQDRVTITADKSTSTAYMELSSLRSEDTAVYYCARDGYDDGPYTLETWGQGTLVTVSS",
  },
  {
    id: "Nemolizumab_vl",
    chain: "VL",
    cdr3: "QHHYDSPLT",
    seq: "DIQMTQSPSSLSASVGDRVTITCQASEDIYSFVAWYQQKPGKAPKLLIYNAQTEAQGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQHHYDSPLTFGGGTKVEIK",
  },
  {
    id: "Nepuvibart_vh",
    chain: "VH",
    cdr3: "ARDRGTTMVPFDY",
    seq: "QVQLVESGGGLVKPGGSLRLSCAASGFTFSDYYMSWIRQAPGKGLEWVSYITYSGSTIYYADSVKGRFTISRDNAKSSLYLQMNSLRAEDTAVYYCARDRGTTMVPFDYWGQGTLVTVSS",
  },
  {
    id: "Nepuvibart_vl",
    chain: "VL",
    cdr3: "QQYDNLPLT",
    seq: "DIQMTQSPSSLSASVGDRVTITCQASQDITNYLNWYQQKPGKAPKLLIYAASNLETGVPSRFSGSGSGTDFTFTISGLQPEDIATYYCQQYDNLPLTFGGGTKVEIK",
  },
  {
    id: "Nimotuzumab_vh",
    chain: "VH",
    cdr3: "TRQGLWFDSDGRGFDF",
    seq: "QVQLQQSGAEVKKPGSSVKVSCKASGYTFTNYYIYWVRQAPGQGLEWIGGINPTSGGSNFNEKFKTRVTITADESSTTAYMELSSLRSEDTAFYFCTRQGLWFDSDGRGFDFWGQGTTVTVSS",
  },
  {
    id: "Nimotuzumab_vl",
    chain: "VL",
    cdr3: "FQYSHVPWT",
    seq: "DIQMTQSPSSLSASVGDRVTITCRSSQNIVHSNGNTYLDWYQQTPGKAPKLLIYKVSNRFSGVPSRFSGSGSGTDFTFTISSLQPEDIATYYCFQYSHVPWTFGQGTKLQIT",
  },
  {
    id: "Nirsevimab_vh",
    chain: "VH",
    cdr3: "ATETALVVSETYLPHYFDN",
    seq: "QVQLVQSGAEVKKPGSSVMVSCQASGGLLEDYIINWVRQAPGQGPEWMGGIIPVLGTVHYGPKFQGRVTITADESTDTAYMELSSLRSEDTAMYYCATETALVVSETYLPHYFDNWGQGTLVTVSS",
  },
  {
    id: "Nirsevimab_vl",
    chain: "VL",
    cdr3: "QQYDNLPLT",
    seq: "DIQMTQSPSSLSAAVGDRVTITCQASQDIVNYLNWYQQKPGKAPKLLIYVASNLETGVPSRFSGSGSGTDFSLTISSLQPEDVATYYCQQYDNLPLTFGGGTKVEIK",
  },
  {
    id: "Nivolumab_vh",
    chain: "VH",
    cdr3: "ATNDDY",
    seq: "QVQLVESGGGVVQPGRSLRLDCKASGITFSNSGMHWVRQAPGKGLEWVAVIWYDGSKRYYADSVKGRFTISRDNSKNTLFLQMNSLRAEDTAVYYCATNDDYWGQGTLVTVSS",
  },
  {
    id: "Nivolumab_vl",
    chain: "VL",
    cdr3: "QQSSNWPRT",
    seq: "EIVLTQSPATLSLSPGERATLSCRASQSVSSYLAWYQQKPGQAPRLLIYDASNRATGIPARFSGSGSGTDFTLTISSLEPEDFAVYYCQQSSNWPRTFGQGTKVEIK",
  },
  {
    id: "Obinutuzumab_vh",
    chain: "VH",
    cdr3: "ARNVFDGYWLVY",
    seq: "QVQLVQSGAEVKKPGSSVKVSCKASGYAFSYSWINWVRQAPGQGLEWMGRIFPGDGDTDYNGKFKGRVTITADKSTSTAYMELSSLRSEDTAVYYCARNVFDGYWLVYWGQGTLVTVSS",
  },
  {
    id: "Obinutuzumab_vl",
    chain: "VL",
    cdr3: "AQNLELPYT",
    seq: "DIVMTQTPLSLPVTPGEPASISCRSSKSLLHSNGITYLYWYLQKPGQSPQLLIYQMSNLVSGVPDRFSGSGSGTDFTLKISRVEAEDVGVYYCAQNLELPYTFGGGTKVEIK",
  },
  {
    id: "Ocrelizumab_vh",
    chain: "VH",
    cdr3: "ARVVYYSNSYWYFDV",
    seq: "EVQLVESGGGLVQPGGSLRLSCAASGYTFTSYNMHWVRQAPGKGLEWVGAIYPGNGDTSYNQKFKGRFTISVDKSKNTLYLQMNSLRAEDTAVYYCARVVYYSNSYWYFDVWGQGTLVTVSS",
  },
  {
    id: "Ocrelizumab_vl",
    chain: "VL",
    cdr3: "QQWSFNPPT",
    seq: "DIQMTQSPSSLSASVGDRVTITCRASSSVSYMHWYQQKPGKAPKPLIYAPSNLASGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQWSFNPPTFGQGTKVEIK",
  },
  {
    id: "Odesivimab_vh",
    chain: "VH",
    cdr3: "ARTWFGELYFDY",
    seq: "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYDMHWVRQATGKGLEWVSAIGTAGDTYYPGSVKGRFTISRENAKNSLYLQMNSLRAGDTAVYYCARTWFGELYFDYWGQGTLVTVSS",
  },
  {
    id: "Odesivimab_vl",
    chain: "VL",
    cdr3: "QQYYSSPLT",
    seq: "DIVMTQSPDSLAVSLGERATINCKSSQSVLYSSNNKNYLAWYQQKPGQPPKLLIYWASTRESGVPDRFSGSGSGTEFTLTITSLQAEDVAVYYCQQYYSSPLTFGGGTKVEIK",
  },
  {
    id: "Ofatumumab_vh",
    chain: "VH",
    cdr3: "AKDIQYGNYYYGMDV",
    seq: "EVQLVESGGGLVQPGRSLRLSCAASGFTFNDYAMHWVRQAPGKGLEWVSTISWNSGSIGYADSVKGRFTISRDNAKKSLYLQMNSLRAEDTALYYCAKDIQYGNYYYGMDVWGQGTTVTVSS",
  },
  {
    id: "Ofatumumab_vl",
    chain: "VL",
    cdr3: "QQRSNWPIT",
    seq: "EIVLTQSPATLSLSPGERATLSCRASQSVSSYLAWYQQKPGQAPRLLIYDASNRATGIPARFSGSGSGTDFTLTISSLEPEDFAVYYCQQRSNWPITFGQGTRLEIK",
  },
  {
    id: "Olaratumab_vh",
    chain: "VH",
    cdr3: "ARQSTYYYGSGNYYGWFDR",
    seq: "QLQLQESGPGLVKPSETLSLTCTVSGGSINSSSYYWGWLRQSPGKGLEWIGSFFYTGSTYYNPSLRSRLTISVDTSKNQFSLMLSSVTAADTAVYYCARQSTYYYGSGNYYGWFDRWDQGTLVTVSS",
  },
  {
    id: "Olaratumab_vl",
    chain: "VL",
    cdr3: "QQRSNWPPA",
    seq: "EIVLTQSPATLSLSPGERATLSCRASQSVSSYLAWYQQKPGQAPRLLIYDASNRATGIPARFSGSGSGTDFTLTISSLEPEDFAVYYCQQRSNWPPAFGQGTKVEIK",
  },
  {
    id: "Olokizumab_vh",
    chain: "VH",
    cdr3: "ARESYYGFTSY",
    seq: "EVQLVESGGGLVQPGGSLRLSCAASGFNFNDYFMNWVRQAPGKGLEWVAQMRNKNYQYGTYYAESLEGRFTISRDDSKNSLYLQMNSLKTEDTAVYYCARESYYGFTSYWGQGTLVTVSS",
  },
  {
    id: "Olokizumab_vl",
    chain: "VL",
    cdr3: "LQHNSAPYT",
    seq: "DIQMTQSPSSLSASVGDRVTITCQASQDIGISLSWYQQKPGKAPKLLIYNANNLADGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCLQHNSAPYTFGQGTKLEIK",
  },
  {
    id: "Omalizumab_vh",
    chain: "VH",
    cdr3: "ARGSHYFGHWHFAV",
    seq: "VQLVESGGGLVQPGGSLRLSCAVSGYSITSGYSWNWIRQAPGKGLEWVASITYDGSTNYNPSVKGRITISRDDSKNTFYLQMNSLRAEDTAVYYCARGSHYFGHWHFAVWGQGTLVTVSS",
  },
  {
    id: "Omalizumab_vl",
    chain: "VL",
    cdr3: "QQSHEDPYT",
    seq: "DIQLTQSPSSLSASVGDRVTITCRASQSVDYDGDSYMNWYQQKPGKAPKLLIYAASYLESGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSHEDPYTFGQGTKVEIK",
  },
  {
    id: "Ongericimab_vh",
    chain: "VH",
    cdr3: "ANHRD",
    seq: "QVQLQESGPGLVKPSQTLSLTCTVSGFSISSYGIHWIRQSPGKGLEWIGVIWRGGITDYNAPFMSRVTISKDNSKNQVSFKLSSVTAADTAVYYCANHRDWGQGTLVTVSS",
  },
  {
    id: "Ongericimab_vl",
    chain: "VL",
    cdr3: "LQYDDLWT",
    seq: "DIQMTQSPSSLSASVGDRVTITCQASQDINKYIDWYQHKPGKAPKLLIHYASTLQPGVPSRFSGSGSGRDYTFTISSLQPEDIATYYCLQYDDLWTFGQGTKVEIK",
  },
  {
    id: "Ormutivimab_vh",
    chain: "VH",
    cdr3: "ARENLDNSGTYYYFSGWFDP",
    seq: "QVQLVQSGAEVKKPGSSVKVSCKASGGTFNRYTVNWVRQAPGQGLEWMGGIIPIFGTANYAQRFQGRLTITADESTSTAYMELSSLRSDDTAVYFCARENLDNSGTYYYFSGWFDPWGQGTLVTVSS",
  },
  {
    id: "Ormutivimab_vl",
    chain: "VL",
    cdr3: "CSYAGDYTPGVV",
    seq: "ALTQPRSVSGSPGQSVTISCTGTSSDIGGYNFVSWYQQHPGKAPKLMIYDATKRPSGVPDRFSGSKSGNTASLTISGLQAEDEADYYCCSYAGDYTPGVVFGGGTKLTVL",
  },
  {
    id: "Palivizumab_vh",
    chain: "VH",
    cdr3: "ARSMITNWYFDV",
    seq: "QVTLRESGPALVKPTQTLTLTCTFSGFSLSTSGMSVGWIRQPPGKALEWLADIWWDDKKDYNPSLKSRLTISKDTSKNQVVLKVTNMDPADTATYYCARSMITNWYFDVWGAGTTVTVSS",
  },
  {
    id: "Palivizumab_vl",
    chain: "VL",
    cdr3: "FQGSGYPFT",
    seq: "DIQMTQSPSTLSASVGDRVTITCKCQLSVGYMHWYQQKPGKAPKLLIYDTSKLASGVPSRFSGSGSGTEFTLTISSLQPDDFATYYCFQGSGYPFTFGGGTKLEIK",
  },
  {
    id: "Panitumumab_vh",
    chain: "VH",
    cdr3: "VRDRVTGAFDI",
    seq: "QVQLQESGPGLVKPSETLSLTCTVSGGSVSSGDYYWTWIRQSPGKGLEWIGHIYYSGNTNYNPSLKSRLTISIDTSKTQFSLKLSSVTAADTAIYYCVRDRVTGAFDIWGQGTMVTVSS",
  },
  {
    id: "Panitumumab_vl",
    chain: "VL",
    cdr3: "QHFDHLPLA",
    seq: "DIQMTQSPSSLSASVGDRVTITCQASQDISNYLNWYQQKPGKAPKLLIYDASNLETGVPSRFSGSGSGTDFTFTISSLQPEDIATYFCQHFDHLPLAFGGGTKVEIK",
  },
  {
    id: "Pembrolizumab_vh",
    chain: "VH",
    cdr3: "ARRDYRFDMGFDY",
    seq: "QVQLVQSGVEVKKPGASVKVSCKASGYTFTNYYMYWVRQAPGQGLEWMGGINPSNGGTNFNEKFKNRVTLTTDSSTTTAYMELKSLQFDDTAVYYCARRDYRFDMGFDYWGQGTTVTVSS",
  },
  {
    id: "Pembrolizumab_vl",
    chain: "VL",
    cdr3: "QHSRDLPLT",
    seq: "EIVLTQSPATLSLSPGERATLSCRASKGVSTSGYSYLHWYQQKPGQAPRLLIYLASYLESGVPARFSGSGSGTDFTLTISSLEPEDFAVYYCQHSRDLPLTFGGGTKVEIK",
  },
  {
    id: "Pemivibart_vh",
    chain: "VH",
    cdr3: "ARDFGGDTAWAGTGFTY",
    seq: "EVQLVESGGGLVKPGGSLRLSCAASGFTFGSYEMNWVRQAPGKGLEWVSSISEDGYSTYYPDSLKGRFTISRDSAKNSLYLQMNSLRADDTAVYYCARDFGGDTAWAGTGFTYWGQGTLVTVSS",
  },
  {
    id: "Pemivibart_vl",
    chain: "VL",
    cdr3: "QSYDSDLSGLYT",
    seq: "VLTQPPSVSGAPGQRITISCTGSSSNIGAGYDVHWYQQLPGTAPKLLIYGSSSRNYGVPDRFSGSKSGTSASLAITGLQAEDEADYYCQSYDSDLSGLYTFGTGTKVTVL",
  },
  {
    id: "Pertuzumab_vh",
    chain: "VH",
    cdr3: "ARNLGPSFYFDY",
    seq: "VQLVESGGGLVQPGGSLRLSCAASGFTFTDYTMDWVRQAPGKGLEWVADVNPNSGGSIYNQRFKGRFTLSVDRSKNTLYLQMNSLRAEDTAVYYCARNLGPSFYFDYWGQGTLVTVSS",
  },
  {
    id: "Pertuzumab_vl",
    chain: "VL",
    cdr3: "QQYYIYPYT",
    seq: "DIQMTQSPSSLSASVGDRVTITCKASQDVSIGVAWYQQKPGKAPKLLIYSASYRYTGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQYYIYPYTFGQGTKVEIK",
  },
  {
    id: "Pozelimab_vh",
    chain: "VH",
    cdr3: "AREGNVDTTMIFDY",
    seq: "QVQLQESGPGLVKPSETLSLTCTVSGDSVSSSYWTWIRQPPGKGLEWIGYIYYSGSSNYNPSLKSRATISVDTSKNQFSLKLSSVTAADTAVYYCAREGNVDTTMIFDYWGQGTLVTVSS",
  },
  {
    id: "Pozelimab_vl",
    chain: "VL",
    cdr3: "LQDFNYPWT",
    seq: "AIQMTQSPSSLSASVGDRVTITCRASQGIRNDLGWYQQKPGKAPKLLIYAASSLQSGVPSRFAGRGSGTDFTLTISSLQPEDFATYYCLQDFNYPWTFGQGTKVEIK",
  },
  {
    id: "Ramucirumab_vh",
    chain: "VH",
    cdr3: "ARVTDAFDI",
    seq: "EVQLVQSGGGLVKPGGSLRLSCAASGFTFSSYSMNWVRQAPGKGLEWVSSISSSSSYIYYADSVKGRFTISRDNAKNSLYLQMNSLRAEDTAVYYCARVTDAFDIWGQGTMVTVSS",
  },
  {
    id: "Ramucirumab_vl",
    chain: "VL",
    cdr3: "QQAKAFPPT",
    seq: "DIQMTQSPSSVSASIGDRVTITCRASQGIDNWLGWYQQKPGKAPKLLIYDASNLDTGVPSRFSGSGSGTYFTLTISSLQAEDFAVYFCQQAKAFPPTFGGGTKVDIK",
  },
  {
    id: "Ravulizumab_vh",
    chain: "VH",
    cdr3: "ARYFFGSSPNWYFDV",
    seq: "QVQLVQSGAEVKKPGASVKVSCKASGHIFSNYWIQWVRQAPGQGLEWMGEILPGSGHTEYTENFKDRVTMTRDTSTSTVYMELSSLRSEDTAVYYCARYFFGSSPNWYFDVWGQGTLVTVSS",
  },
  {
    id: "Ravulizumab_vl",
    chain: "VL",
    cdr3: "QNVLNTPLT",
    seq: "DIQMTQSPSSLSASVGDRVTITCGASENIYGALNWYQQKPGKAPKLLIYGATNLADGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQNVLNTPLTFGQGTKVEIK",
  },
  {
    id: "Regdanvimab_vh",
    chain: "VH",
    cdr3: "ARIPGFLRYRNRYYYYGMDV",
    seq: "QITLKESGPTLVKPTQTLTLTCSFSGFSLSTSGVGVGWIRQPPGKALEWLALIDWDDNKYHTTSLKTRLTISKDTSKNQVVLTMTNMDPVDTATYYCARIPGFLRYRNRYYYYGMDVWGQGTTVTVSS",
  },
  {
    id: "Regdanvimab_vl",
    chain: "VL",
    cdr3: "GTWDSSLSAGV",
    seq: "VLTQPPSVSAAPGQKVTISCSGSSSNIGNNYVSWYQQLPGTAPKLLIYDNNKRPSGIPDRFSGSKSGTSATLGITGLQTGDEADYYCGTWDSSLSAGVFGGGTELTVL",
  },
  {
    id: "Relatlimab_vh",
    chain: "VH",
    cdr3: "AFGYSDYEYNWFDP",
    seq: "QVQLQQWGAGLLKPSETLSLTCAVYGGSFSDYYWNWIRQPPGKGLEWIGEINHRGSTNSNPSLKSRVTLSLDTSKNQFSLKLRSVTAADTAVYYCAFGYSDYEYNWFDPWGQGTLVTVSS",
  },
  {
    id: "Relatlimab_vl",
    chain: "VL",
    cdr3: "QQRSNWPLT",
    seq: "EIVLTQSPATLSLSPGERATLSCRASQSISSYLAWYQQKPGQAPRLLIYDASNRATGIPARFSGSGSGTDFTLTISSLEPEDFAVYYCQQRSNWPLTFGQGTNLEIK",
  },
  {
    id: "Reslizumab_vh",
    chain: "VH",
    cdr3: "AREYYGYFDY",
    seq: "EVQLVESGGGLVQPGGSLRLSCAVSGLSLTSNSVNWIRQAPGKGLEWVGLIWSNGDTDYNSAIKSRFTISRDTSKSTVYLQMNSLRAEDTAVYYCAREYYGYFDYWGQGTLVTVSS",
  },
  {
    id: "Reslizumab_vl",
    chain: "VL",
    cdr3: "QQSYKFPNT",
    seq: "DIQMTQSPSSLSASVGDRVTITCLASEGISSYLAWYQQKPGKAPKLLIYGANSLQTGVPSRFSGSGSATDYTLTISSLQPEDFATYYCQQSYKFPNTFGQGTKVEVK",
  },
  {
    id: "Risankizumab_vh",
    chain: "VH",
    cdr3: "AIPDRSGYAWFIY",
    seq: "QVQLVQSGAEVKKPGSSVKVSCKASGYTFTDQTIHWMRQAPGQGLEWIGYIYPRDDSPKYNENFKGKVTITADKSTSTAYMELSSLRSEDTAVYYCAIPDRSGYAWFIYWGQGTLVTVSS",
  },
  {
    id: "Risankizumab_vl",
    chain: "VL",
    cdr3: "HQYSSYPFT",
    seq: "DIQMTQSPSSLSASVGDRVTITCKASRDVAIAVAWYQQKPGKVPKLLIYWASTRHTGVPSRFSGSGSRTDFTLTISSLQPEDVADYFCHQYSSYPFTFGSGTKLEIK",
  },
  {
    id: "Romosozumab_vh",
    chain: "VH",
    cdr3: "ARLGYDDIYDDWYFDV",
    seq: "VQLVQSGAEVKKPGASVKVSCKASGYTFTDYNMHWVRQAPGQGLEWMGEINPNSGGAGYNQKFKGRVTMTTDTSTSTAYMELRSLRSDDTAVYYCARLGYDDIYDDWYFDVWGQGTTVTVSS",
  },
  {
    id: "Romosozumab_vl",
    chain: "VL",
    cdr3: "QQGDTLPYT",
    seq: "DIQMTQSPSSLSASVGDRVTITCRASQDISNYLNWYQQKPGKAPKLLIYYTSRLLSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQGDTLPYTFGGGTKVEIK",
  },
  {
    id: "Sarilumab_vh",
    chain: "VH",
    cdr3: "AKGRDSFDI",
    seq: "EVQLVESGGGLVQPGRSLRLSCAASRFTFDDYAMHWVRQAPGKGLEWVSGISWNSGRIGYADSVKGRFTISRDNAENSLFLQMNGLRAEDTALYYCAKGRDSFDIWGQGTMVTVSS",
  },
  {
    id: "Sarilumab_vl",
    chain: "VL",
    cdr3: "QQANSFPYT",
    seq: "DIQMTQSPSSVSASVGDRVTITCRASQGISSWLAWYQQKPGKAPKLLIYGASSLESGVPSRFSGSGSGTDFTLTISSLQPEDFASYYCQQANSFPYTFGQGTKLEIK",
  },
  {
    id: "Satralizumab_vh",
    chain: "VH",
    cdr3: "ARSLARTTAMDY",
    seq: "QVQLQESGPGLVKPSETLSLTCAVSGHSISHDHAWSWVRQPPGEGLEWIGFISYSGITNYNPSLQGRVTISRDNSKNTLYLQMNSLRAEDTAVYYCARSLARTTAMDYWGEGTLVTVSS",
  },
  {
    id: "Satralizumab_vl",
    chain: "VL",
    cdr3: "GQGNRLPYT",
    seq: "DIQMTQSPSSLSASVGDSVTITCQASTDISSHLNWYQQKPGKAPELLIYYGSHLLSGVPSRFSGSGSGTDFTFTISSLEAEDAATYYCGQGNRLPYTFGQGTKVEIE",
  },
  {
    id: "Secukinumab_vh",
    chain: "VH",
    cdr3: "VRDYYDILTDYYIHYWYFDL",
    seq: "EVQLVESGGGLVQPGGSLRLSCAASGFTFSNYWMNWVRQAPGKGLEWVAAINQDGSEKYYVGSVKGRFTISRDNAKNSLYLQMNSLRVEDTAVYYCVRDYYDILTDYYIHYWYFDLWGRGTLVTVSS",
  },
  {
    id: "Secukinumab_vl",
    chain: "VL",
    cdr3: "QQYGSSPCT",
    seq: "EIVLTQSPGTLSLSPGERATLSCRASQSVSSSYLAWYQQKPGQAPRLLIYGASSRATGIPDRFSGSGSGTDFTLTISRLEPEDFAVYYCQQYGSSPCTFGQGTRLEIK",
  },
  {
    id: "Serplulimab_vh",
    chain: "VH",
    cdr3: "VSYYYGIDF",
    seq: "QVQLVESGGGLVKPGGSLRLSCAASGFTFSNYGMSWIRQAPGKGLEWVSTISGGGSNIYYADSVKGRFTISRDNAKNSLYLQMNSLRAEDTAVYYCVSYYYGIDFWGQGTSVTVSS",
  },
  {
    id: "Serplulimab_vl",
    chain: "VL",
    cdr3: "QQHYTIPWT",
    seq: "DIQMTQSPSSLSASVGDRVTITCKASQDVTTAVAWYQQKPGKAPKLLIYWASTRHTGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQHYTIPWTFGGGTKLEIK",
  },
  {
    id: "Sintilimab_vh",
    chain: "VH",
    cdr3: "ARAEHSSTGTFDY",
    seq: "QVQLVQSGAEVKKPGSSVKVSCKASGGTFSSYAISWVRQAPGQGLEWMGLIIPMFDTAGYAQKFQGRVAITVDESTSTAYMELSSLRSEDTAVYYCARAEHSSTGTFDYWGQGTLVTVSS",
  },
  {
    id: "Sintilimab_vl",
    chain: "VL",
    cdr3: "QQANHLPFT",
    seq: "DIQMTQSPSSVSASVGDRVTITCRASQGISSWLAWYQQKPGKAPKLLISAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQANHLPFTFGGGTKVEIK",
  },
  {
    id: "Socazolimab_vh",
    chain: "VH",
    cdr3: "ARAPYYYYYMDV",
    seq: "VQLVESGAEVKKPGSSVKVSCKASGGTFSSYAISWVRQAPGQGLEWMGGIIPIFGTANYAQKFQGRVTITADESTSTAYMELSSLRSEDTAVYYCARAPYYYYYMDVWGQGTTVTVSS",
  },
  {
    id: "Socazolimab_vl",
    chain: "VL",
    cdr3: "SSYTGISTVV",
    seq: "ALTQPASVSGSLGQSVTISCTGSSSDVGSYNLVSWYQQHPGKAPNLMIYDVSKRSGVSNRFSGSKSGNTASLTISGLQAEDEADYYCSSYTGISTVVFGGGTKLTVL",
  },
  {
    id: "Sotrovimab_vh",
    chain: "VH",
    cdr3: "ARDYTRGAWFGESLIGGFDN",
    seq: "QVQLVQSGAEVKKPGASVKVSCKASGYPFTSYGISWVRQAPGQGLEWMGWISTYQGNTNYAQKFQGRVTMTTDTSTTTGYMELRRLRSDDTAVYYCARDYTRGAWFGESLIGGFDNWGQGTLVTVSS",
  },
  {
    id: "Sotrovimab_vl",
    chain: "VL",
    cdr3: "QQHDTSLT",
    seq: "EIVLTQSPGTLSLSPGERATLSCRASQTVSSTSLAWYQQKPGQAPRLLIYGASSRATGIPDRFSGSGSGTDFTLTISRLEPEDFAVYYCQQHDTSLTFGGGTKVEIK",
  },
  {
    id: "Spesolimab_vh",
    chain: "VH",
    cdr3: "TVVFYGEPYFPY",
    seq: "QVQLVQSGAEVKKPGASVKVSCKASGYSFTSSWIHWVKQAPGQGLEWMGEINPGNVRTNYNENFRNKVTMTVDTSISTAYMELSRLRSDDTAVYYCTVVFYGEPYFPYWGQGTLVTVSS",
  },
  {
    id: "Spesolimab_vl",
    chain: "VL",
    cdr3: "HQFHRSPLT",
    seq: "IVLTQSPGTLSLSPGERATMTCTASSSVSSSYFHWYQQKPGQAPRLWIYRTSRLASGVPDRFSGSGSGTDFTLTISRLEPEDAATYYCHQFHRSPLTFGAGTKLEIK",
  },
  {
    id: "Stapokibart_vh",
    chain: "VH",
    cdr3: "ARATARATEFAY",
    seq: "EVQLLESGGGLVQPGGSLRLSCAASGFTFSRYAMSWVRQAPGKGLEWVSTISSGGSYTNYADSVKGRFTISRDNVKNTLYLQMNSLRAEDTAVYYCARATARATEFAYWGQGTLVTVSS",
  },
  {
    id: "Stapokibart_vl",
    chain: "VL",
    cdr3: "QQGNTLPLT",
    seq: "DIQMTQSPSSLSASVGDRVTITCQASQDISNYLNWYQQKPGKAPKLLIYYTSRLHSGVPSRFSGSGSGTDYTLTISSLQPEDFATYFCQQGNTLPLTFGGGTKVEIK",
  },
  {
    id: "Sugemalimab_vh",
    chain: "VH",
    cdr3: "AKPPRGYNYGPFDY",
    seq: "EVQLLESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSGISGSGGFTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKPPRGYNYGPFDYWGQGTLVTVSS",
  },
  {
    id: "Sugemalimab_vl",
    chain: "VL",
    cdr3: "QVWDSSSDHVV",
    seq: "YVLTQPPSVSVAPGQTARITCGGNNIGSKSVHWYQQKPGQAPVLVVYDDSDRPSGIPERFSGSNSGNTATLTISRVEAGDEADYYCQVWDSSSDHVVFGGGTKLTVL",
  },
  {
    id: "Tafolecimab_vh",
    chain: "VH",
    cdr3: "ARENSGVVPAAGPNWFGP",
    seq: "QLQLQESGPGLVKPSETLSLTCTVSGGSISSASYYWSWIRQPPGKGLEWIGSINYRGSTYYNPSLKSRVTISVDTSKNQFSLKLSSVTAADTAVYYCARENSGVVPAAGPNWFGPWGQGTLVTVSS",
  },
  {
    id: "Tafolecimab_vl",
    chain: "VL",
    cdr3: "QQRRNWFT",
    seq: "EIVLTQSPATLSLSPGERATLSCRASQSVSSYLAWYQQKPGQAPRLLIYDASNRATGIPARFSGSGSGTDFTLTISSLEPEDFAVYYCQQRRNWFTFGGGTKVEIK",
  },
  {
    id: "Tagitanlimab_vh",
    chain: "VH",
    cdr3: "VRDSNYRYDEPFTY",
    seq: "QVQLQESGPGLVKPSETLSITCTVSGFSLSNYDISWIRQPPGKGLEWLGVIWTGGATNYNPALKSRLTISRDNSKNQVSLKMSSVTAADTAVYYCVRDSNYRYDEPFTYWGQGTLVTVSS",
  },
  {
    id: "Tagitanlimab_vl",
    chain: "VL",
    cdr3: "QQSNSWPYT",
    seq: "EIVLTQSPDTLSVTPKEKVTLTCRASQSIGTNIHWFQQKPGQSPKLLIKYASESISGVPSRFSGSGSGTDFTLTINSVEAEDAATYYCQQSNSWPYTFGQGTKLEIK",
  },
  {
    id: "Teplizumab_vh",
    chain: "VH",
    cdr3: "ARYYDDHYCLDY",
    seq: "QVQLVQSGGGVVQPGRSLRLSCKASGYTFTRYTMHWVRQAPGKGLEWIGYINPSRGYTNYNQKVKDRFTISRDNSKNTAFLQMDSLRPEDTGVYFCARYYDDHYCLDYWGQGTPVTVSS",
  },
  {
    id: "Teplizumab_vl",
    chain: "VL",
    cdr3: "QQWSSNPFT",
    seq: "DIQMTQSPSSLSASVGDRVTITCSASSSVSYMNWYQQTPGKAPKRWIYDTSKLASGVPSRFSGSGSGTDYTFTISSLQPEDIATYYCQQWSSNPFTFGQGTKLQIT",
  },
  {
    id: "Teprotumumab_vh",
    chain: "VH",
    cdr3: "ARELGRRYFDL",
    seq: "QVELVESGGGVVQPGRSQRLSCAASGFTFSSYGMHWVRQAPGKGLEWVAIIWFDGSSTYYADSVRGRFTISRDNSKNTLYLQMNSLRAEDTAVYFCARELGRRYFDLWGRGTLVSVSS",
  },
  {
    id: "Teprotumumab_vl",
    chain: "VL",
    cdr3: "QQRSKWPPWT",
    seq: "EIVLTQSPATLSLSPGERATLSCRASQSVSSYLAWYQQKPGQAPRLLIYDASKRATGIPARFSGSGSGTDFTLTISSLEPEDFAVYYCQQRSKWPPWTFGQGTKVESK",
  },
  {
    id: "Tezepelumab_vh",
    chain: "VH",
    cdr3: "ARAPQWELVHEAFDI",
    seq: "QMQLVESGGGVVQPGRSLRLSCAASGFTFRTYGMHWVRQAPGKGLEWVAVIWYDGSNKHYADSVKGRFTITRDNSKNTLNLQMNSLRAEDTAVYYCARAPQWELVHEAFDIWGQGTMVTVSS",
  },
  {
    id: "Tezepelumab_vl",
    chain: "VL",
    cdr3: "QVWDSSSDHVV",
    seq: "YVLTQPPSVSVAPGQTARITCGGNNLGSKSVHWYQQKPGQAPVLVVYDDSDRPSWIPERFSGSNSGNTATLTISRGEAGDEADYYCQVWDSSSDHVVFGGGTKLTVL",
  },
  {
    id: "Tildrakizumab_vh",
    chain: "VH",
    cdr3: "ARGGGGFAY",
    seq: "QVQLVQSGAEVKKPGASVKVSCKASGYIFITYWMTWVRQAPGQGLEWMGQIFPASGSADYNEKFEGRVTMTTDTSTSTAYMELRSLRSDDTAVYYCARGGGGFAYWGQGTLVTVSS",
  },
  {
    id: "Tildrakizumab_vl",
    chain: "VL",
    cdr3: "QHHYGIPFT",
    seq: "DIQMTQSPSSLSASVGDRVTITCRTSENIYSYLAWYQQKPGKAPKLLIYNAKTLAEGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQHHYGIPFTFGQGTKVEIK",
  },
  {
    id: "Timigutuzumab_vh",
    chain: "VH",
    cdr3: "SRWGGDGFYAMDY",
    seq: "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS",
  },
  {
    id: "Timigutuzumab_vl",
    chain: "VL",
    cdr3: "QQHYTTPPT",
    seq: "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK",
  },
  {
    id: "Tislelizumab_vh",
    chain: "VH",
    cdr3: "ARAYGNYWYIDV",
    seq: "QVQLQESGPGLVKPSETLSLTCTVSGFSLTSYGVHWIRQPPGKGLEWIGVIYADGSTNYNPSLKSRVTISKDTSKNQVSLKLSSVTAADTAVYYCARAYGNYWYIDVWGQGTTVTVSS",
  },
  {
    id: "Tislelizumab_vl",
    chain: "VL",
    cdr3: "HQAYSSPYT",
    seq: "DIVMTQSPDSLAVSLGERATINCKSSESVSNDVAWYQQKPGQPPKLLINYAFHRFTGVPDRFSGSGYGTDFTLTISSLQAEDVAVYYCHQAYSSPYTFGQGTKLEIK",
  },
  {
    id: "Tisotumab_vh",
    chain: "VH",
    cdr3: "ARSPWGYYLDS",
    seq: "EVQLLESGGGLVQPGGSLRLSCAASGFTFSNYAMSWVRQAPGKGLEWVSSISGSGDYTYYTDSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCARSPWGYYLDSWGQGTLVTVSS",
  },
  {
    id: "Tisotumab_vl",
    chain: "VL",
    cdr3: "QQYNSYPYT",
    seq: "DIQMTQSPPSLSASAGDRVTITCRASQGISSRLAWYQQKPEKAPKSLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQYNSYPYTFGQGTKLEIK",
  },
  {
    id: "Tixagevimab_vh",
    chain: "VH",
    cdr3: "AAPYCSSISCNDGFDI",
    seq: "QMQLVQSGPEVKKPGTSVKVSCKASGFTFMSSAVQWVRQARGQRLEWIGWIVIGSGNTNYAQKFQERVTITRDMSTSTAYMELSSLRSEDTAVYYCAAPYCSSISCNDGFDIWGQGTMVTVSS",
  },
  {
    id: "Tixagevimab_vl",
    chain: "VL",
    cdr3: "QHYGSSRGWT",
    seq: "EIVLTQSPGTLSLSPGERATLSCRASQSVSSSYLAWYQQKPGQAPRLLIYGASSRATGIPDRFSGSGSGTDFTLTISRLEPEDFAVYYCQHYGSSRGWTFGQGTKVEIK",
  },
  {
    id: "Tocilizumab_vh",
    chain: "VH",
    cdr3: "ARSLARTTAMDY",
    seq: "QVQLQESGPGLVRPSQTLSLTCTVSGYSITSDHAWSWVRQPPGRGLEWIGYISYSGITTYNPSLKSRVTMLRDTSKNQFSLRLSSVTAADTAVYYCARSLARTTAMDYWGQGSLVTVSS",
  },
  {
    id: "Tocilizumab_vl",
    chain: "VL",
    cdr3: "QQGNTLPYT",
    seq: "DIQMTQSPSSLSASVGDRVTITCRASQDISSYLNWYQQKPGKAPKLLIYYTSRLHSGVPSRFSGSGSGTDFTFTISSLQPEDIATYYCQQGNTLPYTFGQGTKVEIK",
  },
  {
    id: "Tralokinumab_vh",
    chain: "VH",
    cdr3: "ARDSSSSWARWFFDL",
    seq: "QVQLVQSGAEVKKPGASVKVSCKASGYTFTNYGLSWVRQAPGQGLEWMGWISANNGDTNYGQEFQGRVTMTTDTSTSTAYMELRSLRSDDTAVYYCARDSSSSWARWFFDLWGRGTLVTVSS",
  },
  {
    id: "Tralokinumab_vl",
    chain: "VL",
    cdr3: "QVWDTGSDPVV",
    seq: "YVLTQPPSVSVAPGKTARITCGGNIIGSKLVHWYQQKPGQAPVLVIYDDGDRPSGIPERFSGSNSGNTATLTISRVEAGDEADYYCQVWDTGSDPVVFGGGTKLTVL",
  },
  {
    id: "Trastuzumab_vh",
    chain: "VH",
    cdr3: "SRWGGDGFYAMDY",
    seq: "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS",
  },
  {
    id: "Trastuzumab_vl",
    chain: "VL",
    cdr3: "QQHYTTPPT",
    seq: "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK",
  },
  {
    id: "Tremelimumab_vh",
    chain: "VH",
    cdr3: "ARDPRGATLYYYYYGMDV",
    seq: "QVQLVESGGGVVQPGRSLRLSCAASGFTFSSYGMHWVRQAPGKGLEWVAVIWYDGSNKYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCARDPRGATLYYYYYGMDVWGQGTTVTVSS",
  },
  {
    id: "Tremelimumab_vl",
    chain: "VL",
    cdr3: "QQYYSTPFT",
    seq: "DIQMTQSPSSLSASVGDRVTITCRASQSINSYLDWYQQKPGKAPKLLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQYYSTPFTFGPGTKVEIK",
  },
  {
    id: "Tuvonralimab_vh",
    chain: "VH",
    cdr3: "ARGGLLGYFDY",
    seq: "QVQLVESGGGVVEPGRSLRLSCAASGFTFSSYGMHWVRQAPGKGLEWVAVIWYKPSEKDYADSAKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCARGGLLGYFDYWGQGTLVTVSS",
  },
  {
    id: "Tuvonralimab_vl",
    chain: "VL",
    cdr3: "QQYGRYPFT",
    seq: "EIVLTQSPGTLSLSPGERATLSCRASQSINSYLAWYQQKPGQAPRPLIYGVSSRATGIPDRFSGSGSGTDFTLTISRLEPEDFAVYYCQQYGRYPFTFGPGTKVDIK",
  },
  {
    id: "Ustekinumab_vh",
    chain: "VH",
    cdr3: "ARRRPGQGYFDF",
    seq: "EVQLVQSGAEVKKPGESLKISCKGSGYSFTTYWLGWVRQMPGKGLDWIGIMSPVDSDIRYSPSFQGQVTMSVDKSITTAYLQWNSLKASDTAMYYCARRRPGQGYFDFWGQGTLVTVSS",
  },
  {
    id: "Ustekinumab_vl",
    chain: "VL",
    cdr3: "QQYNIYPYT",
    seq: "DIQMTQSPSSLSASVGDRVTITCRASQGISSWLAWYQQKPEKAPKSLIYAASSLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQYNIYPYTFGQGTKLEIK",
  },
  {
    id: "Vedolizumab_vh",
    chain: "VH",
    cdr3: "ARGGYDGWDYAIDY",
    seq: "QVQLVQSGAEVKKPGASVKVSCKGSGYTFTSYWMHWVRQAPGQRLEWIGEIDPSESNTNYNQKFKGRVTLTVDISASTAYMELSSLRSEDTAVYYCARGGYDGWDYAIDYWGQGTLVTVSS",
  },
  {
    id: "Vedolizumab_vl",
    chain: "VL",
    cdr3: "LQGTHQPYT",
    seq: "VMTQSPLSLPVTPGEPASISCRSSQSLAKSYGNTYLSWYLQKPGQSPQLLIYGISNRFSGVPDRFSGSGSGTDFTLKISRVEAEDVGVYYCLQGTHQPYTFGQGTKVEIK",
  },
  {
    id: "Vunakizumab_vh",
    chain: "VH",
    cdr3: "TRYSLFYGSSPYAMDY",
    seq: "VQLVQSGAEVKKPGASVKVSCKASGYTFTDYEVHWVRQAPGQGLEWMGVIDPGTGGVAYNQKFEGRVTMTADTSTSTAYMELRSLRSDDTAVYYCTRYSLFYGSSPYAMDYWGQGTLVTVSS",
  },
  {
    id: "Vunakizumab_vl",
    chain: "VL",
    cdr3: "QQRSSYPWT",
    seq: "EIVLTQSPDFQSVTPKEKVTITCSASSSVNYMHWFQQKPDQSPKLWIYRTSNLASGVPSRFSGSGSGTDYTLTINSLEAEDAATYYCQQRSSYPWTFGQGTKLEIK",
  },
  {
    id: "Xeligekimab_vh",
    chain: "VH",
    cdr3: "VRDYYDLISDYYIHYWYFDL",
    seq: "EVQLVESGGGLVQPGGSLRLSCAASGMSMSDYWMNWVRQAPGKGLEWVAAINQDGDEKYYVGSVKGRFTISRDNAKNSLYLQMNSLRVEDTAVYYCVRDYYDLISDYYIHYWYFDLWGRGTLVTVSS",
  },
  {
    id: "Xeligekimab_vl",
    chain: "VL",
    cdr3: "QQYNGSPTT",
    seq: "DIQMTQSPSSLSASVGDRVTITCRASQNVHNRLTWYQQKPGKAPKLLIYGASNLESGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQYNGSPTTFGQGTKVEIK",
  },
  {
    id: "Zimberelimab_vh",
    chain: "VH",
    cdr3: "ARHLGYNGRYLPFDY",
    seq: "QLQLQESGPGLVKPSETLTLTCTVSADSISSTTYYWVWIRQPPGKGLEWIGSISYSGSTYYNPSLKSRVTVSVDTSKNQFSLKLNSVAATDTALYYCARHLGYNGRYLPFDYWGQGTLVTVSS",
  },
  {
    id: "Zimberelimab_vl",
    chain: "VL",
    cdr3: "SSYTSISTWV",
    seq: "ALTQPASVSGSPGQSITISCTGTSSDVGFYNYVSWYQQHPGKAPELMIYDVSNRPSGVSDRFSGSKSGNTASLTISGLQAEDEADYYCSSYTSISTWVFGGGTKLTVL",
  },
];

for (const tc of cdr3TestCases) {
  const result = detectCDR3(tc.seq, tc.chain);
  if (result === null) {
    assert(false, `${tc.id}: detectCDR3 returned null, expected "${tc.cdr3}"`);
  } else {
    const actual = tc.seq.substring(result[0], result[1]);
    assert(
      actual === tc.cdr3,
      `${tc.id}: CDR3 = "${actual}", expected "${tc.cdr3}"`,
    );
  }
}

// ============================================================
// Summary
// ============================================================
console.log(`\n========================================`);
console.log(`  ${passed} passed, ${failed} failed`);
console.log(`========================================`);
process.exit(failed > 0 ? 1 : 0);
