/**
 * Headless browser smoke test for LAMBS UI (Single + Multiple mAb Analysis tabs).
 *
 * Prerequisites:
 *   pnpm install
 *   pnpm exec playwright install chromium
 *
 * Usage:
 *   node scripts/verify-ui.mjs
 *   LAMBS_UI_OUT=/path/to/dir node scripts/verify-ui.mjs
 */
import { mkdirSync } from "node:fs";
import { join, dirname } from "node:path";
import { fileURLToPath, pathToFileURL } from "node:url";
import { chromium } from "playwright";

const __dirname = dirname(fileURLToPath(import.meta.url));
const REPO_ROOT = join(__dirname, "..");
const INDEX_URL = pathToFileURL(join(REPO_ROOT, "index.html")).href;
const OUT_DIR = process.env.LAMBS_UI_OUT || "/opt/cursor/artifacts/screenshots";

mkdirSync(OUT_DIR, { recursive: true });

function assert(condition, message) {
  if (!condition) throw new Error(message);
}

const browser = await chromium.launch({ headless: true });

try {
  const page = await browser.newPage({ viewport: { width: 1200, height: 1220 } });
  await page.goto(INDEX_URL, { waitUntil: "domcontentloaded" });

  // --- Single mAb Analysis ---
  console.log("=== Single mAb Analysis tab ===");
  await page.evaluate(() => {
    switchTab("single");
    loadExampleSingle();
    analyze();
  });
  await page
    .locator("#results")
    .getByText(/germline|VH|VL|liability/i)
    .first()
    .waitFor({ timeout: 15000 });

  const singleApi = await page.evaluate(() => {
    const r = LAMBS.analyzeMabReportFromRaw(
      document.getElementById("vh-input").value,
      document.getElementById("vl-input").value
    );
    if (!r.ok) return { ok: false, errors: r.errors };
    return {
      ok: true,
      kind: r.report.kind,
      vhVGene: r.report.chains.vh.variable.germline.vGene?.gene,
      vlVGene: r.report.chains.vl.variable.germline.vGene?.gene,
    };
  });
  console.log(JSON.stringify(singleApi, null, 2));
  assert(singleApi.ok && singleApi.kind === "single-mab", "Single tab analysis failed");

  const singleShot = join(OUT_DIR, "lambs-single-analysis-demo.png");
  await page.screenshot({ path: singleShot, fullPage: false });
  console.log("Screenshot:", singleShot);

  // --- Multiple mAb Analysis ---
  console.log("\n=== Multiple mAb Analysis tab ===");
  await page.evaluate(() => {
    switchTab("multi");
    loadExampleCSV();
    runMultiAnalysis();
  });
  await page
    .locator("#multi-results")
    .getByText(/cluster|alignment|consensus/i)
    .first()
    .waitFor({ timeout: 15000 });

  const multiApi = await page.evaluate(() => ({
    basketCount: typeof basket !== "undefined" ? basket.length : 0,
    clusterReport: LAMBS.lastClusterReport
      ? {
          kind: LAMBS.lastClusterReport.kind,
          inputCount: LAMBS.lastClusterReport.inputCount,
          clusterCount: LAMBS.lastClusterReport.clusterCount,
          firstClusterSize: LAMBS.lastClusterReport.clusters?.[0]?.size,
        }
      : null,
  }));
  console.log(JSON.stringify(multiApi, null, 2));
  assert(multiApi.basketCount === 4, "Expected 4 example sequences in basket");
  assert(multiApi.clusterReport?.kind === "cluster", "Multi tab cluster report missing");
  assert(multiApi.clusterReport.clusterCount >= 1, "Expected at least one cluster");

  const multiShot = join(OUT_DIR, "lambs-multiple-analysis-demo.png");
  await page.screenshot({ path: multiShot, fullPage: false });
  console.log("Screenshot:", multiShot);

  console.log("\nUI verification passed (Single + Multiple mAb Analysis).");
} finally {
  await browser.close();
}
