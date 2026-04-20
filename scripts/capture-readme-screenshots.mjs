/**
 * Captures README screenshots at 1200px width and syncs README.md <img height="..."> to PNG dimensions.
 *
 * One-time local setup: pnpm exec playwright install chromium
 */
import { readFileSync, writeFileSync } from "node:fs";
import { dirname, join } from "node:path";
import { fileURLToPath, pathToFileURL } from "node:url";
import { chromium } from "playwright";

const __dirname = dirname(fileURLToPath(import.meta.url));
const REPO_ROOT = join(__dirname, "..");
const INDEX_URL = pathToFileURL(join(REPO_ROOT, "index.html")).href;
const README_PATH = join(REPO_ROOT, "README.md");
const VIEWPORT_WIDTH = 1200;

/** Fixed viewport height (px): README framing; avoids huge full-document shots. Adjust per shot if needed. */
const SHOTS = [
  {
    filename: "screenshot_single_mab_analysis.png",
    viewportHeight: 1022,
    setup: async (page) => {
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
    },
  },
  {
    filename: "screenshot_multiple_mab_analysis.png",
    viewportHeight: 1080,
    setup: async (page) => {
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
    },
  },
];

function pngPixelSize(filePath) {
  const buf = readFileSync(filePath);
  // PNG: 8-byte signature, then 4-byte length, then 4-byte "IHDR", then 13-byte IHDR data (width, height, …).
  if (buf.length < 24) throw new Error("File too small for PNG: " + filePath);
  if (
    buf[0] !== 0x89 ||
    buf[1] !== 0x50 ||
    buf[2] !== 0x4e ||
    buf[3] !== 0x47
  ) {
    throw new Error("Not a PNG: " + filePath);
  }
  if (
    buf[12] !== 0x49 ||
    buf[13] !== 0x48 ||
    buf[14] !== 0x44 ||
    buf[15] !== 0x52
  ) {
    throw new Error("Missing IHDR chunk: " + filePath);
  }
  return {
    width: buf.readUInt32BE(16),
    height: buf.readUInt32BE(20),
  };
}

function patchReadmeHeights(heightsByBasename) {
  let md = readFileSync(README_PATH, "utf8");
  for (const [basename, height] of Object.entries(heightsByBasename)) {
    const escaped = basename.replace(/[.*+?^${}()|[\]\\]/g, "\\$&");
    const re = new RegExp(
      '(<img\\s+src="' +
        escaped +
        '"\\s+alt="[^"]*"\\s+width="' +
        VIEWPORT_WIDTH +
        '"\\s+height=")(\\d+)(")',
      "g",
    );
    md = md.replace(re, "$1" + height + "$3");
  }
  writeFileSync(README_PATH, md);
}

const browser = await chromium.launch({ headless: true });
try {
  const heightsByBasename = {};

  for (const shot of SHOTS) {
    const context = await browser.newContext({
      viewport: { width: VIEWPORT_WIDTH, height: 720 },
      deviceScaleFactor: 1,
    });
    const page = await context.newPage();
    await page.goto(INDEX_URL, { waitUntil: "domcontentloaded" });
    await shot.setup(page);
    await page.evaluate(
      () =>
        new Promise((resolve) => {
          requestAnimationFrame(() => requestAnimationFrame(resolve));
        }),
    );

    await page.evaluate(() => window.scrollTo(0, 0));
    const viewportH = shot.viewportHeight;
    if (!viewportH || viewportH < 200)
      throw new Error("Each shot needs viewportHeight (px)");
    await page.setViewportSize({ width: VIEWPORT_WIDTH, height: viewportH });
    await page.evaluate(
      () =>
        new Promise((resolve) => {
          requestAnimationFrame(() => requestAnimationFrame(resolve));
        }),
    );

    const outPath = join(REPO_ROOT, shot.filename);
    await page.screenshot({ path: outPath, fullPage: false });
    await context.close();

    const { width, height } = pngPixelSize(outPath);
    if (width !== VIEWPORT_WIDTH) {
      throw new Error(
        `Expected ${shot.filename} width ${VIEWPORT_WIDTH}, got ${width}`,
      );
    }
    heightsByBasename[shot.filename] = height;
  }

  patchReadmeHeights(heightsByBasename);
  console.log(
    "Wrote screenshots and updated README heights:",
    heightsByBasename,
  );
} finally {
  await browser.close();
}
