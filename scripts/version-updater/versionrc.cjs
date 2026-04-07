"use strict";

module.exports = {
  packageFiles: [{ filename: "package.json", type: "json" }],
  bumpFiles: [
    { filename: "package.json", type: "json" },
    {
      filename: "index.html",
      updater: "./scripts/version-updater/index-html-version.cjs",
    },
  ],
};
