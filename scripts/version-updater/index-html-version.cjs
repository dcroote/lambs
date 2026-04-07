"use strict";

const SPAN_RE = /<span id="lambs-app-version">([^<]*)<\/span>/;

module.exports.readVersion = function readVersion(contents) {
  const m = contents.match(SPAN_RE);
  if (!m) {
    throw new Error(
      'index.html: expected <span id="lambs-app-version">...</span> for release tooling',
    );
  }
  return m[1].trim();
};

module.exports.writeVersion = function writeVersion(contents, version) {
  if (!SPAN_RE.test(contents)) {
    throw new Error(
      "index.html: could not find lambs-app-version span to update",
    );
  }
  return contents.replace(
    SPAN_RE,
    `<span id="lambs-app-version">${version}</span>`,
  );
};
