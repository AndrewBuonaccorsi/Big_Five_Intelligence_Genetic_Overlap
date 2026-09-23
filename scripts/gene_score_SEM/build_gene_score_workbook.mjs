#!/usr/bin/env node

import fs from "node:fs/promises";
import path from "node:path";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const [resultsPath, workbookPath, previewDirectory] = process.argv.slice(2);
if (!resultsPath || !workbookPath || !previewDirectory) {
  throw new Error(
    "Usage: build_gene_score_workbook.mjs " +
      "<results.json> <output.xlsx> <preview-directory>",
  );
}

const results = JSON.parse(await fs.readFile(resultsPath, "utf8"));
const workbook = Workbook.create();

const colors = {
  navy: "#17324D",
  teal: "#1F7A8C",
  white: "#FFFFFF",
  text: "#1F2937",
  line: "#D6DEE8",
  paleTeal: "#E8F3F4",
  paleRed: "#FDECEC",
  paleGold: "#FFF4D6",
  red: "#B42318",
};

const headerStyle = {
  fill: colors.navy,
  font: { bold: true, color: colors.white },
  verticalAlignment: "center",
  wrapText: true,
  borders: {
    bottom: { style: "medium", color: colors.teal },
  },
};

function columnName(columnIndex) {
  let value = columnIndex + 1;
  let name = "";
  while (value > 0) {
    const remainder = (value - 1) % 26;
    name = String.fromCharCode(65 + remainder) + name;
    value = Math.floor((value - 1) / 26);
  }
  return name;
}

function humanize(value) {
  return String(value)
    .replaceAll("_", " ")
    .replace(/\b\w/g, (character) => character.toUpperCase())
    .replace(/\bCfi\b/g, "CFI")
    .replace(/\bTli\b/g, "TLI")
    .replace(/\bRmsea\b/g, "RMSEA")
    .replace(/\bSrmr\b/g, "SRMR")
    .replace(/\bCi\b/g, "CI")
    .replace(/\bIq\b/g, "IQ")
    .replace(/\bSe\b/g, "SE");
}

function valueMatrix(records, keys) {
  return records.map((record) =>
    keys.map((key) => {
      const value = record[key];
      if (value === undefined || value === null) return null;
      if (typeof value === "string" && value.startsWith("=")) {
        return `'${value}`;
      }
      return value;
    }),
  );
}

function columnWidth(key, matrix) {
  const normalized = key.toLowerCase();
  if (matrix) return normalized === "indicator" ? 30 : 22;
  if (normalized === "item") return 38;
  if (normalized === "value") return 92;
  if (normalized === "model") return 38;
  if (normalized === "comparison") return 44;
  if (normalized === "analysis") return 32;
  if (normalized.includes("conclusion")) return 72;
  if (normalized.includes("loading_type")) return 40;
  if (normalized === "factor" || normalized === "indicator") return 32;
  if (normalized.includes("trait_")) return 30;
  if (normalized.includes("equality_constraint")) return 22;
  if (
    normalized.includes("gradient") ||
    normalized.includes("residual_variance")
  ) {
    return 22;
  }
  if (
    normalized.includes("parameters") ||
    normalized.includes("post_estimation")
  ) {
    return 21;
  }
  return 18;
}

function applyNumericFormats(sheet, keys, rowCount, matrix) {
  if (rowCount < 1) return;
  keys.forEach((key, index) => {
    if (matrix && index === 0) return;
    const normalized = key.toLowerCase();
    const range = sheet.getRangeByIndexes(1, index, rowCount, 1);
    if (
      normalized === "n_genes" ||
      normalized.includes("degrees_of_freedom") ||
      normalized.includes("parameters") ||
      normalized.includes("variances") ||
      normalized.includes("replicates") ||
      normalized === "jackknife_blocks" ||
      normalized === "block"
    ) {
      range.format.numberFormat = "#,##0";
    } else if (
      normalized.includes("p_value") ||
      normalized === "rmsea_p_close"
    ) {
      range.format.numberFormat = "0.000E+00";
    } else if (
      matrix ||
      normalized.includes("correlation") ||
      normalized.includes("standard_error") ||
      normalized.includes("estimate") ||
      normalized.includes("loading") ||
      normalized === "z" ||
      normalized.includes("ci_") ||
      normalized.includes("chi_square") ||
      normalized === "cfi" ||
      normalized === "tli" ||
      normalized === "rmsea" ||
      normalized === "srmr" ||
      normalized.includes("gradient")
    ) {
      range.format.numberFormat = "0.0000";
    }
  });
}

function writeTableSheet({
  name,
  records,
  tableName,
  freezeColumns = 0,
  matrix = false,
}) {
  if (!Array.isArray(records) || records.length === 0) {
    throw new Error(`${name} has no records.`);
  }

  const sheet = workbook.worksheets.add(name);
  const keys = Object.keys(records[0]);
  const finalColumn = columnName(keys.length - 1);
  const finalRow = records.length + 1;

  sheet.showGridLines = false;
  sheet.getRangeByIndexes(0, 0, 1, keys.length).values = [
    keys.map(humanize),
  ];
  sheet.getRangeByIndexes(0, 0, 1, keys.length).format = headerStyle;
  sheet.getRangeByIndexes(0, 0, 1, keys.length).format.rowHeight = 38;
  sheet.getRangeByIndexes(1, 0, records.length, keys.length).values =
    valueMatrix(records, keys);
  sheet.getRangeByIndexes(1, 0, records.length, keys.length).format = {
    font: { color: colors.text },
    verticalAlignment: "center",
    borders: {
      insideHorizontal: { style: "thin", color: colors.line },
    },
  };
  sheet.getRangeByIndexes(1, 0, records.length, keys.length).format.wrapText =
    true;
  applyNumericFormats(sheet, keys, records.length, matrix);

  keys.forEach((key, index) => {
    sheet.getRangeByIndexes(0, index, 1, 1).format.columnWidth =
      columnWidth(key, matrix);
  });

  const table = sheet.tables.add(
    `A1:${finalColumn}${finalRow}`,
    true,
    tableName,
  );
  table.showFilterButton = !matrix;

  if (keys.includes("admissible")) {
    const column = keys.indexOf("admissible");
    const range = sheet.getRangeByIndexes(1, column, records.length, 1);
    range.conditionalFormats.add("cellIs", {
      operator: "equal",
      formula: "FALSE",
      format: {
        fill: colors.paleRed,
        font: { bold: true, color: colors.red },
      },
    });
    range.conditionalFormats.add("cellIs", {
      operator: "equal",
      formula: "TRUE",
      format: {
        fill: colors.paleTeal,
        font: { bold: true, color: colors.teal },
      },
    });
  }

  if (keys.includes("p_value")) {
    const column = keys.indexOf("p_value");
    sheet
      .getRangeByIndexes(1, column, records.length, 1)
      .conditionalFormats.add("cellIs", {
        operator: "lessThan",
        formula: 0.05,
        format: { fill: colors.paleGold },
      });
  }

  sheet.freezePanes.freezeRows(1);
  if (freezeColumns > 0) {
    sheet.freezePanes.freezeColumns(freezeColumns);
  }

  return {
    sheet,
    keys,
    finalColumn,
    finalRow,
  };
}

const sheets = [
  writeTableSheet({
    name: "README",
    records: results.metadata,
    tableName: "GeneScoreMetadata",
  }),
  writeTableSheet({
    name: "Correlation Matrix",
    records: results.correlation_matrix,
    tableName: "GeneScoreCorrelationMatrix",
    freezeColumns: 1,
    matrix: true,
  }),
  writeTableSheet({
    name: "Jackknife SE Matrix",
    records: results.standard_error_matrix,
    tableName: "GeneScoreStandardErrorMatrix",
    freezeColumns: 1,
    matrix: true,
  }),
  writeTableSheet({
    name: "Pair Details",
    records: results.pair_details,
    tableName: "GeneScorePairDetails",
    freezeColumns: 2,
  }),
  writeTableSheet({
    name: "Model Fit",
    records: results.model_fit,
    tableName: "GeneScoreModelFit",
    freezeColumns: 1,
  }),
  writeTableSheet({
    name: "Model Comparison",
    records: results.model_comparison,
    tableName: "GeneScoreModelComparison",
  }),
  writeTableSheet({
    name: "Model Diagnostics",
    records: results.model_diagnostics,
    tableName: "GeneScoreModelDiagnostics",
    freezeColumns: 1,
  }),
  writeTableSheet({
    name: "Implied Trait Correlations",
    records: results.latent_factor_correlations,
    tableName: "GeneScoreImpliedCorrelations",
    freezeColumns: 4,
  }),
  writeTableSheet({
    name: "Hierarchical Loadings",
    records: results.hierarchical_loadings,
    tableName: "GeneScoreHierarchicalLoadings",
    freezeColumns: 3,
  }),
];

await fs.mkdir(path.dirname(workbookPath), { recursive: true });
await fs.mkdir(previewDirectory, { recursive: true });

for (const entry of sheets) {
  const previewLastRow = Math.min(entry.finalRow, 35);
  const preview = await workbook.render({
    sheetName: entry.sheet.name,
    range: `A1:${entry.finalColumn}${previewLastRow}`,
    scale: entry.sheet.name === "Model Diagnostics" ? 0.75 : 1,
    format: "png",
  });
  const safeName = entry.sheet.name.toLowerCase().replaceAll(" ", "_");
  await fs.writeFile(
    path.join(previewDirectory, `${safeName}.png`),
    new Uint8Array(await preview.arrayBuffer()),
  );

  const inspection = await workbook.inspect({
    kind: "table",
    range:
      `${entry.sheet.name}!A1:${entry.finalColumn}` +
      `${Math.min(entry.finalRow, 12)}`,
    include: "values,formulas",
    tableMaxRows: 12,
    tableMaxCols: entry.keys.length,
    maxChars: 8000,
  });
  console.log(inspection.ndjson);
}

const formulaErrors = await workbook.inspect({
  kind: "match",
  searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A",
  options: { useRegex: true, maxResults: 300 },
  summary: "final formula error scan",
});
console.log(formulaErrors.ndjson);

const output = await SpreadsheetFile.exportXlsx(workbook);
await output.save(workbookPath);
await fs.rm(`${workbookPath}.inspect.ndjson`, { force: true });
console.log(`Workbook written to ${workbookPath}`);
