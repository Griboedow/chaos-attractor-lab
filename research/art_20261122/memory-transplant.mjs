import { readFileSync, writeFileSync } from "node:fs";
import { execFileSync } from "node:child_process";
import { fileURLToPath } from "node:url";
import { donorMatches, makeAdvance } from "./memory_transplant_map.mjs";

const results = new URL("../../results/art_20261122/memory-transplant/", import.meta.url);
const archive = new URL("../../results/wow_search_20261122/orbit/case-86780.json", import.meta.url);
const { parameters } = JSON.parse(readFileSync(archive, "utf8"));
const donors = JSON.parse(readFileSync(new URL("donors.json", results), "utf8"));
const { fields } = donors;
const advance = makeAdvance(parameters);
const threshold = 7.5;
const windows = [[6145, 8192], [24577, 32768]];
const doses = Array.from({ length: 193 }, (_, i) => 5.5 + 3.7 * i / 192);

function classify(x, y, memory) {
  let state = [x, y, memory], previous = 0;
  const labels = [];
  for (const [first, last] of windows) {
    const run = advance(state, last - previous, first - previous);
    labels.push(run.high < threshold ? "L" : run.low > threshold ? "U" : "?");
    state = run.state;
    previous = last;
  }
  return labels[0] === labels[1] ? labels[0] : "?";
}

if (process.argv.includes("--check-donors")) {
  const matches = donorMatches(donors, advance);
  console.log(`${matches.filter(Boolean).length}/${fields.length} donors match their published coordinates exactly.`);
  process.exit(matches.every(Boolean) ? 0 : 1);
}

const replayFile = fileURLToPath(new URL("memory_transplant_replay.py", import.meta.url));
const replay = JSON.parse(execFileSync(process.env.PYTHON || "python", [replayFile],
  { encoding: "utf8", windowsHide: true }));
function switchCount(labels) {
  let count = 0;
  for (let i = 1; i < labels.length; i++) {
    if (labels[i - 1] !== "?" && labels[i] !== "?" && labels[i - 1] !== labels[i]) count++;
  }
  return count;
}
const rows = fields.map(({ id, state: [x, y] }, rowIndex) => {
  const labels = doses.map(value => classify(x, y, value));
  const other = replay.doses[rowIndex];
  const agreed = labels.filter((label, i) => label !== "?" && label === other[i]);
  return {
    field: id, doses: labels.join(""), numpy_doses: other, switches: switchCount(labels),
    switches_among_agreed: switchCount(agreed)
  };
});
const swaps = fields.map(({ id, state: [x, y] }, i) => ({
  field: id,
  outcomes: fields.map(({ state }) => classify(x, y, state[2])).join(""),
  numpy_outcomes: replay.swaps[i]
}));
const packet = {
  case: 86780, parameters, threshold, windows, numpy_version: replay.numpy,
  dose_grid: { first: doses[0], last: doses.at(-1), count: doses.length },
  donor_order: fields.map(field => field.id), swaps, rows
};
writeFileSync(new URL("outcomes.json", results), JSON.stringify(packet, null, 2) + "\n");

const ink = "#182e37", paper = "#fbf9f2", grey = "#526269";
const colors = { L: "#146e87", U: "#efaa45", "?": "#a93252" };
const svg = [
  '<svg xmlns="http://www.w3.org/2000/svg" width="1120" height="665" viewBox="0 0 1120 665" role="img" aria-labelledby="title desc">',
  '<title id="title">Transplanting memory into an optical cavity</title>',
  '<desc id="desc">Each row holds one complex field constant. The horizontal axis is the initial detuning. Teal and amber cells indicate low and high late detuning. Red marks a disagreement between the two calculations. The last two columns count switches in the reference scan and among agreeing cells.</desc>'
];
function rect(x, y, w, h, fill) {
  svg.push(`<rect x="${x}" y="${y}" width="${w}" height="${h}" fill="${fill}"/>`);
}
function label(x, y, content, size = 15, color = ink, align = "start", weight = 400) {
  svg.push(`<text x="${x}" y="${y}" font-family="Arial, Helvetica, sans-serif" font-size="${size}" font-weight="${weight}" fill="${color}" text-anchor="${align}">${content}</text>`);
}
function strip(values, first, last, x, y, width, height) {
  const cellWidth = width / (last - first + 1);
  for (let start = first; start <= last;) {
    let end = start + 1;
    while (end <= last && values[end] === values[start]) end++;
    // A small overlap closes viewer seams between adjacent polygons.
    rect(x + (start - first) * cellWidth, y, (end - start) * cellWidth + 0.05,
         height, colors[values[start]]);
    start = end;
  }
}
const whole = { x: 95, width: 243, first: 0, last: 192 };
const zoom = { x: 385, width: 553, first: 51, last: 89 };
const maxSwitches = Math.max(...rows.map(row => row.switches));
rect(0, 0, 1120, 665, paper);
label(35, 49, "Transplanting memory into an optical cavity", 36, ink, "start", 700);
label(36, 83, "One field per row. Only the starting memory changes from left to right.", 20, grey);
rect(36, 104, 18, 12, colors.L);  label(63, 116, "Lower regime", 16);
rect(197, 104, 18, 12, colors.U); label(224, 116, "Upper regime", 16);
rect(359, 104, 18, 12, colors["?"]); label(386, 116, "Binary64 runs disagree", 16);

label(34, 161, "Field", 16, ink, "start", 700);
label(whole.x, 161, "Full dose range", 18, ink, "start", 700);
label(zoom.x, 161, "Boundary in detail", 18, ink, "start", 700);
label(1027, 161, "Binary64 switches", 16, ink, "middle", 700);
label(988, 191, "All", 15, grey, "middle");
label(1060, 191, "Agreed", 15, grey, "middle");
function ypos(i) { return 207 + i * 21 + (i >= 8 ? 14 : 0); }
for (const panel of [whole, zoom]) {
  const n = panel.last - panel.first + 1;
  const ticks = panel === whole ? [0, 96, 192] : [51, 61, 70, 80, 89];
  for (const i of ticks) {
    const x = panel.x + (i - panel.first + 0.5) * panel.width / n;
    label(x, 191, doses[i].toFixed(2), 15, grey, "middle");
  }
}
for (const [i, row] of rows.entries()) {
  const y = ypos(i), bold = row.switches === maxSwitches;
  if (bold) rect(24, y - 2, 1067, 20, "#e0e5dc");
  label(34, y + 13, row.field, 16, ink, "start", bold ? 700 : 400);
  for (const p of [whole, zoom]) {
    strip(row.doses, p.first, p.last, p.x, y, p.width, 16);
    const width = p.width / (p.last - p.first + 1);
    for (let j = p.first; j <= p.last; j++) {
      if (row.doses[j] !== row.numpy_doses[j] || row.doses[j] === "?") {
        rect(p.x + (j - p.first) * width, y + 11, width, 5, colors["?"]);
      }
    }
  }
  label(988, y + 13, row.switches, 16, ink, "middle", bold ? 700 : 400);
  label(1060, y + 13, row.switches_among_agreed, 16, ink, "middle", bold ? 700 : 400);
}
label(36, 598, "L01–L08 come from the lower regime. U01–U08 come from the upper regime.", 16, grey);
label(36, 627, "All counts the switches in the reference run. Agreed counts the switches that remain after omitting red cells.", 16, grey);
svg.push("</svg>");
const image = fileURLToPath(new URL("memory-transplant.png", results));
execFileSync("magick", ["-background", "white", "-density", "132", "svg:-",
  "-strip", "-define", "png:exclude-chunk=date,time", "PNG24:" + image],
  { input: svg.join("\n"), windowsHide: true });
