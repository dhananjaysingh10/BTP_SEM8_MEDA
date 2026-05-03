import fs from 'fs';

// ===================== SEEDED PRNG =====================
function mulberry32(seed) {
  return function () {
    seed |= 0; seed = seed + 0x6D2B79F5 | 0;
    let t = Math.imul(seed ^ seed >>> 15, 1 | seed);
    t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t;
    return ((t ^ t >>> 14) >>> 0) / 4294967296;
  };
}

// ===================== IMPORT ALGORITHMS =====================
const bmCode = fs.readFileSync('benchmark.js', 'utf-8');
const algoCode = bmCode.split('// ========== BENCHMARK SCRIPT ==========')[0];
const cleanCode = algoCode.replace(/^import.*$/gm, '');
const moduleExports = {};
eval(cleanCode + `\nmoduleExports.ratioApprox=ratioApprox;moduleExports.buildTree=buildTree;moduleExports.getStats=getStats;\n`);
const { ratioApprox, buildTree, getStats } = moduleExports;

// ===================== CONSTANTS =====================
const ALGOS = ["rma", "bs", "apdp", "larp", "ilp", "harp"];
const METRICS = ["m", "splits", "l", "maxL", "leaves", "d", "p"];
const LABELS = { m: "Mixes", splits: "Splits", l: "Dilution(l)", maxL: "MaxDilution", leaves: "Leaves", d: "Depth", p: "Parallelism" };
const LOWER_BETTER = new Set(["m", "splits", "leaves", "d", "p"]);
const HIGHER_BETTER = new Set(["l", "maxL"]);

// ===================== TEST CASE GENERATION =====================
function generateTestSuite(seed = 42) {
  const rng = mulberry32(seed);
  const ri = (lo, hi) => Math.floor(rng() * (hi - lo + 1)) + lo;
  const tests = [];
  let id = 0;

  const add = (cat, depth, ratios, note = "") => {
    tests.push({ id: ++id, cat, depth, ratios, note });
  };

  // ============ HAND-CRAFTED EDGE CASES ============
  // From literature: common test cases used in RMA/BS papers

  // Paper examples (Roy et al., Thies et al.)
  add("paper", 7, [2, 3, 5, 7, 11, 13, 87], "Roy2010-Fig3");
  add("paper", 7, [1, 1, 1, 1, 1, 1, 1, 1], "8-equal-d7");
  add("paper", 8, [1, 2, 4, 8, 16, 32, 64, 128], "powers-of-2");
  add("paper", 7, [1, 127], "extreme-2fluid");
  add("paper", 7, [64, 64], "perfect-half");
  add("paper", 8, [1, 1, 1, 1, 252], "1-dominant-d8");
  add("paper", 7, [10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 8], "13-near-equal");
  add("paper", 6, [1, 3, 5, 7, 9, 11], "odd-progression-d6");
  add("paper", 8, [3, 5, 7, 11, 13, 17, 19, 23], "primes-d8");

  // Exact half-sum possible (should be zero-split)
  add("edge", 7, [32, 32, 32, 32], "exact-quarter");
  add("edge", 7, [64, 32, 16, 8, 4, 2, 1, 1], "binary-weights");
  add("edge", 7, [40, 24, 40, 24], "symmetric-pair");
  add("edge", 8, [128, 64, 32, 16, 8, 4, 2, 2], "binary-d8");

  // No exact half-sum (forces splits)
  add("edge", 7, [33, 31, 33, 31], "no-exact-half");
  add("edge", 7, [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3], "32-equal-tiny");
  add("edge", 7, [1, 1, 126], "2-trace-1-bulk");
  add("edge", 7, [63, 63, 1, 1], "2-large-2-trace");

  // Extreme cases
  add("edge", 10, [1, 1023], "extreme-2fluid-d10");
  add("edge", 5, [8, 8, 8, 8], "small-depth-d5");
  add("edge", 6, [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1], "33-equal-d6");
  add("edge", 9, [1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 137], "fibonacci-d9");

  // Stress: many fluids high depth
  add("stress", 10, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15], "15-fluids-d10");
  add("stress", 8, Array.from({ length: 15 }, (_, i) => i + 1), "15-sequential-d8");

  // ============ RANDOM: FEW FLUIDS (2-3) ============
  for (let i = 0; i < 40; i++) {
    const n = ri(2, 3), d = ri(6, 8);
    add("rand-few", d, Array.from({ length: n }, () => ri(1, 100)));
  }

  // ============ RANDOM: MEDIUM FLUIDS (4-7) ============
  for (let i = 0; i < 80; i++) {
    const n = ri(4, 7), d = ri(6, 9);
    add("rand-med", d, Array.from({ length: n }, () => ri(1, 100)));
  }

  // ============ RANDOM: MANY FLUIDS (8-15) ============
  for (let i = 0; i < 60; i++) {
    const n = ri(8, 15), d = ri(7, 10);
    add("rand-many", d, Array.from({ length: n }, () => ri(1, 60)));
  }

  // ============ RANDOM: SKEWED (one dominant) ============
  for (let i = 0; i < 40; i++) {
    const n = ri(3, 8), d = ri(6, 9);
    const ratios = Array.from({ length: n }, () => ri(1, 15));
    ratios[0] = ri(60, 250);
    add("rand-skew", d, ratios);
  }

  // ============ RANDOM: NEAR-EQUAL ============
  for (let i = 0; i < 40; i++) {
    const n = ri(3, 12), d = ri(6, 9);
    const base = ri(5, 40);
    add("rand-equal", d, Array.from({ length: n }, () => base + ri(-3, 3)));
  }

  // ============ RANDOM: VERY LARGE DEPTH ============
  for (let i = 0; i < 20; i++) {
    const n = ri(4, 10), d = ri(9, 12);
    add("rand-deep", d, Array.from({ length: n }, () => ri(1, 100)));
  }

  // ============ RANDOM: TINY DEPTH ============
  for (let i = 0; i < 20; i++) {
    const n = ri(2, 5), d = ri(4, 6);
    add("rand-shallow", d, Array.from({ length: n }, () => ri(1, 50)));
  }

  return tests;
}

// ===================== RUNNER =====================
function runSuite(tests) {
  const results = [];
  for (const t of tests) {
    const adj = ratioApprox(t.ratios, t.depth);
    const fl = adj.map((_, i) => `x${i + 1}`);
    const row = { id: t.id, cat: t.cat, depth: t.depth, nFluids: t.ratios.length, note: t.note || "" };
    for (const algo of ALGOS) {
      row[algo] = getStats(buildTree(algo, adj, fl, t.depth));
    }
    results.push(row);
  }
  return results;
}

// ===================== STATISTICS =====================
function computeStats(results, filter) {
  const data = filter ? results.filter(filter) : results;
  const n = data.length;
  if (n === 0) return null;

  // Per-algo, per-metric: avg, median, std, min, max
  const stats = {};
  ALGOS.forEach(a => {
    stats[a] = {};
    METRICS.forEach(m => {
      const vals = data.map(r => r[a][m] || 0).sort((x, y) => x - y);
      const sum = vals.reduce((s, v) => s + v, 0);
      const avg = sum / n;
      const med = n % 2 === 1 ? vals[Math.floor(n / 2)] : (vals[n / 2 - 1] + vals[n / 2]) / 2;
      const variance = vals.reduce((s, v) => s + (v - avg) ** 2, 0) / n;
      stats[a][m] = { avg: +avg.toFixed(2), med, std: +Math.sqrt(variance).toFixed(2), min: vals[0], max: vals[n - 1] };
    });
  });

  // Winners per metric (by avg)
  const winners = {};
  METRICS.forEach(m => {
    const fn = LOWER_BETTER.has(m) ? Math.min : Math.max;
    const best = fn(...ALGOS.map(a => stats[a][m].avg));
    winners[m] = ALGOS.filter(a => stats[a][m].avg === best);
  });

  // Win counts: per test, who has the best value
  const wins = {};
  ALGOS.forEach(a => { wins[a] = {}; METRICS.forEach(m => wins[a][m] = 0); });
  for (const r of data) {
    for (const m of METRICS) {
      const vals = ALGOS.map(a => r[a][m] || 0);
      const best = LOWER_BETTER.has(m) ? Math.min(...vals) : Math.max(...vals);
      ALGOS.forEach((a, i) => { if (vals[i] === best) wins[a][m]++; });
    }
  }

  // Composite Pareto rank: for each test, rank algos by a composite score
  // Composite = normalized(mixes) + normalized(splits) + normalized(dilution) + normalized(maxL)
  const paretoWins = {};
  ALGOS.forEach(a => paretoWins[a] = 0);
  for (const r of data) {
    const composites = {};
    // Normalize each metric across algos for this test case
    for (const a of ALGOS) {
      let score = 0;
      for (const m of METRICS) {
        const vals = ALGOS.map(a2 => r[a2][m] || 0);
        const min = Math.min(...vals), max = Math.max(...vals);
        const v = r[a][m] || 0;
        if (max === min) score += 1;
        else if (LOWER_BETTER.has(m)) score += (max - v) / (max - min);
        else score += (v - min) / (max - min);
      }
      composites[a] = score;
    }
    const best = Math.max(...Object.values(composites));
    ALGOS.forEach(a => { if (composites[a] === best) paretoWins[a]++; });
  }

  return { stats, winners, wins, paretoWins, n };
}

// ===================== DISPLAY =====================
function pad(s, w) { return String(s).padStart(w); }
const SEP = "─";

function printAvgTable(title, st) {
  if (!st) return;
  console.log(`\n┌${"─".repeat(70)}`);
  console.log(`│ ${title} (n=${st.n})`);
  console.log(`├${"─".repeat(70)}`);
  const hdr = "  Metric".padEnd(16) + ALGOS.map(a => pad(a.toUpperCase(), 9)).join("") + "   Winner";
  console.log(hdr);
  console.log("  " + "─".repeat(hdr.length - 2));
  for (const m of METRICS) {
    let row = ("  " + LABELS[m]).padEnd(16);
    for (const a of ALGOS) row += pad(st.stats[a][m].avg, 9);
    row += "   " + st.winners[m].map(a => a.toUpperCase()).join(",");
    console.log(row);
  }
}

function printMedianTable(st) {
  if (!st) return;
  console.log(`\n  MEDIANS`);
  const hdr = "  Metric".padEnd(16) + ALGOS.map(a => pad(a.toUpperCase(), 9)).join("");
  console.log(hdr);
  console.log("  " + "─".repeat(hdr.length - 2));
  for (const m of METRICS) {
    let row = ("  " + LABELS[m]).padEnd(16);
    for (const a of ALGOS) row += pad(st.stats[a][m].med, 9);
    console.log(row);
  }
}

function printWinCounts(st) {
  if (!st) return;
  console.log(`\n  WIN COUNTS (per test case, ties counted for all)`);
  const hdr = "  Metric".padEnd(16) + ALGOS.map(a => pad(a.toUpperCase(), 9)).join("");
  console.log(hdr);
  console.log("  " + "─".repeat(hdr.length - 2));
  for (const m of METRICS) {
    let row = ("  " + LABELS[m]).padEnd(16);
    for (const a of ALGOS) {
      const pct = Math.round(st.wins[a][m] / st.n * 100);
      row += pad(`${st.wins[a][m]}(${pct}%)`, 9);
    }
    console.log(row);
  }
}

function printParetoWins(st) {
  if (!st) return;
  console.log(`\n  COMPOSITE SCORE WINS (best all-around per test)`);
  let row = "  ";
  for (const a of ALGOS) {
    const pct = Math.round(st.paretoWins[a] / st.n * 100);
    row += `${a.toUpperCase()}: ${st.paretoWins[a]}/${st.n} (${pct}%)  `;
  }
  console.log(row);
}

function printDelta(current, baseline) {
  console.log(`\n┌${"─".repeat(70)}`);
  console.log(`│ DELTA vs BASELINE`);
  console.log(`├${"─".repeat(70)}`);
  const hdr = "  Metric".padEnd(16) + ALGOS.map(a => pad(a.toUpperCase(), 11)).join("");
  console.log(hdr);
  console.log("  " + "─".repeat(hdr.length - 2));
  let anyChange = false;
  for (const m of METRICS) {
    let row = ("  " + LABELS[m]).padEnd(16);
    for (const a of ALGOS) {
      if (!baseline[a] || !baseline[a][m]) {
        row += pad("N/A", 11);
        continue;
      }
      const cur = current[a][m].avg, base = baseline[a][m].avg;
      const diff = +(cur - base).toFixed(2);
      if (diff !== 0) anyChange = true;
      const better = LOWER_BETTER.has(m) ? diff < 0 : diff > 0;
      const worse = LOWER_BETTER.has(m) ? diff > 0 : diff < 0;
      const sign = diff > 0 ? "+" : "";
      const mark = diff === 0 ? " =" : better ? " ✓" : " ✗";
      row += pad(`${sign}${diff}${mark}`, 11);
    }
    console.log(row);
  }
  if (!anyChange) console.log("\n  No changes detected — algorithms are identical to baseline.");
}

function printScoreboard(overall) {
  if (!overall) return;
  // Compute a composite rank: award points per metric
  // 1st = 6pts, 2nd = 4pts, 3rd = 3pts, 4th = 2pts, 5th = 1pt, 6th = 0pts
  const PTS = [6, 4, 3, 2, 1, 0];
  const points = {};
  ALGOS.forEach(a => points[a] = 0);
  for (const m of METRICS) {
    const ranked = [...ALGOS].sort((a, b) => {
      const va = overall.stats[a][m].avg, vb = overall.stats[b][m].avg;
      return LOWER_BETTER.has(m) ? va - vb : vb - va;
    });
    ranked.forEach((a, i) => points[a] += PTS[i]);
  }
  const sorted = [...ALGOS].sort((a, b) => points[b] - points[a]);

  console.log(`\n╔${"═".repeat(50)}`);
  console.log(`║  FINAL SCOREBOARD (${METRICS.length} metrics × [5,3,2,1,0] pts)`);
  console.log(`╠${"═".repeat(50)}`);
  sorted.forEach((a, i) => {
    const bar = "█".repeat(Math.round(points[a] / 2));
    const medal = i === 0 ? "🥇" : i === 1 ? "🥈" : i === 2 ? "🥉" : "  ";
    console.log(`║  ${medal} #${i + 1}  ${a.toUpperCase().padEnd(6)} ${String(points[a]).padStart(3)} pts  ${bar}`);
  });
  console.log(`╚${"═".repeat(50)}`);
}

// ===================== MAIN =====================
const BASELINE_FILE = "benchmark_baseline.json";
const RESULTS_FILE = "benchmark_latest.json";
const CSV_FILE = "benchmark_detailed.csv";

console.log("╔══════════════════════════════════════════════════╗");
console.log("║        MIXING TREE BENCHMARK FRAMEWORK          ║");
console.log("╚══════════════════════════════════════════════════╝");

const tests = generateTestSuite(42);
const cats = [...new Set(tests.map(t => t.cat))];
console.log(`\nTests: ${tests.length} total`);
for (const c of cats) console.log(`  ${c.padEnd(14)} ${tests.filter(t => t.cat === c).length} cases`);

console.log("\nRunning all algorithms...");
const t0 = Date.now();
const results = runSuite(tests);
const elapsed = Date.now() - t0;
console.log(`Done in ${(elapsed / 1000).toFixed(1)}s`);

// Overall stats
const overall = computeStats(results);
printAvgTable("OVERALL AVERAGES", overall);
printMedianTable(overall);
printWinCounts(overall);
printParetoWins(overall);

// Per category
for (const cat of cats) {
  const st = computeStats(results, r => r.cat === cat);
  printAvgTable(`CATEGORY: ${cat}`, st);
}

// Scoreboard
printScoreboard(overall);

// Save
const saveData = {
  timestamp: new Date().toISOString(), elapsed, n: overall.n,
  stats: overall.stats, wins: overall.wins, paretoWins: overall.paretoWins
};
fs.writeFileSync(RESULTS_FILE, JSON.stringify(saveData, null, 2));

let csv = "id,cat,depth,nFluids,note," + ALGOS.flatMap(a => METRICS.map(m => `${a}_${m}`)).join(",") + "\n";
for (const r of results) {
  csv += `${r.id},${r.cat},${r.depth},${r.nFluids},${r.note},`;
  csv += ALGOS.flatMap(a => METRICS.map(m => r[a][m] || 0)).join(",") + "\n";
}
fs.writeFileSync(CSV_FILE, csv);
console.log(`\nSaved: ${RESULTS_FILE}, ${CSV_FILE}`);

// Baseline comparison
if (fs.existsSync(BASELINE_FILE)) {
  const bl = JSON.parse(fs.readFileSync(BASELINE_FILE, 'utf-8'));
  console.log(`Baseline: ${bl.timestamp}`);
  printDelta(overall.stats, bl.stats);
} else {
  fs.writeFileSync(BASELINE_FILE, JSON.stringify(saveData, null, 2));
  console.log(`First run — saved as baseline (${BASELINE_FILE}).`);
  console.log(`After making changes, re-run to see deltas.`);
  console.log(`To reset baseline: copy benchmark_latest.json benchmark_baseline.json`);
}

console.log("\nDone.");
