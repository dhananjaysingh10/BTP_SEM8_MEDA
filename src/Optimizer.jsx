import { useState, useEffect } from "react";

const METRICS = [
  { key: "m", label: "Mixes (mix/split cycles)", lowerBetter: true, primary: true },
  { key: "splits", label: "Splits (fluid contamination)", lowerBetter: true, primary: true },
  { key: "l", label: "Dilution length (layout reuse)", lowerBetter: false, primary: true },
  { key: "maxL", label: "Max dilution chain", lowerBetter: false, primary: false },
  { key: "leaves", label: "Leaves (dispensing steps)", lowerBetter: true, primary: false },
  { key: "d", label: "Depth", lowerBetter: true, primary: false },
  { key: "p", label: "Parallelism needed", lowerBetter: true, primary: false },
];

const PRESETS = {
  cost: { m: 90, splits: 40, l: 20, maxL: 10, leaves: 80, d: 10, p: 10 },
  quality: { m: 30, splits: 80, l: 90, maxL: 70, leaves: 20, d: 10, p: 10 },
  balanced: { m: 60, splits: 60, l: 60, maxL: 40, leaves: 40, d: 20, p: 20 },
};

const ALGOS = ["rma", "bs", "apdp", "larp", "ilp"];
const ALGO_NAME = { rma: "RMA", bs: "BS", apdp: "AP-DP", larp: "LARP", ilp: "ILP" };

function normalize(statsMap, weights) {
  const norm = {};
  ALGOS.forEach(a => norm[a] = {});

  METRICS.forEach(mc => {
    const vals = ALGOS.map(a => statsMap[a]?.[mc.key] ?? 0);
    const min = Math.min(...vals), max = Math.max(...vals);
    ALGOS.forEach(a => {
      const v = statsMap[a]?.[mc.key] ?? 0;
      if (max === min) norm[a][mc.key] = 1;
      else if (mc.lowerBetter) norm[a][mc.key] = (max - v) / (max - min);
      else norm[a][mc.key] = (v - min) / (max - min);
    });
  });

  const wSum = Object.values(weights).reduce((a, b) => a + b, 0) || 1;
  const scores = {};
  ALGOS.forEach(a => {
    let s = 0;
    METRICS.forEach(mc => s += weights[mc.key] * norm[a][mc.key]);
    scores[a] = s / wSum;
  });

  const ranked = [...ALGOS].sort((a, b) => scores[b] - scores[a]);
  return { norm, scores, ranked, best: ranked[0] };
}

export default function OptimizerPage({ raw, setRaw, depth, setDepth, buildTreeFn, getStatsFn, ratioApproxFn, layoutTreeFn, TreeViewCmp }) {
  const [weights, setWeights] = useState({ m: 70, splits: 50, l: 80, maxL: 30, leaves: 30, d: 10, p: 10 });
  const [showAdvanced, setShowAdvanced] = useState(false);
  const [results, setResults] = useState(null);
  const [allF, setAllF] = useState([]);
  const [winnerTree, setWinnerTree] = useState(null);

  const setW = (key, val) => setWeights(w => ({ ...w, [key]: val }));

  const run = () => {
    const r = raw.split(",").map(Number).filter(n => !isNaN(n) && n > 0);
    if (r.length < 2) return;
    const a = ratioApproxFn(r, depth);
    const fl = a.map((_, i) => `x${i + 1}`);
    setAllF(fl);
    const stats = {}, trees = {};
    ALGOS.forEach(algo => {
      const root = buildTreeFn(algo, a, fl, depth);
      if (root) { stats[algo] = getStatsFn(root); trees[algo] = root; }
    });
    const res = normalize(stats, weights);
    res.stats = stats;
    res.trees = trees;
    setResults(res);
    setWinnerTree(layoutTreeFn(trees[res.best]));
  };

  useEffect(run, []);

  const primary = METRICS.filter(m => m.primary);
  const advanced = METRICS.filter(m => !m.primary);

  return (
    <div style={{ maxWidth: 900, margin: "0 auto", padding: "20px 24px" }}>
      <h2 style={{ margin: "0 0 4px", fontSize: 18 }}>Smart Optimizer</h2>
      <p style={{ margin: "0 0 16px", color: "#666", fontSize: 13 }}>
        Set importance weights for each metric. The tool compares all 5 algorithms and recommends the best one.
      </p>

      {/* Inputs */}
      <div style={{ display: "flex", gap: 12, alignItems: "end", marginBottom: 16, flexWrap: "wrap" }}>
        <div>
          <label style={{ fontSize: 12, display: "block", marginBottom: 2, color: "#666" }}>Fluid ratios (comma-separated)</label>
          <input value={raw} onChange={e => setRaw(e.target.value)}
            style={{ padding: "6px 10px", border: "1px solid #ccc", borderRadius: 3, fontSize: 13, width: 200 }} />
        </div>
        <div>
          <label style={{ fontSize: 12, display: "block", marginBottom: 2, color: "#666" }}>Depth (d)</label>
          <input type="number" value={depth} min={2} max={14} onChange={e => setDepth(+e.target.value)}
            style={{ padding: "6px 10px", border: "1px solid #ccc", borderRadius: 3, fontSize: 13, width: 60 }} />
        </div>
        <button onClick={run} style={{ padding: "7px 20px", background: "#2563eb", color: "#fff", border: "none", borderRadius: 3, fontSize: 13, fontWeight: 600 }}>
          Find Best Algorithm
        </button>
      </div>

      {/* Presets */}
      <div style={{ marginBottom: 12 }}>
        <span style={{ fontSize: 12, color: "#666", marginRight: 8 }}>Quick presets:</span>
        {Object.entries(PRESETS).map(([id, p]) => (
          <button key={id} onClick={() => setWeights({ ...p })}
            style={{ marginRight: 6, padding: "4px 12px", border: "1px solid #ccc", borderRadius: 3, background: "#f9f9f9", fontSize: 12 }}>
            {id.charAt(0).toUpperCase() + id.slice(1)}
          </button>
        ))}
      </div>

      {/* Weights */}
      <fieldset style={{ border: "1px solid #ddd", borderRadius: 4, padding: "12px 16px", marginBottom: 12 }}>
        <legend style={{ fontSize: 13, fontWeight: 600, padding: "0 6px" }}>Priority Weights</legend>
        {primary.map(m => (
          <div key={m.key} style={{ display: "flex", alignItems: "center", gap: 10, marginBottom: 8 }}>
            <span style={{ width: 240, fontSize: 13 }}>
              {m.label} <span style={{ color: "#999", fontSize: 11 }}>({m.lowerBetter ? "↓ lower better" : "↑ higher better"})</span>
            </span>
            <input type="range" min={0} max={100} value={weights[m.key]} onChange={e => setW(m.key, +e.target.value)}
              style={{ flex: 1 }} />
            <span style={{ width: 30, textAlign: "right", fontSize: 13, fontWeight: 600 }}>{weights[m.key]}</span>
          </div>
        ))}

        <button onClick={() => setShowAdvanced(s => !s)}
          style={{ background: "none", border: "none", color: "#2563eb", fontSize: 12, padding: 0, marginTop: 4 }}>
          {showAdvanced ? "▼ Hide" : "▶ Show"} advanced metrics ({advanced.length})
        </button>

        {showAdvanced && advanced.map(m => (
          <div key={m.key} style={{ display: "flex", alignItems: "center", gap: 10, marginTop: 8 }}>
            <span style={{ width: 240, fontSize: 13 }}>
              {m.label} <span style={{ color: "#999", fontSize: 11 }}>({m.lowerBetter ? "↓ lower" : "↑ higher"})</span>
            </span>
            <input type="range" min={0} max={100} value={weights[m.key]} onChange={e => setW(m.key, +e.target.value)}
              style={{ flex: 1 }} />
            <span style={{ width: 30, textAlign: "right", fontSize: 13, fontWeight: 600 }}>{weights[m.key]}</span>
          </div>
        ))}
      </fieldset>

      {/* Results */}
      {results && (
        <>
          {/* Winner */}
          <div style={{ padding: "12px 16px", border: "2px solid #2563eb", borderRadius: 4, marginBottom: 16, background: "#f0f5ff" }}>
            <strong>Recommended: {ALGO_NAME[results.best]}</strong>
            <span style={{ marginLeft: 12, color: "#555" }}>
              Score: {(results.scores[results.best] * 100).toFixed(1)}%
            </span>
          </div>

          {/* Comparison table */}
          <h3 style={{ fontSize: 14, margin: "0 0 8px" }}>Comparison</h3>
          <table style={{ marginBottom: 16 }}>
            <thead>
              <tr>
                <th style={{ textAlign: "left" }}>Algorithm</th>
                {METRICS.map(m => (
                  <th key={m.key} title={m.label}>
                    {m.key}<br />
                    <span style={{ fontSize: 10, fontWeight: 400, color: "#888" }}>w={weights[m.key]}</span>
                  </th>
                ))}
                <th>Score</th>
                <th>Rank</th>
              </tr>
            </thead>
            <tbody>
              {results.ranked.map((algo, i) => (
                <tr key={algo} style={{ background: i === 0 ? "#f0f5ff" : "transparent", fontWeight: i === 0 ? 600 : 400 }}>
                  <td style={{ textAlign: "left", fontWeight: 600 }}>
                    {i === 0 ? "★ " : ""}{ALGO_NAME[algo]}
                  </td>
                  {METRICS.map(m => {
                    const raw = results.stats[algo]?.[m.key] ?? 0;
                    const n = results.norm[algo][m.key];
                    return (
                      <td key={m.key}>
                        {raw}
                        <br />
                        <span style={{ fontSize: 10, color: n > 0.7 ? "#16a34a" : n < 0.3 ? "#dc2626" : "#999" }}>
                          {(n * 100).toFixed(0)}%
                        </span>
                      </td>
                    );
                  })}
                  <td style={{ fontWeight: 700 }}>{(results.scores[algo] * 100).toFixed(1)}%</td>
                  <td>#{i + 1}</td>
                </tr>
              ))}
            </tbody>
          </table>

          {/* Winner tree */}
          <h3 style={{ fontSize: 14, margin: "0 0 8px" }}>Best Tree — {ALGO_NAME[results.best]}</h3>
          <p style={{ margin: "0 0 6px", fontSize: 11, color: "#999" }}>Scroll to zoom, drag to pan</p>
          <div style={{ height: 500, border: "1px solid #ddd", borderRadius: 4, background: "#fff" }}>
            <TreeViewCmp treeData={winnerTree} allFluids={allF} />
          </div>
        </>
      )}
    </div>
  );
}
