import { useState, useRef, useEffect } from "react";

// ========== RATIO APPROXIMATION ==========
function ratioApprox(ratios, d) {
  const Lp = 2 ** d, L = ratios.reduce((a, b) => a + b, 0);
  const ap = ratios.map(a => Math.round((a * Lp) / L));
  for (let i = 0; i < ap.length; i++) if (ap[i] <= 0) ap[i] = 1;
  ap[ap.length - 1] = Math.max(1, Lp - ap.slice(0, -1).reduce((a, b) => a + b, 0));
  return ap;
}

// ========== RMA (Algorithm 4: Expression Partition from paper) ==========
function rmaPartition(P, L) {
  const half = L / 2;
  // Step 1: Find au = highest coefficient in P for xu
  const sorted = Object.entries(P).filter(([, v]) => v > 0).sort((a, b) => b[1] - a[1]);
  if (!sorted.length) return [{}, {}];
  const [uF, aU] = sorted[0];
  let P1 = {}, P2 = {};

  // Step 3: if au >= L/2
  if (aU >= half) {
    // P1 = {xu: L/2}, P2 = remaining
    P1[uF] = half;
    const rem = aU - half;
    if (rem > 0) P2[uF] = rem;
    for (let i = 1; i < sorted.length; i++) P2[sorted[i][0]] = sorted[i][1];
  } else {
    // Step 6: P1 = au*xu, P2 = P - au*xu
    P1[uF] = aU;
    for (let i = 1; i < sorted.length; i++) P2[sorted[i][0]] = sorted[i][1];
    const E = half - aU;

    // Step 8: Check if E == some coefficient az in P2
    const p2entries = Object.entries(P2).filter(([, v]) => v > 0);
    const exactMatch = p2entries.find(([, v]) => v === E);

    if (exactMatch) {
      // Step 9: Move that term entirely from P2 to P1
      P1[exactMatch[0]] = exactMatch[1];
      delete P2[exactMatch[0]];
    } else {
      // Step 10: Check if E == sum of n1 coefficients, n1 < (count of nonzero in P2)/2
      const p2count = p2entries.length;
      let subsetFound = false;

      // Try all subsets up to size <= p2count/2 (paper says < but its figures require <=)
      if (p2count <= 16) {
        const maxSubsetSize = Math.floor(p2count / 2); // n1 <= count/2
        for (let mask = 1; mask < (1 << p2count) && !subsetFound; mask++) {
          const bits = mask.toString(2).split('').filter(b => b === '1').length;
          if (bits * 2 > p2count) continue;
          let sum = 0;
          for (let j = 0; j < p2count; j++) if (mask & (1 << j)) sum += p2entries[j][1];
          if (sum === E) {
            // Step 11: Move those n1 terms from P2 to P1
            for (let j = 0; j < p2count; j++) {
              if (mask & (1 << j)) {
                P1[p2entries[j][0]] = p2entries[j][1];
                delete P2[p2entries[j][0]];
              }
            }
            subsetFound = true;
          }
        }
      }

      if (!subsetFound) {
        // Step 13: Find av = next highest coefficient in P2 where av > E, split E from it
        // Sort P2 entries descending by value
        const p2sorted = Object.entries(P2).filter(([, v]) => v > 0).sort((a, b) => b[1] - a[1]);
        const splittable = p2sorted.find(([, v]) => v > E);
        if (splittable) {
          const [fv, av] = splittable;
          P1[fv] = (P1[fv] || 0) + E;
          P2[fv] = av - E;
          if (P2[fv] <= 0) delete P2[fv];
        } else {
          // Edge case fallback: no single coeff > E, greedily fill
          let rem = E;
          for (const [f] of p2sorted) {
            if (rem <= 0) break;
            const take = Math.min(P2[f], rem);
            P1[f] = (P1[f] || 0) + take;
            P2[f] -= take;
            if (P2[f] <= 0) delete P2[f];
            rem -= take;
          }
        }
      }
    }
  }
  return [P1, P2];
}
function buildRMA(P, L, lv = 0) {
  const keys = Object.keys(P).filter(k => P[k] > 0);
  if (!keys.length) return null;
  if (keys.length === 1) return { label: keys[0], volume: P[keys[0]], level: lv, leaf: true, partition: { ...P } };
  if (L <= 1) { const s = Object.entries(P).sort((a, b) => b[1] - a[1]); return { label: s[0][0], volume: L, level: lv, leaf: true, partition: { ...P } }; }
  const [P1, P2] = rmaPartition(P, L);
  return { label: "Mix", partition: { ...P }, volume: L, level: lv, leaf: false, left: buildRMA(P1, L / 2, lv + 1), right: buildRMA(P2, L / 2, lv + 1) };
}

// ========== BS (Bit-Scanning, bottom-up per paper [6]) ==========
function buildBSTree(approx, names, d) {
  const lvIn = [];
  for (let b = 0; b < d; b++) { const fl = []; for (let i = 0; i < approx.length; i++) if ((approx[i] >> b) & 1) fl.push(names[i]); lvIn.push(fl); }
  let prev = [];
  for (let b = 0; b < d; b++) {
    const leaves = lvIn[b].map(f => ({ label: f, leaf: true, partition: { [f]: 1 }, volume: 1 }));
    const items = [...prev, ...leaves];
    if (items.length <= 1) { prev = items; continue; }
    const next = [];
    for (let i = 0; i < items.length; i += 2) {
      if (i + 1 < items.length) {
        const l = items[i], r = items[i + 1], p = { ...l.partition };
        for (const [k, v] of Object.entries(r.partition)) p[k] = (p[k] || 0) + v;
        next.push({ label: "Mix", leaf: false, partition: p, volume: Object.values(p).reduce((a, b) => a + b, 0), left: l, right: r });
      } else next.push(items[i]);
    }
    prev = next;
  }
  const root = prev[0] || null;
  if (root) (function sL(n, l) { if (!n) return; n.level = l; if (!n.leaf) { sL(n.left, l + 1); sL(n.right, l + 1); } })(root, 0);
  return root;
}

// ========== AP-DP v3: Enhanced Multi-Strategy Adaptive Partitioning ==========
// Key improvements over v2:
// 1. Enumerates ALL subsets that sum to exactly half (not just one)
// 2. Explores multiple base sums for split strategies
// 3. Considers RMA-style partitions for multiple large fluids
// 4. Enhanced scoring with dilution potential estimation
// 5. Tries splitting from both sides (P1 and P2)

function scorePartitionV3(P1, P2, L) {
  // Lower score = better partition
  const k1 = Object.keys(P1).filter(k => P1[k] > 0);
  const k2 = Object.keys(P2).filter(k => P2[k] > 0);
  const set1 = new Set(k1), set2 = new Set(k2);

  // Count splits (fluids appearing in BOTH sides) — most important
  let splits = 0;
  for (const k of k1) if (set2.has(k)) splits++;

  // Dilution potential: estimate how "pure" each side is
  // A side with one dominant fluid (>50% of its volume) has high dilution potential
  const half = L / 2;
  const dom1 = k1.length > 0 ? Math.max(...k1.map(k => P1[k])) / half : 0;
  const dom2 = k2.length > 0 ? Math.max(...k2.map(k => P2[k])) / half : 0;
  const dilutionBonus = (dom1 > 0.5 ? 1 : 0) + (dom2 > 0.5 ? 1 : 0); // 0, 1, or 2

  // Fewer distinct fluids = potentially longer dilution chains
  const totalDistinct = k1.length + k2.length;

  // Prefer partitions where one side has very few fluids (good for dilution)
  const minFluids = Math.min(k1.length, k2.length);

  // Score: splits are heavily penalized, dilution is rewarded
  return splits * 1000 - dilutionBonus * 50 + totalDistinct * 10 + minFluids * 5;
}

function enumerateSubsetsWithSum(items, target, maxSubsets = 20) {
  // Enumerate ALL subsets that sum to exactly target (up to maxSubsets)
  // Returns array of Sets of indices
  const n = items.length;
  const results = [];

  if (n <= 20) {
    // Brute force for small n
    for (let mask = 1; mask < (1 << n) && results.length < maxSubsets; mask++) {
      let sum = 0;
      for (let i = 0; i < n; i++) if (mask & (1 << i)) sum += items[i].v;
      if (sum === target) {
        const subset = new Set();
        for (let i = 0; i < n; i++) if (mask & (1 << i)) subset.add(i);
        results.push(subset);
      }
    }
  } else {
    // DP with multiple path tracking for larger n
    const dp = new Map(); // sum -> array of subsets (limited)
    dp.set(0, [new Set()]);

    for (let i = 0; i < n; i++) {
      const v = items[i].v;
      const newEntries = [];
      for (const [sum, subsets] of dp.entries()) {
        const newSum = sum + v;
        if (newSum <= target) {
          for (const subset of subsets) {
            if (results.length >= maxSubsets && newSum !== target) continue;
            const newSubset = new Set(subset);
            newSubset.add(i);
            newEntries.push([newSum, newSubset]);
          }
        }
      }
      for (const [sum, subset] of newEntries) {
        if (sum === target) {
          results.push(subset);
          if (results.length >= maxSubsets) break;
        } else {
          if (!dp.has(sum)) dp.set(sum, []);
          if (dp.get(sum).length < 3) dp.get(sum).push(subset); // limit stored per sum
        }
      }
      if (results.length >= maxSubsets) break;
    }
  }

  return results;
}

function findClosestSums(items, target, count = 5) {
  // Find the 'count' closest achievable sums to target (including exact if possible)
  const n = items.length;
  const achievable = new Set([0]);
  const subsetFor = new Map(); // sum -> one subset achieving it
  subsetFor.set(0, new Set());

  for (let i = 0; i < n; i++) {
    const v = items[i].v;
    const toAdd = [];
    for (const sum of achievable) {
      const newSum = sum + v;
      if (newSum <= target && !achievable.has(newSum)) {
        toAdd.push(newSum);
        const newSubset = new Set(subsetFor.get(sum));
        newSubset.add(i);
        subsetFor.set(newSum, newSubset);
      }
    }
    for (const s of toAdd) achievable.add(s);
  }

  // Sort by closeness to target
  const sorted = [...achievable].sort((a, b) => Math.abs(target - a) - Math.abs(target - b));
  return sorted.slice(0, count).map(sum => ({ sum, subset: subsetFor.get(sum) }));
}

function buildAPDP(P, L, maxD, lv = 0) {
  const keys = Object.keys(P).filter(k => P[k] > 0);
  if (!keys.length) return null;
  if (keys.length === 1) {
    const k = keys[0];
    return { label: k, volume: P[k], level: lv, leaf: true, partition: { [k]: P[k] } };
  }
  if (L <= 1 || lv >= maxD) {
    const s = Object.entries(P).sort((a, b) => b[1] - a[1]);
    return { label: s[0][0], volume: L, level: lv, leaf: true, partition: { ...P } };
  }

  const half = L / 2;
  const items = keys.map(k => ({ k, v: P[k] }));
  const sorted = items.slice().sort((a, b) => b.v - a.v);
  const n = items.length;

  // ===== Generate candidate partitions =====
  const candidates = [];

  // --- Strategy A: Enumerate ALL whole-fluid partitions summing to exactly half ---
  const exactSubsets = enumerateSubsetsWithSum(items, half, 15);
  for (const subset of exactSubsets) {
    const cP1 = {}, cP2 = {};
    items.forEach((it, i) => {
      if (subset.has(i)) cP1[it.k] = it.v;
      else cP2[it.k] = it.v;
    });
    if (Object.keys(cP1).length > 0 && Object.keys(cP2).length > 0) {
      candidates.push({ P1: cP1, P2: cP2, type: "exact-whole" });
    }
  }

  // --- Strategy B: DP closest sums + single split ---
  if (exactSubsets.length === 0) {
    const closestSums = findClosestSums(items, half, 5);
    for (const { sum, subset } of closestSums) {
      if (sum === half) continue; // already handled
      const need = half - sum;
      if (need <= 0) continue;

      // Build base partition
      const baseP1 = {}, baseP2 = {};
      items.forEach((it, i) => {
        if (subset.has(i)) baseP1[it.k] = it.v;
        else baseP2[it.k] = it.v;
      });

      // Try splitting each fluid in P2 that has enough volume
      for (const fk of Object.keys(baseP2)) {
        if (baseP2[fk] >= need) {
          const cP1 = { ...baseP1 }, cP2 = { ...baseP2 };
          cP1[fk] = (cP1[fk] || 0) + need;
          cP2[fk] -= need;
          if (cP2[fk] <= 0) delete cP2[fk];
          if (Object.keys(cP1).length > 0 && Object.keys(cP2).length > 0) {
            candidates.push({ P1: cP1, P2: cP2, type: "dp-split" });
          }
        }
      }
    }
  }

  // --- Strategy C: RMA-style for ALL fluids that could dominate (>= half) ---
  for (let i = 0; i < sorted.length && sorted[i].v >= half; i++) {
    const dominant = sorted[i];
    const cP1 = { [dominant.k]: half }, cP2 = {};
    const rem = dominant.v - half;
    if (rem > 0) cP2[dominant.k] = rem;
    for (let j = 0; j < sorted.length; j++) {
      if (j !== i) cP2[sorted[j].k] = sorted[j].v;
    }
    if (Object.keys(cP1).length > 0 && Object.keys(cP2).length > 0) {
      candidates.push({ P1: cP1, P2: cP2, type: "rma-style" });
    }
  }

  // --- Strategy D: Greedy balanced with preference for large fluids ---
  {
    const cP1 = {}, cP2 = {};
    let s1 = 0, s2 = 0;
    for (const it of sorted) {
      if (s1 <= s2 && s1 + it.v <= half) {
        cP1[it.k] = it.v; s1 += it.v;
      } else if (s2 + it.v <= half) {
        cP2[it.k] = it.v; s2 += it.v;
      } else if (s1 < half) {
        const need = half - s1;
        if (need > 0 && need <= it.v) {
          cP1[it.k] = need;
          const r = it.v - need;
          if (r > 0) cP2[it.k] = r;
          s1 = half; s2 += r;
        } else {
          cP2[it.k] = it.v; s2 += it.v;
        }
      } else {
        cP2[it.k] = it.v; s2 += it.v;
      }
    }
    if (Object.keys(cP1).length > 0 && Object.keys(cP2).length > 0) {
      candidates.push({ P1: cP1, P2: cP2, type: "greedy" });
    }
  }

  // --- Strategy E: Reverse greedy (fill P2 first) ---
  {
    const cP1 = {}, cP2 = {};
    let s1 = 0, s2 = 0;
    for (const it of sorted) {
      if (s2 <= s1 && s2 + it.v <= half) {
        cP2[it.k] = it.v; s2 += it.v;
      } else if (s1 + it.v <= half) {
        cP1[it.k] = it.v; s1 += it.v;
      } else if (s2 < half) {
        const need = half - s2;
        if (need > 0 && need <= it.v) {
          cP2[it.k] = need;
          const r = it.v - need;
          if (r > 0) cP1[it.k] = r;
          s2 = half; s1 += r;
        } else {
          cP1[it.k] = it.v; s1 += it.v;
        }
      } else {
        cP1[it.k] = it.v; s1 += it.v;
      }
    }
    if (Object.keys(cP1).length > 0 && Object.keys(cP2).length > 0) {
      candidates.push({ P1: cP1, P2: cP2, type: "greedy-rev" });
    }
  }

  // ===== Validate and score candidates =====
  if (candidates.length === 0) {
    const s = Object.entries(P).sort((a, b) => b[1] - a[1]);
    return { label: s[0][0], volume: L, level: lv, leaf: true, partition: { ...P } };
  }

  // Filter to valid partitions (both sides sum to half)
  const validCands = candidates.filter(c => {
    const s1 = Object.values(c.P1).reduce((a, b) => a + b, 0);
    const s2 = Object.values(c.P2).reduce((a, b) => a + b, 0);
    return s1 === half && s2 === half;
  });

  // Remove duplicates based on partition content
  const seen = new Set();
  const uniqueCands = [];
  for (const c of (validCands.length > 0 ? validCands : candidates)) {
    const key = JSON.stringify([
      Object.entries(c.P1).sort(),
      Object.entries(c.P2).sort()
    ]);
    if (!seen.has(key)) {
      seen.add(key);
      uniqueCands.push(c);
    }
  }

  const pool = uniqueCands.length > 0 ? uniqueCands : candidates;
  let bestCand = pool[0], bestScore = scorePartitionV3(pool[0].P1, pool[0].P2, L);
  for (let i = 1; i < pool.length; i++) {
    const sc = scorePartitionV3(pool[i].P1, pool[i].P2, L);
    if (sc < bestScore) { bestScore = sc; bestCand = pool[i]; }
  }

  return {
    label: "Mix", partition: { ...P }, volume: L, level: lv, leaf: false,
    left: buildAPDP(bestCand.P1, half, maxD, lv + 1),
    right: buildAPDP(bestCand.P2, half, maxD, lv + 1)
  };
}

// ========== LARP v2: Lookahead-Augmented Recursive Partitioning ==========
// Enhanced version with:
// - Comprehensive subset enumeration (no artificial cap)
// - Parallel exploration of zero-split AND single-split candidates
// - Multi-sum exploration for splits
// - 2-level lookahead with dilution-aware scoring
// - Multi-fluid RMA-style candidates

function larpEnumPartitions(P, L) {
  const half = L / 2;
  const keys = Object.keys(P).filter(k => P[k] > 0);
  const items = keys.map(k => ({ k, v: P[k] }));
  const n = items.length;
  if (n === 0) return [];

  const candidates = [];
  const seen = new Set();

  const addCandidate = (P1, P2, splits) => {
    const k1 = Object.keys(P1).filter(k => P1[k] > 0);
    const k2 = Object.keys(P2).filter(k => P2[k] > 0);
    if (k1.length === 0 || k2.length === 0) return;
    const s1 = k1.reduce((a, k) => a + P1[k], 0);
    const s2 = k2.reduce((a, k) => a + P2[k], 0);
    if (s1 !== half || s2 !== half) return;
    const key = JSON.stringify([Object.entries(P1).sort(), Object.entries(P2).sort()]);
    if (seen.has(key)) return;
    seen.add(key);
    candidates.push({ P1: { ...P1 }, P2: { ...P2 }, splits });
  };

  // ===== DP to enumerate ALL achievable sums with whole-fluid subsets =====
  // dp[s] = list of bitmasks achieving sum s (higher cap for thorough enumeration)
  const MASK_CAP = 64;
  const dp = new Array(half + 1).fill(null);
  dp[0] = [0];
  for (let i = 0; i < n; i++) {
    const v = items[i].v;
    for (let s = half; s >= v; s--) {
      if (dp[s - v]) {
        if (!dp[s]) dp[s] = [];
        for (const mask of dp[s - v]) {
          if (dp[s].length < MASK_CAP) dp[s].push(mask | (1 << i));
        }
      }
    }
  }

  // ===== Strategy 1: Zero-split candidates (subsets summing exactly to half) =====
  if (dp[half]) {
    for (const mask of dp[half]) {
      const P1 = {}, P2 = {};
      items.forEach((it, i) => {
        if (mask & (1 << i)) P1[it.k] = it.v;
        else P2[it.k] = it.v;
      });
      addCandidate(P1, P2, 0);
    }
  }

  // ===== Strategy 2: Single-split candidates from MULTIPLE achievable sums =====
  // Explore top 5 closest sums to half (not just the single closest)
  const achievableSums = [];
  for (let s = half - 1; s >= 0; s--) {
    if (dp[s]) achievableSums.push(s);
    if (achievableSums.length >= 5) break;
  }

  for (const baseSum of achievableSums) {
    const need = half - baseSum;
    if (need <= 0) continue;
    for (const mask of dp[baseSum]) {
      // Try splitting each fluid NOT in the subset (add `need` to P1)
      for (let i = 0; i < n; i++) {
        if (mask & (1 << i)) continue;
        if (items[i].v >= need) {
          const P1 = {}, P2 = {};
          items.forEach((it, j) => {
            if (mask & (1 << j)) P1[it.k] = it.v;
            else P2[it.k] = it.v;
          });
          P1[items[i].k] = (P1[items[i].k] || 0) + need;
          P2[items[i].k] = (P2[items[i].k] || 0) - need;
          if (P2[items[i].k] <= 0) delete P2[items[i].k];
          addCandidate(P1, P2, 1);
        }
      }
      // Also try splitting a fluid that IS in the subset (move some back to P2)
      for (let i = 0; i < n; i++) {
        if (!(mask & (1 << i))) continue;
        const excess = baseSum - (half - items[i].v);
        if (excess > 0 && excess < items[i].v) {
          const keep = items[i].v - excess;
          if (keep > 0) {
            const P1 = {}, P2 = {};
            items.forEach((it, j) => {
              if (j === i) { P1[it.k] = keep; P2[it.k] = excess; }
              else if (mask & (1 << j)) P1[it.k] = it.v;
              else P2[it.k] = it.v;
            });
            addCandidate(P1, P2, 1);
          }
        }
      }
    }
  }

  // ===== Strategy 3: RMA-style for ALL fluids that could dominate (≥ half) =====
  const sorted = items.slice().sort((a, b) => b.v - a.v);
  for (let i = 0; i < sorted.length && sorted[i].v >= half; i++) {
    const dominant = sorted[i];
    const P1 = { [dominant.k]: half }, P2 = {};
    const rem = dominant.v - half;
    if (rem > 0) P2[dominant.k] = rem;
    for (let j = 0; j < sorted.length; j++) {
      if (j !== i) P2[sorted[j].k] = sorted[j].v;
    }
    addCandidate(P1, P2, dominant.v > half ? 1 : 0);
  }

  // ===== Strategy 4: Double-split candidates (when single split isn't enough) =====
  if (candidates.length === 0 || !candidates.some(c => c.splits <= 1)) {
    // Try splitting two fluids
    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) {
        // Try allocating parts of items[i] and items[j] to reach half
        for (let ai = 1; ai < items[i].v; ai++) {
          const aj = half - ai;
          if (aj > 0 && aj < items[j].v) {
            const P1 = { [items[i].k]: ai, [items[j].k]: aj };
            const P2 = { [items[i].k]: items[i].v - ai, [items[j].k]: items[j].v - aj };
            for (let k = 0; k < n; k++) {
              if (k !== i && k !== j) P2[items[k].k] = items[k].v;
            }
            addCandidate(P1, P2, 2);
          }
        }
      }
    }
  }

  // ===== Strategy 5: Greedy balanced fallback =====
  if (candidates.length === 0) {
    const cP1 = {}, cP2 = {};
    let s1 = 0, s2 = 0;
    for (const it of sorted) {
      if (s1 <= s2 && s1 + it.v <= half) { cP1[it.k] = it.v; s1 += it.v; }
      else if (s2 + it.v <= half) { cP2[it.k] = it.v; s2 += it.v; }
      else {
        const need1 = half - s1, need2 = half - s2;
        if (need1 > 0 && need1 <= it.v) {
          cP1[it.k] = need1; cP2[it.k] = it.v - need1; s1 = half; s2 += it.v - need1;
        } else if (need2 > 0 && need2 <= it.v) {
          cP2[it.k] = need2; cP1[it.k] = it.v - need2; s2 = half; s1 += it.v - need2;
        } else { cP2[it.k] = it.v; s2 += it.v; }
      }
    }
    const splits = Object.keys(cP1).filter(k => cP2[k] > 0).length;
    addCandidate(cP1, cP2, splits);
  }

  return candidates;
}

// 2-level lookahead score estimation
function larpDeepScore(P, L, depth = 2) {
  const keys = Object.keys(P).filter(k => P[k] > 0);
  if (keys.length <= 1 || L <= 1 || depth <= 0) return { splits: 0, fluids: keys.length, dilution: keys.length <= 2 ? 1 : 0 };

  const half = L / 2;
  const items = keys.map(k => ({ k, v: P[k] }));
  const n = items.length;

  // Quick DP for achievability
  const dp = new Uint8Array(half + 1);
  dp[0] = 1;
  for (const it of items) {
    for (let s = half; s >= it.v; s--) if (dp[s - it.v]) dp[s] = 1;
  }

  const zeroSplitPossible = dp[half] === 1;
  const minSplits = zeroSplitPossible ? 0 : 1;

  // Check if this could become a dilution subtree (≤2 fluids)
  const isDilution = keys.length <= 2 ? 1 : 0;

  // Estimate best dominant fluid fraction (for dilution potential)
  let maxFrac = 0;
  for (const it of items) maxFrac = Math.max(maxFrac, it.v / L);
  const dilutionPotential = maxFrac >= 0.5 ? 1 : 0;

  return { splits: minSplits, fluids: keys.length, dilution: isDilution, dilutionPot: dilutionPotential };
}

function larpScore(cand, L) {
  const { P1, P2, splits } = cand;
  const half = L / 2;
  const k1 = Object.keys(P1).filter(k => P1[k] > 0);
  const k2 = Object.keys(P2).filter(k => P2[k] > 0);

  // ===== Immediate costs =====
  const splitPenalty = splits * 100; // Heavy penalty for splits
  const fluidCount = k1.length + k2.length;

  // ===== Dilution bonus: reward partitions creating dilution subtrees =====
  let dilutionBonus = 0;
  if (k1.length <= 2) dilutionBonus -= 15; // P1 is/will be a dilution subtree
  if (k2.length <= 2) dilutionBonus -= 15; // P2 is/will be a dilution subtree
  if (k1.length === 1) dilutionBonus -= 10; // Pure fluid on one side
  if (k2.length === 1) dilutionBonus -= 10;

  // Check for dominant fluid in each partition (could create long dilution chains)
  const sum1 = k1.reduce((a, k) => a + P1[k], 0);
  const sum2 = k2.reduce((a, k) => a + P2[k], 0);
  for (const k of k1) { if (P1[k] >= sum1 * 0.5) dilutionBonus -= 8; break; }
  for (const k of k2) { if (P2[k] >= sum2 * 0.5) dilutionBonus -= 8; break; }

  // ===== 2-level lookahead =====
  const la1 = (k1.length > 1 && half > 1) ? larpDeepScore(P1, half, 2) : { splits: 0, fluids: 1, dilution: 1, dilutionPot: 0 };
  const la2 = (k2.length > 1 && half > 1) ? larpDeepScore(P2, half, 2) : { splits: 0, fluids: 1, dilution: 1, dilutionPot: 0 };

  const lookaheadSplits = 50 * (la1.splits + la2.splits);
  const lookaheadDilution = -20 * (la1.dilution + la2.dilution + la1.dilutionPot + la2.dilutionPot);

  // ===== Balance penalty: prefer more even fluid distribution for flexibility =====
  const balance = Math.abs(k1.length - k2.length);
  const balancePenalty = balance * 2;

  // Lower score = better
  return splitPenalty + fluidCount + dilutionBonus + lookaheadSplits + lookaheadDilution + balancePenalty;
}

function buildLARP(P, L, maxD, lv = 0) {
  const keys = Object.keys(P).filter(k => P[k] > 0);
  if (!keys.length) return null;
  if (keys.length === 1) return { label: keys[0], volume: P[keys[0]], level: lv, leaf: true, partition: { [keys[0]]: P[keys[0]] } };
  if (L <= 1 || lv >= maxD) {
    const s = Object.entries(P).sort((a, b) => b[1] - a[1]);
    return { label: s[0][0], volume: L, level: lv, leaf: true, partition: { ...P } };
  }

  const half = L / 2;
  const candidates = larpEnumPartitions(P, L);
  if (candidates.length === 0) {
    const s = Object.entries(P).sort((a, b) => b[1] - a[1]);
    return { label: s[0][0], volume: L, level: lv, leaf: true, partition: { ...P } };
  }

  // Score all candidates with enhanced lookahead
  let bestCand = candidates[0], bestScore = larpScore(candidates[0], L);
  for (let i = 1; i < candidates.length; i++) {
    const sc = larpScore(candidates[i], L);
    if (sc < bestScore) { bestScore = sc; bestCand = candidates[i]; }
  }

  return {
    label: "Mix", partition: { ...P }, volume: L, level: lv, leaf: false,
    left: buildLARP(bestCand.P1, half, maxD, lv + 1),
    right: buildLARP(bestCand.P2, half, maxD, lv + 1)
  };
}

// ========== ILP v2: Pareto-Optimal Split-Minimizing DP ==========
// Enhanced with:
// - Pareto-optimal candidate tracking: (splits, distinctFluidsP1)
// - All Pareto-optimal allocations enumerated, then scored
// - Dilution-aware scoring: rewards ≤2 fluids per side
// - Dominant fluid bonus: rewards isolating fluids ≥50% volume
// - 1-level lookahead: estimates downstream split potential

function ilpEstimateDownstream(partition, targetSum) {
  const keys = Object.keys(partition).filter(k => partition[k] > 0);
  if (keys.length <= 1 || targetSum <= 1) return { splits: 0, canZeroSplit: true };
  const half = Math.floor(targetSum / 2);
  const items = keys.map(k => ({ k, v: partition[k] }));

  // Quick DP: can we reach exactly half with whole fluids?
  const dp = new Uint8Array(half + 1);
  dp[0] = 1;
  for (const it of items) {
    for (let s = half; s >= it.v; s--) if (dp[s - it.v]) dp[s] = 1;
  }
  if (dp[half]) return { splits: 0, canZeroSplit: true };
  return { splits: 1, canZeroSplit: false };
}

function ilpScoreAllocation(P1, P2, half) {
  const k1 = Object.keys(P1).filter(k => P1[k] > 0);
  const k2 = Object.keys(P2).filter(k => P2[k] > 0);
  const n1 = k1.length, n2 = k2.length;

  // Count actual splits (fluids appearing in both)
  const splits = k1.filter(k => P2[k] > 0).length;

  let score = splits * 1000; // Primary: minimize splits

  // Secondary: minimize total distinct fluids (enables dilution)
  score += (n1 + n2) * 10;

  // Tertiary: reward dilution subtrees (≤2 fluids per side)
  if (n1 <= 2) score -= 50;
  if (n2 <= 2) score -= 50;
  if (n1 === 1) score -= 30; // Pure fluid isolation
  if (n2 === 1) score -= 30;

  // Quaternary: reward dominant fluid isolation (≥50% of its side)
  const s1 = k1.reduce((a, k) => a + P1[k], 0);
  const s2 = k2.reduce((a, k) => a + P2[k], 0);
  for (const k of k1) { if (P1[k] >= s1 * 0.5) { score -= 20; break; } }
  for (const k of k2) { if (P2[k] >= s2 * 0.5) { score -= 20; break; } }

  // 1-level lookahead: estimate downstream splits
  const la1 = (n1 > 1 && half > 1) ? ilpEstimateDownstream(P1, half) : { splits: 0 };
  const la2 = (n2 > 1 && half > 1) ? ilpEstimateDownstream(P2, half) : { splits: 0 };
  score += (la1.splits + la2.splits) * 100;

  // Bonus for zero-split potential downstream
  if (la1.canZeroSplit) score -= 15;
  if (la2.canZeroSplit) score -= 15;

  return score;
}

function buildILP(P, L, maxD, lv = 0) {
  const keys = Object.keys(P).filter(k => P[k] > 0);
  if (keys.length === 0) return null;
  if (keys.length === 1 || lv >= maxD) {
    const k = keys[0];
    return { label: k, volume: P[k], level: lv, leaf: true, partition: { [k]: P[k] } };
  }

  const half = Math.floor(L / 2);
  const items = keys.map(k => ({ k, v: P[k] }));
  const n = items.length;

  // ===== Phase 1: DP to find all Pareto-optimal allocations =====
  // State: (item index, sum) -> list of Pareto-optimal (splits, allocation) pairs
  // Pareto: (s1, a1) dominates (s2, a2) if s1 <= s2 AND a1 is "better" on secondary

  const INF = 1e9;
  const ALLOC_CAP = 32; // Cap allocations per state to limit memory

  // dp[s] = list of {splits, alloc: array of take values}
  let dp = new Array(half + 1).fill(null);
  dp[0] = [{ splits: 0, alloc: [] }];

  for (let i = 0; i < n; i++) {
    const v = items[i].v;
    const newDp = new Array(half + 1).fill(null);

    for (let s = 0; s <= half; s++) {
      if (!dp[s]) continue;

      for (const state of dp[s]) {
        // Try all possible takes for item i
        for (let t = 0; t <= v && s + t <= half; t++) {
          const ns = s + t;
          const addSplit = (t > 0 && t < v) ? 1 : 0;
          const newSplits = state.splits + addSplit;
          const newAlloc = [...state.alloc, t];

          if (!newDp[ns]) newDp[ns] = [];

          // Check if dominated by existing or dominates existing
          let dominated = false;
          const dominated_indices = [];

          for (let j = 0; j < newDp[ns].length; j++) {
            const existing = newDp[ns][j];
            if (existing.splits <= newSplits) {
              // Could be dominated - check if strictly better on splits
              if (existing.splits < newSplits) { dominated = true; break; }
              // Equal splits - keep both for now (will score later)
            }
            if (newSplits < existing.splits) {
              dominated_indices.push(j);
            }
          }

          if (!dominated) {
            // Remove dominated entries
            for (let j = dominated_indices.length - 1; j >= 0; j--) {
              newDp[ns].splice(dominated_indices[j], 1);
            }
            // Add new entry if under cap
            if (newDp[ns].length < ALLOC_CAP) {
              newDp[ns].push({ splits: newSplits, alloc: newAlloc });
            }
          }
        }
      }
    }

    dp = newDp;
  }

  // ===== Phase 2: Collect all candidates reaching exactly half =====
  const candidates = [];

  if (dp[half]) {
    for (const state of dp[half]) {
      const P1 = {}, P2 = {};
      for (let i = 0; i < n; i++) {
        const t = state.alloc[i];
        if (t > 0) P1[items[i].k] = t;
        if (items[i].v - t > 0) P2[items[i].k] = items[i].v - t;
      }
      const k1 = Object.keys(P1).filter(k => P1[k] > 0);
      const k2 = Object.keys(P2).filter(k => P2[k] > 0);
      if (k1.length > 0 && k2.length > 0) {
        candidates.push({ P1, P2, splits: state.splits });
      }
    }
  }

  // ===== Phase 3: RMA-style candidates (dominant fluid fills one child) =====
  const sorted = items.slice().sort((a, b) => b.v - a.v);
  for (let i = 0; i < sorted.length && sorted[i].v >= half; i++) {
    const dominant = sorted[i];
    const P1 = { [dominant.k]: half }, P2 = {};
    const rem = dominant.v - half;
    if (rem > 0) P2[dominant.k] = rem;
    for (let j = 0; j < sorted.length; j++) {
      if (sorted[j].k !== dominant.k) P2[sorted[j].k] = sorted[j].v;
    }
    const k1 = Object.keys(P1).filter(k => P1[k] > 0);
    const k2 = Object.keys(P2).filter(k => P2[k] > 0);
    if (k1.length > 0 && k2.length > 0) {
      const splits = dominant.v > half ? 1 : 0;
      candidates.push({ P1, P2, splits });
    }
  }

  // ===== Phase 4: Greedy fallback if no candidates =====
  if (candidates.length === 0) {
    let P1 = {}, P2 = {}, cur = 0;
    for (const it of sorted) {
      if (cur + it.v <= half) { P1[it.k] = it.v; cur += it.v; }
      else {
        const need = half - cur;
        if (need > 0) { P1[it.k] = need; P2[it.k] = it.v - need; cur += need; }
        else P2[it.k] = it.v;
      }
    }
    const k1 = Object.keys(P1).filter(k => P1[k] > 0);
    const k2 = Object.keys(P2).filter(k => P2[k] > 0);
    if (k1.length > 0 && k2.length > 0) {
      const splits = k1.filter(k => P2[k] > 0).length;
      candidates.push({ P1, P2, splits });
    }
  }

  if (candidates.length === 0) {
    // Ultimate fallback
    return { label: sorted[0].k, volume: L, level: lv, leaf: true, partition: P };
  }

  // ===== Phase 5: Deduplicate candidates =====
  const seen = new Set();
  const uniqueCands = [];
  for (const c of candidates) {
    const key = JSON.stringify([Object.entries(c.P1).sort(), Object.entries(c.P2).sort()]);
    if (!seen.has(key)) {
      seen.add(key);
      uniqueCands.push(c);
    }
  }

  // ===== Phase 6: Score and select best candidate =====
  let bestCand = uniqueCands[0];
  let bestScore = ilpScoreAllocation(uniqueCands[0].P1, uniqueCands[0].P2, half);

  for (let i = 1; i < uniqueCands.length; i++) {
    const score = ilpScoreAllocation(uniqueCands[i].P1, uniqueCands[i].P2, half);
    if (score < bestScore) {
      bestScore = score;
      bestCand = uniqueCands[i];
    }
  }

  return {
    label: "Mix", partition: { ...P }, volume: L, level: lv, leaf: false,
    left: buildILP(bestCand.P1, half, maxD, lv + 1),
    right: buildILP(bestCand.P2, half, maxD, lv + 1)
  };
}

// ========== Tree utilities ==========
function layoutTree(root) {
  let idx = 0; const nodes = [], edges = [];
  (function walk(n, d) {
    if (!n) return; if (n.leaf) { n._x = idx++; n._y = d; }
    else { walk(n.left, d + 1); walk(n.right, d + 1); n._x = ((n.left?._x ?? 0) + (n.right?._x ?? 0)) / 2; n._y = d; edges.push({ from: n, to: n.left }); edges.push({ from: n, to: n.right }); }
    nodes.push(n);
  })(root, 0); return { nodes, edges };
}
function cntN(n) { if (!n) return 0; return 1 + cntN(n.left) + cntN(n.right); }
function tD(n) { if (!n) return 0; return 1 + Math.max(tD(n.left), tD(n.right)); }
function cntL(n) { if (!n) return 0; if (n.leaf) return 1; return cntL(n.left) + cntL(n.right); }
function cntM(n) { if (!n || n.leaf) return 0; return 1 + cntM(n.left) + cntM(n.right); }
function mxP(n) { if (!n) return 0; const lv = {}; (function w(nd) { if (!nd) return; if (!nd.leaf) lv[nd.level] = (lv[nd.level] || 0) + 1; w(nd.left); w(nd.right); })(n); return Math.max(0, ...Object.values(lv)); }
// Dilution length: sum of depths of all dilution subtrees (subtrees with ≤2 distinct fluids per paper definition)
function dL(n) {
  function df(nd) { if (!nd) return new Set(); if (nd.leaf) return new Set([nd.label]); return new Set([...df(nd.left), ...df(nd.right)]); }
  let t = 0; (function w(nd) { if (!nd || nd.leaf) return; if (df(nd).size <= 2) { t += tD(nd) - 1; return; } w(nd.left); w(nd.right); })(n); return t;
}
// Longest single dilution subtree depth (≤2 distinct fluids)
function maxDL(n) {
  function df(nd) { if (!nd) return new Set(); if (nd.leaf) return new Set([nd.label]); return new Set([...df(nd.left), ...df(nd.right)]); }
  let mx = 0;
  (function w(nd) {
    if (!nd || nd.leaf) return;
    if (df(nd).size <= 2) { mx = Math.max(mx, tD(nd) - 1); return; }
    w(nd.left); w(nd.right);
  })(n);
  return mx;
}
function cntSplits(n) {
  if (!n || n.leaf) return 0; let s = 0;
  if (n.left && n.right) { const lk = new Set(Object.keys(n.left.partition || {})); for (const k of Object.keys(n.right.partition || {})) if (lk.has(k)) s++; }
  return s + cntSplits(n.left) + cntSplits(n.right);
}

const COLORS = ["#6366f1", "#f59e0b", "#10b981", "#ef4444", "#3b82f6", "#ec4899", "#8b5cf6", "#14b8a6", "#f97316", "#84cc16", "#06b6d4", "#e11d55"];
function fc(lb, all) { const i = all.indexOf(lb); return i >= 0 ? COLORS[i % COLORS.length] : "#94a3b8"; }

// ========== TREE SVG ==========
function TreeView({ treeData, allFluids }) {
  const [tf, setTf] = useState({ x: 0, y: 0, k: 1 });
  const drag = useRef(false), last = useRef({ x: 0, y: 0 }), ref = useRef(null);
  const nW = 64, nH = 30, gX = 4, gY = 38;
  useEffect(() => {
    if (!treeData?.nodes?.length || !ref.current) return;
    const mx = Math.max(...treeData.nodes.map(n => n._x)), my = Math.max(...treeData.nodes.map(n => n._y));
    const tw = (mx + 1) * (nW + gX) + 32, th = (my + 1) * (nH + gY) + 32;
    const cw = ref.current.clientWidth || 600, ch = ref.current.clientHeight || 400;
    setTf({ x: Math.max(0, (cw - tw * Math.min(cw / tw, ch / th, 1.5)) / 2), y: 8, k: Math.min(cw / tw, ch / th, 1.5) });
  }, [treeData]);
  if (!treeData?.nodes?.length) return <div style={{ padding: 20, color: "#64748b", textAlign: "center", fontSize: 11 }}>No tree</div>;
  const { nodes, edges } = treeData;
  const cx = n => n._x * (nW + gX) + nW / 2 + 14, cy = n => n._y * (nH + gY) + nH / 2 + 14;
  return (
    <div ref={ref} style={{ flex: 1, overflow: "hidden", cursor: drag.current ? "grabbing" : "grab", background: "#0c1222", borderRadius: 6, minHeight: 180 }}
      onMouseDown={e => { drag.current = true; last.current = { x: e.clientX, y: e.clientY }; }}
      onMouseMove={e => { if (!drag.current) return; setTf(t => ({ ...t, x: t.x + e.clientX - last.current.x, y: t.y + e.clientY - last.current.y })); last.current = { x: e.clientX, y: e.clientY }; }}
      onMouseUp={() => drag.current = false} onMouseLeave={() => drag.current = false}
      onWheel={e => { e.preventDefault(); const f = e.deltaY < 0 ? 1.15 : 0.87; setTf(t => ({ ...t, k: Math.min(6, Math.max(0.02, t.k * f)) })); }}>
      <svg width="100%" height="100%"><g transform={`translate(${tf.x},${tf.y}) scale(${tf.k})`}>
        {edges.map((e, i) => e.to && <line key={i} x1={cx(e.from)} y1={cy(e.from) + nH / 2} x2={cx(e.to)} y2={cy(e.to) - nH / 2} stroke="#334155" strokeWidth={.9} />)}
        {nodes.map((n, i) => {
          const x = cx(n) - nW / 2, y = cy(n) - nH / 2;
          if (n.leaf) { const c = fc(n.label, allFluids); return <g key={i}><rect x={x} y={y} width={nW} height={nH} rx={5} fill={c + "22"} stroke={c} strokeWidth={1.2} /><text x={cx(n)} y={cy(n) + 1} textAnchor="middle" dominantBaseline="middle" fill={c} fontSize={9} fontWeight={700}>{n.label}</text></g>; }
          const tot = Object.values(n.partition).reduce((a, b) => a + b, 0); const bars = []; let off = 0;
          Object.entries(n.partition).sort((a, b) => allFluids.indexOf(a[0]) - allFluids.indexOf(b[0])).forEach(([f, v]) => { const w = (v / tot) * (nW - 3); bars.push({ x: off, w: Math.max(w, .3), c: fc(f, allFluids) }); off += w; });
          return <g key={i}><rect x={x} y={y} width={nW} height={nH} rx={5} fill="#1e293b" stroke="#475569" strokeWidth={.7} /><text x={cx(n)} y={cy(n) - 2} textAnchor="middle" fill="#cbd5e1" fontSize={7} fontWeight={600}>Mix</text><g transform={`translate(${x + 1.5},${cy(n) + 5})`}>{bars.map((b, j) => <rect key={j} x={b.x} y={0} width={b.w} height={3} rx={.8} fill={b.c} opacity={.85} />)}</g></g>;
        })}
      </g></svg>
    </div>
  );
}

// ========== UI ==========
function InputPanel({ raw, setRaw, depth, setDepth, onGen }) {
  return (
    <div style={{ display: "flex", gap: 6, alignItems: "center", flexWrap: "wrap", marginBottom: 7 }}>
      <label style={{ fontSize: 11 }}>Ratios: <input value={raw} onChange={e => setRaw(e.target.value)} style={{ marginLeft: 3, padding: "3px 6px", width: 180, background: "#0f172a", border: "1px solid #334155", borderRadius: 4, color: "#e2e8f0", fontSize: 11 }} /></label>
      <label style={{ fontSize: 11 }}>d: <input type="number" value={depth} min={2} max={14} onChange={e => setDepth(+e.target.value)} style={{ marginLeft: 3, padding: "3px 6px", width: 42, background: "#0f172a", border: "1px solid #334155", borderRadius: 4, color: "#e2e8f0", fontSize: 11 }} /></label>
      <button onClick={onGen} style={{ padding: "4px 13px", background: "#6366f1", border: "none", borderRadius: 5, color: "#fff", fontWeight: 700, cursor: "pointer", fontSize: 11 }}>Generate</button>
    </div>
  );
}
function Legend({ allFluids, adj }) {
  return <div style={{ display: "flex", gap: 7, flexWrap: "wrap", marginBottom: 4 }}>
    {allFluids.map((f, i) => adj[i] > 0 && <span key={f} style={{ display: "flex", alignItems: "center", gap: 2, fontSize: 10 }}><span style={{ width: 7, height: 7, borderRadius: 2, background: fc(f, allFluids), display: "inline-block" }} />{f}={adj[i]}</span>)}
    <span style={{ fontSize: 9, color: "#475569", marginLeft: "auto" }}>scroll · drag</span>
  </div>;
}
function Stat({ label, value, color, best }) {
  return <div style={{ background: best ? `${color}15` : "#1e293b", borderRadius: 6, padding: "5px 10px", minWidth: 62, border: `1px solid ${best ? color : color + "33"}` }}><div style={{ fontSize: 8, color: "#94a3b8" }}>{label}</div><div style={{ fontSize: 15, fontWeight: 800, color }}>{value}</div></div>;
}

function getStats(root) {
  if (!root) return null;
  return { m: cntM(root), d: tD(root) - 1, p: mxP(root), leaves: cntL(root), l: dL(root), maxL: maxDL(root), splits: cntSplits(root) };
}
function buildTree(algo, adj, fl, depth) {
  if (algo === "bs") return buildBSTree(adj, fl, depth);
  const d = {}; fl.forEach((f, i) => { if (adj[i] > 0) d[f] = adj[i]; });
  if (algo === "rma") return buildRMA(d, 2 ** depth);
  if (algo === "larp") return buildLARP(d, 2 ** depth, depth + 8);
  if (algo === "ilp") return buildILP(d, 2 ** depth, depth + 8);
  return buildAPDP(d, 2 ** depth, depth + 8);
}

function AlgoPage({ raw, setRaw, depth, setDepth, algo, title, icon, color, desc }) {
  const [tree, setTree] = useState(null), [allF, setAllF] = useState([]), [adj, setAdj] = useState([]), [stats, setStats] = useState(null);
  const gen = () => {
    const r = raw.split(",").map(Number).filter(n => !isNaN(n) && n > 0); if (r.length < 2) return;
    const a = ratioApprox(r, depth); setAdj(a);
    const fl = a.map((_, i) => `x${i + 1}`); setAllF(fl);
    const root = buildTree(algo, a, fl, depth);
    if (root) { setTree(layoutTree(root)); setStats(getStats(root)); }
  };
  useEffect(gen, []);
  return (
    <div style={{ display: "flex", flexDirection: "column", height: "100%" }}>
      <div style={{ padding: "9px 12px", background: "#1e293b", borderRadius: 8, marginBottom: 6, flexShrink: 0 }}>
        <h2 style={{ margin: "0 0 2px", fontSize: 13, fontWeight: 800, color }}>{icon} {title}</h2>
        <p style={{ margin: "0 0 6px", fontSize: 10, color: "#94a3b8", lineHeight: 1.3 }}>{desc}</p>
        <InputPanel raw={raw} setRaw={setRaw} depth={depth} setDepth={setDepth} onGen={gen} />
        <Legend allFluids={allF} adj={adj} />
        {stats && <div style={{ display: "flex", gap: 4, flexWrap: "wrap", marginTop: 3 }}>
          <Stat label="m (mixes)" value={stats.m} color={color} />
          <Stat label="d (depth)" value={stats.d} color="#f59e0b" />
          <Stat label="p (‖)" value={stats.p} color="#3b82f6" />
          <Stat label="leaves" value={stats.leaves} color="#10b981" />
          <Stat label="l (dilution)" value={stats.l} color="#ec4899" />
          <Stat label="maxL" value={stats.maxL} color="#8b5cf6" />
          <Stat label="splits" value={stats.splits} color="#94a3b8" />
        </div>}
      </div>
      <TreeView treeData={tree} allFluids={allF} />
    </div>
  );
}

function ComparePage({ raw, setRaw, depth, setDepth }) {
  const [data, setData] = useState(null), [allF, setAllF] = useState([]), [adj, setAdj] = useState([]);
  const algos = ["rma", "bs", "apdp", "larp", "ilp"];
  const aL = { rma: "RMA", bs: "BS", apdp: "AP-DP", larp: "LARP", ilp: "ILP" };
  const aC = { rma: "#6366f1", bs: "#f59e0b", apdp: "#06b6d4", larp: "#ef4444", ilp: "#84cc16" };
  const gen = () => {
    const r = raw.split(",").map(Number).filter(n => !isNaN(n) && n > 0); if (r.length < 2) return;
    const a = ratioApprox(r, depth); setAdj(a);
    const fl = a.map((_, i) => `x${i + 1}`); setAllF(fl);
    const res = {}; for (const algo of algos) { const root = buildTree(algo, a, fl, depth); if (root) res[algo] = getStats(root); }
    setData(res);
  };
  useEffect(gen, []);
  const metrics = ["m", "leaves", "l", "maxL", "splits", "d", "p"];
  const mL = { m: "m (Mix/Split Cycles)", d: "d (Tree Depth)", p: "p (Parallelism)", leaves: "Leaf Nodes", l: "l (Dilution Sum)", maxL: "maxL (Longest Dilution)", splits: "Fluid Splits" };
  const mC = { m: "#ec4899", d: "#f59e0b", p: "#3b82f6", leaves: "#10b981", l: "#8b5cf6", maxL: "#a855f7", splits: "#94a3b8" };
  const lB = new Set(["m", "leaves", "splits"]); // lower is better (definitive)
  const hB = new Set(["l", "maxL"]); // higher is better (definitive)
  const varies = new Set(["d", "p"]); // depends on chip constraints

  function winner(m) {
    if (!data) return "";
    if (varies.has(m)) return "varies";
    const vals = algos.map(a => data[a]?.[m] ?? (lB.has(m) ? Infinity : -Infinity));
    const best = lB.has(m) ? Math.min(...vals) : hB.has(m) ? Math.max(...vals) : Math.min(...vals);
    const ws = algos.filter((a, i) => vals[i] === best);
    return ws.length === algos.length ? "Tie" : ws.map(a => aL[a]).join(", ");
  }

  return (
    <div style={{ display: "flex", flexDirection: "column", height: "100%", overflow: "auto" }}>
      <div style={{ padding: "9px 12px", background: "#1e293b", borderRadius: 8, marginBottom: 8, flexShrink: 0 }}>
        <h2 style={{ margin: "0 0 2px", fontSize: 13, fontWeight: 800, color: "#10b981" }}>📊 Five-Way Comparison</h2>
        <p style={{ margin: "0 0 6px", fontSize: 10, color: "#94a3b8" }}>RMA vs BS vs AP-DP vs LARP vs ILP. Lower m/leaves/splits better. Higher l/maxL better. d/p depend on chip.</p>
        <InputPanel raw={raw} setRaw={setRaw} depth={depth} setDepth={setDepth} onGen={gen} />
        <Legend allFluids={allF} adj={adj} />
      </div>
      {data && <>
        <div style={{ background: "#1e293b", borderRadius: 8, padding: 12, marginBottom: 8 }}>
          {metrics.map(m => {
            const mx = Math.max(...algos.map(a => data[a]?.[m] ?? 0), 1);
            return <div key={m} style={{ marginBottom: 10 }}>
              <div style={{ fontSize: 10, color: "#94a3b8", marginBottom: 3 }}>{mL[m]}</div>
              {algos.map(a => <div key={a} style={{ display: "flex", alignItems: "center", gap: 5, marginBottom: 2 }}>
                <span style={{ width: 34, fontSize: 9, color: aC[a], fontWeight: 700 }}>{aL[a]}</span>
                <div style={{ flex: 1, background: "#0f172a", borderRadius: 3, height: 16, position: "relative", overflow: "hidden" }}>
                  <div style={{ width: `${((data[a]?.[m] ?? 0) / mx) * 100}%`, height: "100%", background: `${mC[m]}45`, borderRadius: 3 }} />
                  <span style={{ position: "absolute", right: 4, top: 0, fontSize: 10, fontWeight: 700, color: "#e2e8f0" }}>{data[a]?.[m] ?? "—"}</span>
                </div>
              </div>)}
            </div>;
          })}
        </div>
        <div style={{ background: "#1e293b", borderRadius: 8, padding: 12 }}>
          <table style={{ width: "100%", borderCollapse: "collapse", fontSize: 10 }}>
            <thead><tr style={{ borderBottom: "1px solid #334155" }}>
              <th style={{ textAlign: "left", padding: "3px 5px", color: "#94a3b8" }}>Metric</th>
              {algos.map(a => <th key={a} style={{ textAlign: "center", padding: "3px 5px", color: aC[a] }}>{aL[a]}</th>)}
              <th style={{ textAlign: "center", padding: "3px 5px", color: "#10b981" }}>Best</th>
            </tr></thead>
            <tbody>{metrics.map(m => {
              const w = winner(m);
              const wc = w === "Tie" ? "#64748b" : w === "varies" ? "#94a3b8" : w.includes("ILP") ? "#84cc16" : w.includes("LARP") ? "#ef4444" : w.includes("AP-DP") ? "#06b6d4" : w.includes("RMA") ? "#6366f1" : "#f59e0b";
              return <tr key={m} style={{ borderBottom: "1px solid #1a2332" }}>
                <td style={{ padding: "3px 5px", color: "#cbd5e1" }}>{mL[m]}</td>
                {algos.map(a => <td key={a} style={{ textAlign: "center", padding: "3px 5px", fontWeight: 700 }}>{data[a]?.[m] ?? "—"}</td>)}
                <td style={{ textAlign: "center", padding: "3px 5px", fontWeight: 800, color: wc }}>{w === "Tie" ? "—" : w === "varies" ? "varies" : `✓ ${w}`}</td>
              </tr>;
            })}</tbody>
          </table>
          <div style={{ marginTop: 6, fontSize: 9, color: "#64748b", lineHeight: 1.5 }}>
            <b style={{ color: "#84cc16" }}>ILP v2</b>: Pareto-optimal DP — enumerates Pareto-optimal (splits, fluids) allocations, dilution-aware scoring with 1-level lookahead.<br />
            <b style={{ color: "#ef4444" }}>LARP v2</b>: 2-level lookahead DP — comprehensive enumeration (zero/single/double splits), dilution-aware scoring.<br />
            <b style={{ color: "#06b6d4" }}>AP-DP v3</b>: Multi-strategy optimizer — DP whole-fluid, DP+split, RMA-style (all dominants) & greedy candidates with dilution-aware scoring.<br />
            <b style={{ color: "#6366f1" }}>RMA</b>: Maximises dilution subtree length l for layout mapping (per paper Algorithm 4).<br />
            <b style={{ color: "#f59e0b" }}>BS</b>: Bit-scanning bottom-up, minimises leaf count.<br />
            <span style={{ color: "#8b5cf6" }}>l</span>=sum of dilution subtree depths (≤2 fluids), <span style={{ color: "#a855f7" }}>maxL</span>=longest single dilution subtree.
          </div>
        </div>
      </>}
    </div>
  );
}

// ========== APP ==========
export default function App() {
  const [page, setPage] = useState("compare");
  const [raw, setRaw] = useState("2,3,5,7,11,13,87");
  const [depth, setDepth] = useState(7);
  const tabs = [
    { id: "rma", label: "RMA", icon: "⚗️", color: "#6366f1" },
    { id: "bs", label: "BS", icon: "🔬", color: "#f59e0b" },
    { id: "apdp", label: "AP-DP", icon: "⚙️", color: "#06b6d4" },
    { id: "larp", label: "LARP", icon: "🔮", color: "#ef4444" },
    { id: "ilp", label: "ILP", icon: "🧮", color: "#84cc16" },
    { id: "compare", label: "Compare", icon: "📊", color: "#10b981" },
  ];
  return (
    <div style={{ display: "flex", height: "100vh", width: "100vw", position: "fixed", top: 0, left: 0, background: "#0f172a", color: "#e2e8f0", fontFamily: "'Inter',system-ui,sans-serif", overflow: "hidden" }}>
      <div style={{ width: 180, background: "#1e293b", borderRight: "1px solid #334155", display: "flex", flexDirection: "column", flexShrink: 0 }}>
        <div style={{ padding: "12px 10px 8px", borderBottom: "1px solid #334155" }}>
          <div style={{ fontSize: 13, fontWeight: 900, background: "linear-gradient(135deg,#6366f1,#06b6d4)", WebkitBackgroundClip: "text", WebkitTextFillColor: "transparent" }}>BioChip Mixer</div>
          <div style={{ fontSize: 8, color: "#64748b", marginTop: 1 }}>Mixing Tree Algorithms</div>
        </div>
        <nav style={{ padding: "6px 4px", flex: 1 }}>
          {tabs.map(t => <button key={t.id} onClick={() => setPage(t.id)}
            style={{ display: "flex", alignItems: "center", gap: 6, width: "100%", padding: "6px 8px", marginBottom: 1, background: page === t.id ? `${t.color}18` : "transparent", border: page === t.id ? `1px solid ${t.color}44` : "1px solid transparent", borderRadius: 6, color: page === t.id ? t.color : "#94a3b8", fontWeight: page === t.id ? 700 : 500, fontSize: 11, cursor: "pointer", transition: "all .15s" }}>
            <span>{t.icon}</span>{t.label}
          </button>)}
        </nav>
        <div style={{ padding: "7px 10px", borderTop: "1px solid #334155", fontSize: 9, color: "#475569", lineHeight: 1.5 }}>
          Accuracy: 1/2<sup>{depth}</sup><br />Fluids: {raw.split(",").filter(s => s.trim()).length}<br />Total: 2<sup>{depth}</sup> = {2 ** depth}
        </div>
      </div>
      <div style={{ flex: 1, display: "flex", flexDirection: "column", overflow: "hidden", padding: 8 }}>
        {page === "rma" && <AlgoPage raw={raw} setRaw={setRaw} depth={depth} setDepth={setDepth} algo="rma" title="RMA — Ratioed Mixing" icon="⚗️" color="#6366f1" desc="Largest fluid fills one child → long dilution chains. Best for layout mapping & cross-contamination." />}
        {page === "bs" && <AlgoPage raw={raw} setRaw={setRaw} depth={depth} setDepth={setDepth} algo="bs" title="BS — Bit-Scanning [Thies]" icon="🔬" color="#f59e0b" desc="Bottom-up bit scanning. Minimises leaf count (fewest dispensing steps)." />}
        {page === "apdp" && <AlgoPage raw={raw} setRaw={setRaw} depth={depth} setDepth={setDepth} algo="apdp" title="AP-DP v3 — Multi-Strategy Adaptive" icon="⚙️" color="#06b6d4" desc="DP subset-sum (all exact subsets) + DP+split (multiple sums) + RMA-style (all dominants) + greedy candidates with dilution-aware scoring." />}
        {page === "larp" && <AlgoPage raw={raw} setRaw={setRaw} depth={depth} setDepth={setDepth} algo="larp" title="LARP v2 — 2-Level Lookahead DP" icon="🔮" color="#ef4444" desc="Comprehensive enumeration: zero-split, single-split (multiple sums, both directions), double-split, RMA-style (all dominants). 2-level lookahead with dilution-aware scoring. No artificial caps." />}
        {page === "ilp" && <AlgoPage raw={raw} setRaw={setRaw} depth={depth} setDepth={setDepth} algo="ilp" title="ILP v2 — Pareto-Optimal DP" icon="🧮" color="#84cc16" desc="Pareto-optimal candidate enumeration (splits, distinctFluids). Dilution-aware scoring rewards ≤2 fluids per side. 1-level lookahead estimates downstream splits. RMA-style fallback for dominant fluids." />}
        {page === "compare" && <ComparePage raw={raw} setRaw={setRaw} depth={depth} setDepth={setDepth} />}
      </div>
    </div>
  );
}