# LARP — Lookahead-Augmented Recursive Partitioning

A novel algorithm for constructing mixing trees on digital microfluidic (DMF) biochips.

## Motivation

Two established algorithms for automatic solution preparation exist:

| Algorithm | Strategy | Strength | Weakness |
|-----------|----------|----------|----------|
| **RMA** (Roy et al.) | Greedy: always assigns dominant fluid first | High dilution subtree length (*l*) → better layout mapping | Higher mix count (*m*); greedy may miss globally better splits |
| **BS** (Thies et al.) | Bottom-up bit-scanning | Minimum leaf count → fewest dispensing steps | Poor fluid separation; no optimization for cross-contamination |

**Key insight**: Both algorithms are *myopic* — they make locally optimal decisions without considering how the current partition affects future levels of the tree. LARP addresses this by introducing **lookahead scoring**.

## Algorithm Overview

```
LARP(P, L):
  Input:  Partition P = {x₁:a₁, ..., xₙ:aₙ}, Volume L
  Output: Mixing tree root node

  1. Base case: if |P| = 1 or L ≤ 1, return leaf node

  2. ENUMERATE all valid partitions (P₁, P₂) where sum(P₁) = sum(P₂) = L/2:
     a. Use DP subset-sum to find all subsets of whole fluids summing to L/2
        → These are "zero-split" candidates (no fluid is divided between sides)
     b. If no exact subset: find closest sum S < L/2, then for each fluid
        in the remainder with coefficient ≥ (L/2 − S), create a "one-split"
        candidate by splitting exactly (L/2 − S) units from that fluid
     c. Fallback: RMA-style dominant-fluid partition

  3. SCORE each candidate with lookahead:
     For candidate (P₁, P₂):
       immediateCost  = splits × 20 + minDistinct × 3
       lookaheadCost  = quickScore(P₁, L/2) + quickScore(P₂, L/2)
       totalCost      = immediateCost + 8 × lookaheadCost

  4. SELECT candidate with minimum totalCost

  5. RECURSE: left ← LARP(P₁, L/2),  right ← LARP(P₂, L/2)
```

## Key Components

### 1. Exhaustive Partition Enumeration (`larpEnumPartitions`)

Unlike AP-DP which generates only 4 fixed strategies, LARP uses DP to find **all** achievable subset sums:

```
dp[s] = list of bitmasks achieving sum s using whole-fluid subsets
```

Each bitmask represents which fluids go entirely to P₁. This guarantees finding zero-split partitions whenever they exist — something greedy approaches may miss.

### 2. Lookahead Scoring (`larpScore` + `larpQuickScore`)

The scoring function has two components:

**Immediate Cost** (current level quality):
- `splits × 20`: Heavily penalizes fluid splits (splitting a fluid between both children increases cross-contamination)
- `minDistinct × 3`: Prefers partitions where at least one child has few distinct fluids (promotes dilution subtrees)

**Lookahead Cost** (next level estimate):
- For each child P₁ and P₂, runs a quick DP check: "can the *next* level achieve a zero-split partition?"
- If yes → cost 0; if no → cost 1
- Weighted by factor 8 in the total score

This means LARP **prefers partitions that lead to better partitions at the next level**, even if the immediate quality is similar.

### 3. Quick Score Estimator (`larpQuickScore`)

A lightweight DP that only checks *reachability* of sum = L/4 using whole fluids:
```
dp[0] = 1
for each fluid f with coefficient v:
  for s = target down to v:
    dp[s] |= dp[s - v]
return dp[target] ? 0 : 1
```
Time: O(n × L/4), runs once per candidate child — negligible overhead.

## Why LARP Is Different

| Feature | RMA | BS | AP-DP | **LARP** |
|---------|-----|-----|-------|----------|
| Partition generation | Greedy (1 candidate) | Fixed (bit-based) | 4 fixed strategies | **All DP subsets** |
| Scoring considers future levels? | No | No | No | **Yes (1-level)** |
| Guarantees zero-split if possible? | No | No | Sometimes | **Always** |
| Complexity per level | O(n) | O(n) | O(n·L) | O(n·L + C·n·L/4) |

Where C = number of candidates (capped at ~8 per sum), n = fluids, L = volume.

## Expected Performance

| Metric | vs RMA | vs BS | Reasoning |
|--------|--------|-------|-----------|
| m (mixes) | ≤ or better | Similar | Zero-split partitions reduce unnecessary fluid splitting |
| l (dilution) | Similar | Higher | Prefers pure single-fluid children |
| splits (fluid splits) | **Lower** | **Lower** | Primary optimization objective |
| p (parallelism) | Similar | Similar | Structural property depends on input |

## Complexity Analysis

- **Time per partition step**: O(n · L/2) for DP + O(C · n · L/4) for lookahead
- **Total**: O(d · n · L) where d = tree depth
- **Space**: O(L/2) for DP arrays
- **Practical**: Runs in <50ms for typical inputs (n ≤ 12, d ≤ 14)

## Example

For ratio **2:3:5:7:11:13:87** (d=7, L=128):

At the root, L=128, half=64. LARP enumerates all subsets summing to 64:
- Subset {x7:87} → sum=87 > 64, skip
- Subset {x1:2, x5:11, x6:13} → sum=26, too low
- Subset {x7:87} can't be split, so LARP tries single-split candidates
- When combined with splitting x7, LARP finds `P1={x7:64}` and `P2={x7:23, x1:2, ...}`
- The **lookahead** then checks: "can P₂ be split without splitting any fluid at the next level?"
- This guides LARP to choose partitions that cascade cleanly down the tree
