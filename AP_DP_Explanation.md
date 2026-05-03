# AP-DP v3: Adaptive Partitioning via Dynamic Programming
## A First-Principles Explanation with Dry Run & Algorithm Comparison

---

## 1. First Principles — What Are We Solving?

### The Physical Context: Microfluidic Biochips

A Digital Microfluidic Biochip (DMF biochip) manipulates tiny droplets of fluid on a chip. The key operation is **mixing**: two equal-volume droplets are merged. Every mix on the chip is **exactly 50:50** — you cannot mix in any other ratio directly.

```
  Droplet A  +  Droplet B  →  Mix (50% A, 50% B)
   (volume V)    (volume V)     (volume 2V)
```

### The Problem

Given a desired fluid **ratio** like `A:B:C = 1:1:1` (equal parts of three reagents), design the smallest binary tree of mixing operations whose **leaf outputs** produce the correct proportions.  

**Key constraints:**
- Every internal node mixes its two children's output **50:50**
- A tree of depth `d` produces a total of `2^d` units of output
- Leaf nodes are "dispensing" operations — they output one unit of a pure fluid

**Goal:** Build this tree **minimising cost metrics** like:
- `m` — number of mix operations (internal nodes)
- `splits` — number of fluids that appear on **both** sides of any single split (a "fluid split" wastes reagent — the same fluid must be dispensed multiple times)
- `l` — length of **dilution subtrees** (subtrees with ≤ 2 distinct fluids produce via serial dilution, which is cheap)

---

## 2. Pre-processing: `ratioApprox(ratios, d)`

Before building the tree, the continuous ratios are converted to integers that sum to `L = 2^d`.

### Formula
```
approx[i] = round(ratio[i] / sum(ratios) × 2^d)
```
The last element is adjusted to ensure `sum(approx) == 2^d` exactly.

### Example

```
Input ratios:  [1, 1, 1]     d = 3
L = 2^3 = 8
Each proportion = 1/3

approx[0] = round(1/3 × 8) = round(2.667) = 3
approx[1] = round(1/3 × 8) = round(2.667) = 3
approx[2] = 8 - (3 + 3) = 2     ← adjusted last element

Result: [3, 3, 2]  (sum = 8 ✓)
Fluids: A=3, B=3, C=2
```

> **Important:** The approximation introduces a small error. With `d=3`, accuracy is `1/8 = 12.5%` per slot. Larger `d` gives better precision at the cost of more operations.

---

## 3. The Core Idea: Recursive Binary Halving

The tree is built **top-down**. At each node:
1. We have a polynomial `P` — a map from fluid name → integer volume (e.g., `{A:3, B:3, C:2}`)
2. The total volume at this node is `L`
3. We must split `P` into **two sub-polynomials** `P1` and `P2`, each summing to exactly `L/2`
4. We recurse: left child gets `P1`, right child gets `P2`

```
Node with {A:3, B:3, C:2}, L=8
          ↓  split so each child gets L/2 = 4
   Left: {A:3, B:1}    Right: {B:2, C:2}
   (sum=4)               (sum=4)
```

**Why must each child sum to L/2?**
Because the chip mixes children 50:50. If left outputs 4 units and right outputs 4 units, the parent produces 8 units in the ratio of their combined composition.

**Base cases:**
- If `P` has only **1 fluid** → this is a **leaf** (pure dispensing step)
- If `L ≤ 1` → forced leaf

---

## 4. AP-DP v3: The Algorithm in Detail

AP-DP (Adaptive Partitioning via Dynamic Programming) improves on greedy algorithms by **generating multiple candidate partitions** through 5 different strategies, then **scoring and selecting** the best one.

### Phase 1: Generate Candidate Partitions

#### ▸ Strategy A — Exact Whole-Fluid Subset (Zero-Split)

**Question:** Can we put some fluids entirely in P1 and the rest entirely in P2, such that P1 sums to exactly `half`?

This is the **Subset Sum Problem**. AP-DP enumerates ALL subsets (up to 15) that sum exactly to `half`.

```
For P={A:3, B:3, C:2}, half=4:
  Try {A}=3≠4, {B}=3≠4, {C}=2≠4
  Try {A,B}=6≠4, {A,C}=5≠4, {B,C}=5≠4
  → NO exact subset found
```

No zero-split partition exists here (3+3+2=8, no subset sums to 4).

#### ▸ Strategy B — DP Closest Sum + Single Fluid Split

Since no exact subset exists, find the **closest achievable whole-fluid sum** to `half`, then **split exactly one fluid** to close the gap.

```
Achievable whole-fluid sums ≤ 4:
  {}  → 0
  {C} → 2
  {A} → 3
  {B} → 3
  (6, 5, 5 all exceed 4 → excluded)

Closest to 4: sum=3 (via {A} or {B})
  Need = half - sum = 4 - 3 = 1

Try splitting each fluid NOT in the subset:
  Subset {A=3}: remaining fluids are B=3 and C=2
    - Take 1 from B: P1={A:3, B:1}, P2={B:2, C:2}  [B is split ⚠]
    - Take 1 from C: P1={A:3, C:1}, P2={B:3, C:1}  [C is split ⚠]
  Subset {B=3}: symmetric → {B:3,A:1}|{A:2,C:2} and {B:3,C:1}|{A:3,C:1}
  
  sum=2 via {C=2}:
    Need = 4 - 2 = 2
    - Take 2 from A: P1={C:2, A:2}, P2={A:1, B:3}  [A split ⚠]
    - Take 2 from B: P1={C:2, B:2}, P2={A:3, B:1}  [B split ⚠]

Candidates generated: 6 total (all with exactly 1 split)
```

#### ▸ Strategy C — RMA-Style (For Dominant Fluids ≥ half)

If any fluid has volume ≥ `half`, fill one child entirely with that fluid (truncated to `half`).

```
For {A:3, B:3, C:2}, half=4:
  A=3 < 4, B=3 < 4, C=2 < 4
  → No dominant fluid. Strategy C adds nothing.
```

If A were 5 (≥4), Strategy C would give: `P1={A:4}, P2={A:1, B:3, C:2}` (A split).

#### ▸ Strategy D — Greedy Balanced (Large Fluids First)

Fill P1 greedily from the largest fluid down, stopping when P1 reaches `half`.

```
Sorted: A=3, B=3, C=2.  half=4
  A: s1(0)≤s2(0), s1+3=3≤4 → P1={A:3}, s1=3
  B: s2(0)+3=3≤4 → P2={B:3}, s2=3
  C: s1(3)<4, need=1, C=2≥1 → P1={A:3,C:1}=4, P2={B:3,C:1}=4  [C split ⚠]
Candidate D: P1={A:3,C:1} | P2={B:3,C:1}
```

#### ▸ Strategy E — Reverse Greedy (Fill P2 First)

Same as D, but P2 is filled first.

```
Sorted: A=3, B=3, C=2.  half=4
  A: s2(0)≤s1(0), s2+3=3≤4 → P2={A:3}, s2=3
  B: s1(0)+3=3≤4 → P1={B:3}, s1=3
  C: s2(3)<4, need=1, C=2≥1 → P2={A:3,C:1}=4, P1={B:3,C:1}=4  [C split ⚠]
Candidate E: P1={B:3,C:1} | P2={A:3,C:1}  (same as D, different labeling)
```

---

### Phase 2: Filter, Deduplicate, Score

#### Validation
Only keep candidates where BOTH P1 and P2 sum exactly to `half`.

#### Deduplication
Remove mirror images (P1↔P2 is same partition).

#### Scoring — `scorePartitionV3(P1, P2, L)`

This is the **key differentiator** of AP-DP. The score function is:

```
score = splits × 1000
       - dilutionBonus × 50
       + totalDistinct × 10
       + minFluids × 5
```

Where:
- `splits` = number of fluids appearing in **both** P1 and P2 (heavily penalised!)
- `dilutionBonus` = 1 for each side where the **dominant fluid > 50%** of that side's volume (rewards dilution potential)
- `totalDistinct` = |fluids in P1| + |fluids in P2| (fewer is better)
- `minFluids` = min(|P1 fluids|, |P2 fluids|) (prefer asymmetric splits)

**Lower score = better partition.**

---

### Phase 3: Dry Run — Scoring All Candidates

Let's score our candidates for `{A:3, B:3, C:2}, L=8, half=4`:

| Candidate | P1 | P2 | splits | dom1 | dom2 | dilBon | totDist | minF | **Score** |
|---|---|---|---|---|---|---|---|---|---|
| B: `{A:3,B:1}` \| `{B:2,C:2}` | B in both | B in both | **1** | 3/4=0.75✓ | 2/4=0.5✗ | 1 | 4 | 2 | 1000−50+40+10=**1000** |
| B: `{A:3,C:1}` \| `{B:3,C:1}` | C in both | C in both | **1** | 3/4=0.75✓ | 3/4=0.75✓ | 2 | 4 | 2 | 1000−100+40+10=**950** |
| B: `{C:2,A:2}` \| `{A:1,B:3}` | A in both | A in both | **1** | 2/4=0.5✗ | 3/4=0.75✓ | 1 | 4 | 2 | 1000−50+40+10=**1000** |
| B: `{C:2,B:2}` \| `{A:3,B:1}` | B in both | B in both | **1** | 2/4=0.5✗ | 3/4=0.75✓ | 1 | 4 | 2 | 1000−50+40+10=**1000** |
| D: `{A:3,C:1}` \| `{B:3,C:1}` | (duplicate) | — | — | — | — | — | — | — | dup of row 2 |

**Winner: `{A:3, C:1}` | `{B:3, C:1}` with score 950**

This wins because **both sides have a dominant fluid** (A=3/4=75% on left, B=3/4=75% on right), giving the best dilution potential.

---

### Phase 4: Recurse

```
Root: {A:3, B:3, C:2}, L=8
  Best split: P1={A:3, C:1} | P2={B:3, C:1}
```

**Left subtree:** `{A:3, C:1}, L=4, half=2`

```
Strategy A: {A}=3≠2, {C}=1≠2, {A,C}=4≠2 → No exact subset
Strategy B: Achievable ≤2: {0, 1(C)}
  sum=1(C), need=1:
    Take 1 from A: P1={C:1,A:1}=2, P2={A:2}=2  [A split ⚠]
  sum=0, need=2:
    Take 2 from A: P1={A:2}=2, P2={A:1,C:1}=2  [A split ⚠]
Strategy C: A=3≥2 → P1={A:2}, P2={A:1,C:1}           [A split ⚠]

Score P1={C:1,A:1}|P2={A:2}:
  splits=1, dom1=1/2=0.5✗, dom2=2/2=1.0✓ → dilBon=1, totDist=3, minF=1
  score = 1000−50+30+5 = 985

Score P1={A:2}|P2={A:1,C:1}:
  splits=1, dom1=2/2=1.0✓, dom2=1/2=0.5✗ → dilBon=1, totDist=3, minF=1
  score = 1000−50+30+5 = 985  (same)

Pick first: P1={A:2} | P2={A:1,C:1}
  Left child: {A:2} → single fluid → LEAF (A)
  Right child: {A:1,C:1}, L=2, half=1
    Strategy A: {A}=1=half ✓ → P1={A:1}, P2={C:1}  [0 splits!]
    Score: splits=0 → wins easily
    → LEAF(A) and LEAF(C)
```

**Left subtree result:**
```
    Mix(A:3,C:1)
     /          \
  Leaf(A)     Mix(A:1,C:1)
              /           \
           Leaf(A)       Leaf(C)
```

**Right subtree:** `{B:3, C:1}, L=4` — symmetric to left, replace A→B

```
    Mix(B:3,C:1)
     /          \
  Leaf(B)     Mix(B:1,C:1)
              /           \
           Leaf(B)       Leaf(C)
```

---

### Final AP-DP Tree

```
                   Mix(A:3, B:3, C:2)       L=8
                  /                    \
       Mix(A:3, C:1)  L=4        Mix(B:3, C:1)  L=4
        /          \               /           \
    Leaf(A)   Mix(A:1,C:1)  Leaf(B)     Mix(B:1,C:1)
                /        \                /         \
            Leaf(A)    Leaf(C)       Leaf(B)      Leaf(C)
```

**Statistics:**
- `m` (mixes) = 5
- `d` (depth) = 3
- `splits` = 1 (C appears in both children of root) + 1 (A appears in both left-left) + 1 (B appears in both right-right) = **3**
- Leaves: A×3, B×3, C×2 → ratio 3:3:2 ✓

> **Note:** The unavoidable consequence is that C must be dispensed in both subtrees. For an exact 1:1:1 ratio at depth 3, some fluid *must* be split — you cannot pack three equal parts into a binary tree of depth 3 without splitting.

---

## 5. Side-by-Side Comparison: Same Example on All Algorithms

### Setup: Ratios [1, 1, 1], d=3 → approximated to A=3, B=3, C=2

---

### 5.1 RMA (Ratioed Mixing Algorithm)

RMA is a **greedy rule-based algorithm** from the paper. It always builds around the **dominant fluid**.

```
rmaPartition({A:3, B:3, C:2}, L=8, half=4):
  Sorted: A=3, B=3, C=2  (tied; A wins by sort order)
  aU = A = 3  < half(4)
    → P1 = {A:3},  P2 = {B:3, C:2},  E = half - aU = 1

  Is E=1 equal to any P2 coefficient? B=3≠1, C=2≠1 → No
  Try subsets of P2 of size ≤ |P2|/2 = 1 that sum to E=1? B=3≠1, C=2≠1 → No
  Find next largest in P2 where value > E=1: B=3>1 → split B!
    P1 += {B:1} → P1={A:3,B:1}=4✓
    P2[B] = 3-1 = 2 → P2={B:2,C:2}=4✓

Split: P1={A:3, B:1} | P2={B:2, C:2}
```

RMA is **deterministic** — it always starts with the largest fluid and follows a fixed procedure. It does NOT explore other candidates.

**Recursing RMA:**
- `{A:3, B:1}, L=4, half=2`:
  - Sorted: A=3. aU=3≥2 → P1={A:2}, P2={A:1,B:1} [A split]
  - Recurse: `{A:2}→Leaf(A)`, `{A:1,B:1},L=2→P1={A:1}|P2={B:1}→Leaf(A),Leaf(B)`
- `{B:2, C:2}, L=4, half=2`:
  - aU=B=2=half → P1={B:2}, P2={C:2} [0 splits]
  - Both become Leaf(B), Leaf(C)

```
RMA Tree:
                  Mix(A:3, B:3, C:2)
                 /                   \
        Mix(A:3, B:1)             Mix(B:2, C:2)
         /          \              /            \
      Leaf(A)   Mix(A:1,B:1)   Leaf(B)        Leaf(C)
                /          \
            Leaf(A)       Leaf(B)
```

**RMA Statistics:**
- `m` = 3, `d` = 3, `splits at root` = 1 (B split), A split below, total = **2 splits**

---

### 5.2 BS (Bit-Scanning, Thies et al.)

BS is **bottom-up** and purely mechanical. It reads the **binary representation** of each fluid's count and builds the tree bit by bit.

```
approx = [A:3, B:3, C:2]  d=3

Build bit table:
  Bit 0 (value 1): A= 3 → 011 → bit0=1 ✓,  B=3 → 011 → bit0=1 ✓,  C=2 → 010 → bit0=0
  Bit 1 (value 2): A= 3 → 011 → bit1=1 ✓,  B=3 → 011 → bit1=1 ✓,  C=2 → 010 → bit1=1 ✓
  Bit 2 (value 4): A= 3 → 011 → bit2=0,     B=3 → 011 → bit2=0,     C=2 → 010 → bit2=0

lvIn[0] = [A, B]   ← fluids with bit 0 set
lvIn[1] = [A, B, C] ← fluids with bit 1 set
lvIn[2] = []

Process bottom-up:
  b=0: prev=[], items=[Leaf(A), Leaf(B)]
         → Pair them: Mix(A,B). prev = [Mix(A,B)]

  b=1: prev=[Mix(A,B)], new leaves=[Leaf(A), Leaf(B), Leaf(C)]
       items = [Mix(A,B), Leaf(A), Leaf(B), Leaf(C)]
         → Pair 0,1: Mix(Mix(A,B), Leaf(A))
         → Pair 2,3: Mix(Leaf(B), Leaf(C))
       prev = [Mix(Mix(A,B),A),  Mix(B,C)]

  b=2: prev=[Mix(Mix(A,B),A), Mix(B,C)], no new leaves
       items = [Mix(Mix(A,B),A), Mix(B,C)]
         → Pair 0,1: Mix(Mix(Mix(A,B),A), Mix(B,C))
       prev = [root]
```

```
BS Tree:
           Mix(A:3,B:3,C:2)
          /                 \
 Mix(Mix(A,B), A)         Mix(B,C)
  /              \         /      \
Mix(A,B)       Leaf(A)  Leaf(B) Leaf(C)
/        \
Leaf(A)  Leaf(B)
```

**BS Statistics:**
- `m` = 4, `d` = 4, `splits` = A appears in two branches + B appears in two branches = **2 splits**
- BS has **deeper** trees due to bottom-up construction

> **Note:** BS minimises leaf count in principle but here has more mixes and deeper depth than RMA/AP-DP.

---

### 5.3 AP-DP v3 ← (explained above)

Best split chosen was `{A:3,C:1}|{B:3,C:1}` due to dilution scoring.

**AP-DP Statistics:**
- `m` = 5, `d` = 3, `splits` = **3** (worse! — C is split at root, A is split in left, B in right)

> 🔑 **Key insight:** AP-DP's greedy-local scoring sometimes chooses a partition that looks good locally (both sides have a dominant fluid) but creates **more total splits globally**. This is why LARP uses 2-level lookahead.

---

### 5.4 LARP v2 (Lookahead-Augmented Recursive Partitioning)

LARP generates the same candidates as AP-DP but scores them with **2-level lookahead** — it estimates how many splits each child subtree will incur downstream.

```
For {A:3, C:1} | {B:3, C:1}:
  Look ahead into {A:3,C:1}:
    Quick DP: achievable sums ≤ 2 from {A:3,C:1} → {0,1(C),3>2}
    Cannot reach 2 exactly → minSplits=1 downstream

For {A:3, B:1} | {B:2, C:2}:
  Look ahead into {A:3,B:1}:
    Achievable ≤ 2 → {0,1(B),3>2} → cannot reach 2 → minSplits=1
  Look ahead into {B:2,C:2}:
    Achievable ≤ 2 → {0,2(B),2(C)} → exact! → minSplits=0 ✓

LARP score for {A:3,C:1}|{B:3,C:1}: +2 children splits downstream (1+1=2)
LARP score for {A:3,B:1}|{B:2,C:2}: +1 child split downstream (1+0=1) ← WINS
```

LARP would pick `{A:3,B:1}|{B:2,C:2}`, matching RMA but reaching it via DP enumeration.

**LARP Tree:**
```
              Mix(A:3,B:3,C:2)
             /                 \
     Mix(A:3, B:1)         Mix(B:2, C:2)
      /           \           /          \
  Leaf(A)    Mix(A:1,B:1)  Leaf(B)     Leaf(C)
             /            \
           Leaf(A)        Leaf(B)
```

**LARP Statistics:** Same as RMA here — `m=3, d=3, splits=2`.

---

### 5.5 ILP v2 (Pareto-Optimal DP)

ILP also uses DP but tracks **Pareto-optimal** allocation states — keeping all allocations that are non-dominated in (splits, fluid count) space.

For this example it also converges to `{A:3,B:1}|{B:2,C:2}` — same as RMA and LARP.

---

### Summary Comparison Table

| Algorithm | Strategy | Candidate Generation | Scoring | `m` | `d` | `splits` |
|---|---|---|---|---|---|---|
| **RMA** | Greedy (dominant fluid first) | 1 candidate (fixed rules) | None | 3 | 3 | 2 |
| **BS** | Bottom-up bit scanning | Deterministic (no choice) | None | 4 | 4 | 2 |
| **AP-DP v3** | Multi-strategy + DP | 5 strategies, up to 15+ candidates | Local dilution score | 5 | 3 | **3** ← worse here |
| **LARP v2** | Multi-strategy + DP | Same as AP-DP, more thorough | **2-level lookahead** | 3 | 3 | 2 |
| **ILP v2** | Pareto-optimal DP | Pareto-optimal allocations | 1-level lookahead | 3 | 3 | 2 |

---

## 6. Where Each Algorithm Deviates — Step-by-Step

```
All algorithms: Start with {A:3, B:3, C:2}, L=8

STEP 1 — Choose the root split:
  RMA:   P1={A:3,B:1} | P2={B:2,C:2}   ← deterministic, follows fixed rule
  BS:    NOT top-down! Builds from bit-level, no explicit split at root
  AP-DP: Generates 6 candidates, picks {A:3,C:1}|{B:3,C:1} ← local score best
  LARP:  Same 6 candidates, picks {A:3,B:1}|{B:2,C:2}      ← lookahead best
  ILP:   Pareto DP → picks {A:3,B:1}|{B:2,C:2}              ← pareto + score

STEP 2 — Left subtree {A:3,B:1} or {A:3,C:1}:
  (RMA/LARP/ILP) {A:3,B:1}, L=4:
    A=3≥half=2 → P1={A:2}|P2={A:1,B:1} → A split, then two leaves
  (AP-DP) {A:3,C:1}, L=4:
    No exact subset, split A: P1={A:2}|P2={A:1,C:1} → A split, then two leaves

STEP 3 — Right subtree {B:2,C:2} or {B:3,C:1}:
  (RMA/LARP/ILP) {B:2,C:2}, L=4:
    B=2=half → P1={B:2}|P2={C:2} → PERFECT SPLIT, no fluid split ✓
  (AP-DP) {B:3,C:1}, L=4:
    Must split B just like left subtree → additional split ✗

KEY DIVERGENCE: AP-DP's root-level choice creates a worse right subtree.
```

---

## 7. Visualising the Difference

```
RMA / LARP / ILP:                    AP-DP:
                                        
  Mix(3:3:2)                          Mix(3:3:2)
  /         \                         /           \
Mix(3:1)   Mix(2:2)               Mix(3:1)     Mix(3:1)
 / \        /    \                 / \           / \
A  Mix(1:1) B    C               A  Mix(1:1)   B  Mix(1:1)
    / \                               / \           / \
   A   B                             A   C         B   C

m=3, splits=2                      m=5, splits=3 ← additional split!
```

---

## 8. Understanding the Split Metric

A **"fluid split"** means the same fluid is dispensed from two different subtrees. This is physically costly because:
1. The fluid must be loaded at two locations on the chip
2. More droplets of that fluid must be dispensed in total
3. Contamination risk increases

For the ratio `3:3:2` at `d=3`, the **minimum possible splits** is **2** (A must be split because in the RMA/LARP tree it appears in two leaves under different sub-mixes). No algorithm can do better.

AP-DP scores **3 splits** because its root choice forces **C** to also be split at the root, unnecessarily.

---

## 9. Suggested Improvements to AP-DP

### 🔧 Improvement 1: Add 2-Level Lookahead to Scoring

**Problem:** The current `scorePartitionV3` is purely local — it doesn't know that choosing `{A:3,C:1}|{B:3,C:1}` will force C to split downstream in BOTH children.

**Fix:** Borrow LARP's `larpDeepScore` idea — after generating candidates, estimate the minimum additional splits each child will need.

```js
// Pseudo-code: enhanced score
function scoreWithLookahead(P1, P2, L) {
  const baseScore = scorePartitionV3(P1, P2, L);
  const la1 = estimateDownstreamSplits(P1, L/2);  // 1-level lookahead
  const la2 = estimateDownstreamSplits(P2, L/2);
  return baseScore + (la1 + la2) * 500;            // penalise downstream splits
}
```

### 🔧 Improvement 2: Strategy B Should Include Subset-then-Split More Thoroughly

Currently, Strategy B only tries the **top 5 closest sums** × each candidate fluid to split. For problems with many fluids, this misses valid combinations.

**Fix:** In Strategy B, try ALL achievable sums (not just top 5) when the fluid count is manageable.

```js
// More thorough: try ALL achievable sums, not just top 5
const closestSums = findClosestSums(items, half, items.length * 2);
```

### 🔧 Improvement 3: Track Global Splits in Scoring

Rather than only counting whether a fluid is "split at this node", accumulate a **global split counter** as the tree is built and use it to tie-break.

```js
function buildAPDP(P, L, maxD, lv = 0, globalSplits = {}) {
  // Pass globalSplits object down and up the recursion
  // Penalise candidates that re-split an already-split fluid
}
```

### 🔧 Improvement 4: Memoise Subproblems

AP-DP rebuilds identical sub-polynomials multiple times. A memoisation cache on `(P, L)` would avoid redundant computation.

```js
const memo = new Map();
function buildAPDPMemo(P, L, ...) {
  const key = JSON.stringify([Object.entries(P).sort(), L]);
  if (memo.has(key)) return memo.get(key);
  const result = buildAPDPCore(P, L, ...);
  memo.set(key, result);
  return result;
}
```

### 🔧 Improvement 5: Parallelism-Aware Scoring

For chip layouts where `p` (max nodes at any level) matters, add a balance penalty that rewards symmetric trees.

```js
const balancePenalty = Math.abs(
  Object.keys(P1).length - Object.keys(P2).length
) * 3;
```

---

## 10. Quick Reference: Algorithm Philosophies

| | RMA | BS | AP-DP |
|---|---|---|---|
| **Direction** | Top-down | Bottom-up | Top-down |
| **Choices explored** | 1 (deterministic) | 0 (fixed) | Many (up to 30+) |
| **Lookahead** | None | None | None (local only) |
| **Strength** | Dilution chains (l metric) | Few leaves | Dilution score, flexible |
| **Weakness** | Greedy, can miss global optima | Deep trees, no split awareness | Local scoring can misfire |
| **Complexity** | O(n log n) | O(d·n) | O(2^n × d) |

---

## 11. Key Takeaways

1. **Every split in the tree is a 50:50 mix** — there is no other operation.
2. **ratioApprox snaps ratios to powers of 2** — this introduces approximation error that shrinks as `d` grows.
3. **AP-DP generates candidates via 5 strategies** (exact subset, DP+split, RMA-style, greedy, reverse-greedy) and scores them locally.
4. **The local scoring heuristic is AP-DP's key weakness** — it can choose root-level partitions that look good locally but cause more splits downstream.
5. **LARP and ILP address this** with lookahead mechanisms that estimate downstream costs.
6. **For 3:3:2, the minimum split count is 2** — RMA, LARP, and ILP all achieve this; AP-DP gets 3 due to its local scoring.
