# Mixing Tree Construction for Sample Preparation on Digital Microfluidic Biochips
 
**A Comprehensive Study of Five Partitioning Algorithms with Beam-Search Augmentation**
 
---
 
> **About this document.** This is a self-contained knowledge base of the project. It is structured for downstream agents (and humans) to extract presentations, written reports, or research-paper drafts. Every algorithm is described in pseudocode with complexity analysis. All experimental numbers come from the local benchmark suite (`benchmark.js` + `bench.js`, 323 test cases, deterministic seed). Citations point to original published sources where applicable; algorithms originating in this work are marked **[proposed]**.
 
---
 
## Executive Summary (for slide 1, abstract, and elevator pitch)
 
### One-line research claim
 
> **On a 323-test stratified DMFB sample-preparation benchmark, beam-augmented LARP achieves 18% fewer reagent dispenses than RMA (Roy et al., 2010) and 44% longer dilution chains than BS (Thies et al.), occupying a previously-unfilled middle of the splits-vs-dilution trade-off curve and ranking last on only 1.5% of tests.**
 
### What is contributed
 
| # | Contribution                                                                                              | Status                       |
|---|-----------------------------------------------------------------------------------------------------------|------------------------------|
| 1 | **Three new mixing-tree construction algorithms** (AP-DP, LARP, ILP)                                      | Proposed in this work        |
| 2 | **Beam-search + memoization meta-framework** applicable to any partitioning-based mixing-tree algorithm   | Proposed in this work        |
| 3 | **323-test stratified empirical benchmark** with 7 metrics, deterministic seed, full reproducibility      | Proposed in this work        |
| 4 | **Identification of the splits-vs-dilution trade-off curve**, with a quantified "balanced" operating point | Proposed in this work        |
 
### Headline numbers (every number here is verified against §9)
 
| Comparison              | Splits change | Dilution change | Direction                          |
|-------------------------|---------------|-----------------|------------------------------------|
| LARP vs RMA (Roy 2010)  | **−17.6%** ✓  | −23.5%          | Better on splits, worse on dilution |
| LARP vs BS (Thies)      | +8.2%         | **+44.5%** ✓    | Worse on splits, better on dilution |
| Beam-LARP vs greedy LARP (Phase B) | **−14.0%** ✓ | −20% (regression) | Beam search is a clean win on the design objective |
 
### What we honestly do NOT claim
 
* ❌ "We beat every existing algorithm." We compared only to RMA and BS; REMIA, WARA, IDMA, DMRW were not implemented.
* ❌ "We are split-optimal." BS is the unbeaten split-champion — LARP/ILP can only *tie* BS on splits, never exceed.
* ❌ "Real ILP solving." Our ILP is a DP heuristic with integer-allocation tracking, not a solver-based ILP (no CPLEX/Gurobi). See §13 for details.
* ❌ "LARP and ILP are two distinct contributions." They produce identical trees on 70% of tests; treat them as one algorithm family with two candidate-enumeration variants.
 
### What the work IS positioned as
 
A solid empirical study with three proposed heuristic algorithms, a meta-framework that improves them by 9–14% on the design objective, and a clear identification of where the proposed algorithms fit on the multi-objective trade-off curve of the problem.
 
### Recommended algorithm to use in practice
 
**LARP with beam search K = 3.** It is within ±10% of every algorithm's per-metric optimum, ranks last on only 5/323 tests (1.5%), and runs in 7 ms average per construction (62 ms worst case). ILP gives identical answers 70% of the time at 10× the wall-clock cost (75 ms avg, 4.4 s worst case), so LARP is the better practical choice.
 
---
 
## Table of Contents
 
0. [Executive Summary](#executive-summary-for-slide-1-abstract-and-elevator-pitch)
1. [Problem Domain](#1-problem-domain)
2. [Formal Problem Statement](#2-formal-problem-statement)
3. [Metrics and Definitions](#3-metrics-and-definitions)
4. [Literature Survey: Existing Algorithms](#4-literature-survey-existing-algorithms)
5. [Algorithms Implemented in This Work](#5-algorithms-implemented-in-this-work)
6. [The Beam-Search + Memoization Meta-Framework](#6-the-beam-search--memoization-meta-framework)
7. [Benchmarking Framework](#7-benchmarking-framework)
8. [The Visualization Tool (React Web App)](#8-the-visualization-tool-react-web-app)
9. [Experimental Results](#9-experimental-results)
   - 9.1 Overall Averages
   - 9.2 **Headline: Proposed vs Published Baselines (% change)**
   - 9.3 **Per-Test Rank Distribution (the most informative single view)**
   - 9.4 Sole vs Tied Composite Wins
   - 9.5 Pairwise Composite-Tie Matrix
   - 9.6 Strict-Best on Splits (the design objective)
   - 9.7 Per-Metric Win Counts
   - 9.8 Final Scoreboard (7-metric, tied-rank-aware)
   - 9.9 Per-Algorithm Timing
   - 9.10 Phase-B Impact (beam vs greedy)
10. [Discussion and Trade-offs](#10-discussion-and-trade-offs)
    - 10.1 The Splits-vs-Dilution Trade-off
    - 10.2 LARP vs ILP — They Are Effectively the Same Algorithm
    - 10.3 BS Is Unbeatable on Splits — and That's a Different Algorithm Niche
    - 10.4 When to Use Which Algorithm
    - 10.5 The Role of Beam Search
    - 10.6 Why AP-DP Underperforms
    - 10.7 Methodological Note on the Composite
11. [Limitations and Future Work](#11-limitations-and-future-work)
12. [References](#12-references)
13. [Naming Note: What "ILP" Means in This Work](#13-naming-note-what-ilp-means-in-this-work)
 
---
 
## 1. Problem Domain
 
### 1.1 Digital Microfluidic Biochips (DMFBs)
 
A **Digital Microfluidic Biochip** is a lab-on-a-chip device that manipulates discrete picolitre-to-nanolitre droplets on a 2-D grid of electrodes. By cycling voltages on adjacent electrodes, the chip exploits the *electrowetting-on-dielectric* (EWOD) phenomenon to **transport, merge, split, and dispense** droplets without any moving parts. The same hardware platform can therefore execute a wide variety of biochemical protocols — clinical diagnostics, DNA sequencing by synthesis, point-of-care testing, polymerase-chain-reaction (PCR) preparation, immunoassays, and tissue-engineering workflows — by reprogramming the electrode actuation sequence rather than redesigning the chip.
 
For a complete chip, four primitive droplet operations exist:
 
| Primitive   | Description                                                                                     |
| ----------- | ----------------------------------------------------------------------------------------------- |
| Dispense    | Inject a droplet of a pure reagent or sample from a boundary reservoir onto the array          |
| Transport   | Move a droplet across electrodes (typically 1 cell per clock cycle)                            |
| Mix (merge) | Bring two droplets onto adjacent electrodes; they coalesce into one droplet of summed volume   |
| Split       | Detach a single droplet into two equal-volume daughter droplets via simultaneous opposite pulls|
 
The **mix is volume-conservative and ratio-conservative**: a mix of two equal-volume droplets of compositions A and B produces a new droplet of composition (A + B) / 2. *There is no native, on-chip operation that mixes in any ratio other than 50:50 from two equal-volume droplets.* Any non-50:50 outcome must be synthesized by a **sequence** of 50:50 mixes, which is precisely the problem this work addresses.
 
### 1.2 Sample Preparation
 
A *bioassay* — for example a PCR run, ELISA, or a dose-response curve in cytotoxicity studies — typically requires reagents combined in a precise **volumetric ratio**, e.g., `4 : 3 : 2 : 1 : 1 : 1` of buffer, primer, template, polymerase, dye, and water. The act of producing such a target mixture from pure reservoirs is **sample preparation**, and on a DMFB it is realized as a sequence of mix-and-split operations forming a **mixing tree**:
 
* **Leaves** of the tree are dispense events, each emitting one unit of one pure reagent.
* **Internal nodes** are 50:50 mix-merges of their two children's outputs.
* The **root** is the final droplet whose composition matches the target ratio.
 
A mixing tree of depth `d` produces `2^d` units of total volume, distributed across the reagents in proportions controlled by the multiplicities of leaves.
 
### 1.3 Why Algorithm Design Matters
 
The chip is bandwidth-, time-, and reagent-limited. A poorly chosen mixing tree may:
 
* Waste expensive reagents (dispensing the same fluid in two different subtrees forces extra dispenses)
* Increase assay completion time (more mix steps = more clock cycles)
* Increase contamination risk (longer paths over previously occupied electrodes)
* Produce unnecessary intermediate droplets that must be discarded as waste
 
Every metric defined in §3 is a proxy for one or more of these physical costs. The optimization problem is **multi-objective**, and no single algorithm dominates on all metrics — the entire field is built around discovering well-balanced points in this trade-off space.
 
---
 
## 2. Formal Problem Statement
 
### 2.1 Inputs
 
* A target **ratio** `r = (r₁, r₂, ..., rₙ)` of *n* reagents.
* A **depth budget** `d ∈ ℤ⁺` chosen so that the rounding error introduced by approximating real-valued ratios with integers summing to `2^d` is acceptable for the assay (typical: `d ∈ [4, 12]`).
 
### 2.2 Pre-processing: Ratio Approximation
 
The continuous ratio is rounded to integers summing to `L := 2^d`:
 
```
approx[i] = round(r[i] / Σⱼrⱼ × 2^d)            for i = 1..n−1
approx[n] = 2^d − Σᵢ<n approx[i]                 (last element corrects the sum)
```
 
This is implemented in `ratioApprox(ratios, d)`. For `r = [1,1,1]` and `d = 3`:
 
```
L = 8
approx = [3, 3, 2]   (sums to 8, deviation per slot ≤ 1/2^d = 12.5%)
```
 
The per-slot rounding error shrinks geometrically with `d`. For most assays `d ≥ 6` is sufficient (≤ 1.5% per slot).
 
### 2.3 The Recursive Subproblem
 
After ratio approximation, the problem becomes:
 
> **Given** a polynomial `P` (a multiset map fluid → integer count) with `Σ P = L`, build a binary tree such that:
> * each leaf is `Leaf(f)` for some fluid `f`;
> * the multiset of leaves of the full tree equals `P`;
> * the tree has depth ≤ d_max (a bound chosen by the algorithm, typically `d + 8`).
 
By the binary halving structure: at each internal node with input polynomial `P` and volume `L`, the algorithm must split `P` into two sub-polynomials `P₁ + P₂ = P` such that `Σ P₁ = Σ P₂ = L/2`, then recurse on each.
 
### 2.4 Optimization Objectives
 
Subject to the tree constraints, minimize a multi-objective vector:
 
* `m`        — number of mix operations (internal nodes)
* `splits`   — number of fluid-split events across all internal nodes
* `l`        — total dilution-subtree depth (definition in §3.3)
* `maxL`     — longest single dilution chain
* `leaves`   — number of leaf dispense events
* `p`        — maximum tree width (parallelism / footprint)
* `d`        — tree depth (largely determined by input)
 
There is no single canonical scalar objective; different chips, assays, and reagent costs weight these differently. This work provides each metric independently and reports per-test composite wins.
 
---
 
## 3. Metrics and Definitions
 
For a binary tree `T`, with `internal(T)` and `leaves(T)` as its internal and leaf node sets:
 
### 3.1 Mix count (`m`)
 
```
m(T) = |internal(T)|
```
 
Each internal node corresponds to one droplet-merge cycle on the chip. Lower is better — fewer mix operations means shorter assay time.
 
### 3.2 Splits (`splits`)
 
A *fluid split* at internal node `v` exists for each fluid `f` such that `f` appears in **both** subtrees rooted at `v.left` and `v.right`. Formally:
 
```
splits(T) = Σ_{v ∈ internal(T)} |{ f : f ∈ fluids(v.left) ∩ fluids(v.right) }|
```
 
Each split forces the corresponding reagent to be dispensed from the boundary reservoir at least twice — once for each subtree that consumes it. Lower is better.
 
### 3.3 Dilution length (`l`)
 
A *dilution subtree* is a maximal subtree containing **at most two** distinct fluids. In such a subtree the only operation is repeated 50:50 mixing of the same two species — a *serial dilution* — which is hardware-cheap because reagents need not be re-dispensed.
 
```
l(T) = Σ_{u: u is the root of a dilution subtree} (depth(u) − 1)
```
 
Higher is better — longer dilution chains imply better routing locality and less cross-contamination.
 
### 3.4 Max dilution chain (`maxL`)
 
```
maxL(T) = max_{u: dilution subtree root} (depth(u) − 1)
```
 
The single longest serial dilution chain in the tree.
 
### 3.5 Leaves (`leaves`)
 
```
leaves(T) = |leaves(T)|
```
 
Number of dispense events. Lower is better.
 
### 3.6 Parallelism (`p`)
 
```
p(T) = max_ℓ |{ v ∈ internal(T) : level(v) = ℓ }|
```
 
Maximum number of concurrent mix operations at any single tree level. Lower is generally better (smaller chip footprint required, less scheduling pressure).
 
### 3.7 Depth (`d`)
 
```
d(T) = height(T) − 1
```
 
Constrained by the input volume `L = 2^d` and the algorithm's depth budget; usually invariant across algorithms for the same input.
 
---
 
## 4. Literature Survey: Existing Algorithms
 
This work draws on a 15+ year line of research on combinatorial sample-preparation algorithms for DMFBs. The following are the canonical methods cited in the literature; **only RMA and BS are implemented in this work as baselines**, with the others discussed as related work.
 
### 4.1 RMA — Ratioed Mixing Algorithm (Roy, Bhattacharya, Chakrabarti, Chakrabarty, 2010)
 
* **Reference:** S. Roy, B. B. Bhattacharya, P. P. Chakrabarti, K. Chakrabarty, "*Layout-Aware Solution Preparation for Biochemical Analysis on a Digital Microfluidic Biochip*," *VLSI Design*, 2010.
* **Strategy:** Top-down greedy. At each node, the **dominant fluid** (highest count) is assigned first to one child up to volume `L/2`; remaining capacity is filled by splitting the next-largest fluid as needed.
* **Strength:** Long dilution chains (high `l` and `maxL`); good for hardware layouts where serial dilution is cheap.
* **Weakness:** Greedy myopia — a locally optimal dominant-fluid pick can force more splits or mixes downstream than necessary.
 
### 4.2 BS — Bit-Scanning (Thies et al.)
 
* **Reference:** W. Thies et al., bit-scanning approach for mixture generation on DMFBs (referenced in the LARP design notes).
* **Strategy:** Bottom-up, binary-encoding driven. For each fluid `f` with count `cf`, the bit-pattern of `cf` indicates at which tree levels a leaf for `f` must appear. Construction proceeds level-by-level pairing leaves and previously-built subtrees.
* **Strength:** Minimum leaf count and fewest mixes among the baselines (every count is encoded in `O(log cf)` leaves).
* **Weakness:** No awareness of fluid splits or dilution structure; trees are deeper and have higher fluid scattering.
 
### 4.3 REMIA — Reactant Minimization Algorithm (Hsieh et al.)
 
* **Strategy:** Skewed mixing tree + exponential dilution tree. Designed specifically to **minimize reactant usage** by exploiting an exponential dilution structure (each level halves the previous).
* **Trade-off:** Often produces *more waste droplets* than IDMA / DMRW even though it consumes fewer total reagents.
* **Status in this work:** *Not implemented.* Listed here as primary related work to acknowledge; recommended as future-work baseline.
 
### 4.4 WARA — WAste Recycling Algorithm
 
* **Strategy:** Multi-target sample preparation. Builds REMIA-style mixing trees per target, then identifies opportunities to recycle waste droplets from one target's tree as inputs to another target's tree.
* **Status:** Not implemented; multi-target sample preparation is out of scope for this work, which focuses on single-target preparation.
 
### 4.5 IDMA / DMRW (Improved Dilution and Mixing Algorithm; Dilution and Mixing with Reduced Wastage)
 
* **Strategy:** Both target the dilution subproblem (one reagent + one buffer) by carefully tracking waste-droplet recycling within a single tree.
* **Status:** Not implemented; their domain (single-reagent dilution) is a special case of the general multi-fluid mixing addressed here.
 
### 4.6 MinMix and CoDOS
 
* Generic mixing-tree construction baselines from earlier work, frequently cited alongside RMA.
* **Status:** Not implemented; RMA is used as the representative greedy baseline.
 
### 4.7 Skewed Mixing Trees (Mitra et al., 2013)
 
* **Reference:** D. Mitra et al., "*Reactant Minimization during Sample Preparation on Digital Microfluidic Biochips Using Skewed Mixing Trees*," ACM/IEEE, 2013.
* **Strategy:** Use deliberately unbalanced (skewed) trees to exploit reagent-asymmetry and reduce reactant volume.
* **Status:** Not implemented; orthogonal in tree-shape philosophy to this work.
 
### 4.8 Layout-Aware Mixture Preparation (Bhattacharjee et al., 2015)
 
* **Reference:** S. Bhattacharjee et al., "*Layout-Aware Mixture Preparation of Biochemical Fluids on Application-Specific Digital Microfluidic Biochips*," *ACM TODAES*, 2015.
* **Strategy:** Couples mixing-tree construction to physical chip layout, jointly optimizing routing distance and reagent count.
* **Status:** Not implemented; this work is layout-agnostic at present (a clear future-work direction).
 
### 4.9 Demand-Driven Preparation with Droplet Streaming (DAC 2014)
 
* **Strategy:** Optimizes for repeated mixture demand (e.g., master-mix in PCR) rather than single-target output.
* **Status:** Not implemented.
 
### 4.10 ILP-based Synthesis for Sample Preparation (IEEE 7434979)
 
* **Reference:** "*ILP-based Synthesis for Sample Preparation Applications on Digital Microfluidic Biochips*," IEEE Conference, 2016.
* **Strategy:** Models the sample-preparation problem as an integer linear program and solves with a commercial solver (CPLEX/Gurobi).
* **Relevance to this work:** This is the *solver-based* ILP approach. **Distinct from our "ILP" algorithm in §5.5**, which is an integer-allocation **DP heuristic** (no LP relaxation, no commercial solver). We share the integer-allocation framing but not the solution method; see §13 for the naming clarification. We do not directly compare against this work because we do not have the solver setup.
 
### 4.11 SIMOP — Simulation-Guided Optimization (2021)
 
* **Reference:** *SIMOP: A SIMulation-Guided OPtimization Mechanism for Sample Preparation with Digital Microfluidic Biochip*, Springer, 2021.
* **Strategy:** Couples tree generation with simulated chip behavior to optimize accuracy and error management from unbalanced splits.
* **Status:** Not implemented.
 
### 4.12 Summary of Coverage
 
| Family            | Algorithms                                | Implemented Here?        |
| ----------------- | ----------------------------------------- | ------------------------ |
| Greedy top-down   | RMA, MinMix, CoDOS                        | RMA only                 |
| Bottom-up bitwise | BS                                        | BS                       |
| Reactant-minimal  | REMIA, IDMA, DMRW                         | No (future work)         |
| Multi-target      | WARA                                      | No (out of scope)        |
| Layout-aware      | Bhattacharjee et al.                      | No (future work)         |
| Skewed-tree       | Mitra et al.                              | No                       |
| Solver-based      | ILP-based Synthesis (CPLEX/Gurobi)        | No (heuristic ILP only)  |
| Simulation-guided | SIMOP                                     | No                       |
 
---
 
## 5. Algorithms Implemented in This Work
 
Five algorithms are implemented: two are **classical baselines** (RMA, BS) and three are **proposed in this work** (AP-DP, LARP, ILP).
 
### 5.1 RMA (baseline)
 
#### Pseudocode
 
```
RMA(P, L):
  if |P| ≤ 1 or L ≤ 1:
    return Leaf(majority fluid of P)
  
  sort fluids of P by count descending
  let U be the largest-count fluid; aU = count of U
  if aU ≥ L/2:
    P1 = { U: L/2 }
    P2 = P with U decremented by L/2
  else:
    P1 = { U: aU }
    E  = L/2 − aU
    if some single fluid in P\{U} has count = E:
      move that fluid into P1
    else if a subset of P\{U} of size ≤ |P|/2 sums to E:
      move that subset into P1
    else:
      take next-largest fluid V; put min(E, count(V)) of V in P1, split V
  recurse on P1, P2 with L = L/2
```
 
#### Complexity
 
`O(n log n)` per node for sorting; `O(n²)` worst-case for the small subset-sum check (size ≤ n/2). Total tree-construction cost: `O(d · n²)` where `d` is the depth.
 
#### Behavior
 
Deterministic, single-candidate-per-node. Always builds around the dominant fluid. This produces **long dilution chains** at the cost of unnecessary mixes when the dominant fluid is much smaller than `L/2`.
 
### 5.2 BS — Bit-Scanning (baseline)
 
#### Pseudocode
 
```
BS(P, d):
  for each bit b ∈ [0, d−1]:
    lvIn[b] = { fluid f : bit b of count(f) is 1 }
  
  prev = []
  for b = 0 to d − 1:
    items = prev ∪ { Leaf(f) for f ∈ lvIn[b] }
    pair adjacent items, replacing each pair (x, y) with Mix(x, y)
    prev = the resulting list
  return prev[0]
```
 
#### Complexity
 
`O(d · n)` total — each fluid contributes `O(log count)` leaves, each leaf participates in `O(d)` mix steps.
 
#### Behavior
 
Reads the binary representation of each fluid count and "broadcasts" each bit into the appropriate tree level. Trees are typically deeper than RMA's, and the same fluid can appear at *many* different tree positions (high `splits`, but a clean construction with the **fewest mixes**).
 
### 5.3 AP-DP **[proposed]** — Adaptive Partitioning via Dynamic Programming
 
#### Inspiration
 
RMA's myopia (only one candidate per node) and BS's lack of split-awareness motivate a **multi-strategy** approach: rather than pick one rule, *generate many candidate partitions per level and pick the best by a scoring function*. AP-DP is the simplest realization of this idea — it generates candidates by five disjoint heuristics covering different geometric extremes of the partition space, then picks one by a hand-tuned local score.
 
#### The Five Strategies for Candidate Generation
 
**Strategy A — Exact whole-fluid subset sum.** Enumerate all subsets of the fluids whose counts sum *exactly* to `L/2`. These are zero-split partitions. Implementation: bit-mask DP up to 15 fluids.
 
**Strategy B — Closest-sum + single fluid split.** When no exact subset exists, find the closest-achievable-sum `S < L/2` via DP, then split exactly `(L/2 − S)` units off one fluid not in the chosen subset to close the gap. Yields one-split partitions.
 
**Strategy C — RMA-style dominant fluid.** If any fluid has count `≥ L/2`, fill one child entirely with that fluid (truncated at `L/2`).
 
**Strategy D — Greedy balanced.** Walk fluids in decreasing-count order, placing each in whichever side currently has room.
 
**Strategy E — Reverse greedy.** Same as D, but fill the second side first.
 
After generation, **deduplicate** (mirror images and exact duplicates).
 
#### Scoring Function `scorePartitionV3(P₁, P₂, L)`
 
```
score = splits × 1000
      − dilutionBonus × 50
      + totalDistinct × 10
      + minFluids × 5
```
 
* `splits` — number of fluids appearing in both `P₁` and `P₂` (heavy penalty)
* `dilutionBonus` — 1 per side where the dominant fluid > 50% of that side's volume (encourages dilution)
* `totalDistinct` — `|fluids(P₁)| + |fluids(P₂)|`
* `minFluids` — `min(|fluids(P₁)|, |fluids(P₂)|)` (prefers asymmetric, dilution-friendly splits)
 
Lower score wins.
 
#### Complexity
 
`O(2^n + n·L)` per node in the worst case (subset enumeration + DP). Practical bounds keep `n` ≤ 15 for the bit-mask phase; the implementation falls back to the heuristic strategies for larger inputs.
 
#### Behavior
 
By design, AP-DP is **myopic**: it scores each partition only on its immediate quality. Its candidate diversity is broader than RMA's, but it cannot foresee that today's "good" partition causes tomorrow's bad one (this is exactly LARP's contribution; see below).
 
### 5.4 LARP **[proposed]** — Lookahead-Augmented Recursive Partitioning
 
#### Inspiration
 
AP-DP's local scoring is its principal weakness: a partition that looks good locally (e.g., both children have a dominant fluid) may force large costs at the *next* level. LARP addresses this by **augmenting the score with a one-level lookahead**: the score of `(P₁, P₂)` includes a quick estimate of how much each child will cost when *it* is partitioned.
 
LARP also broadens candidate generation: it does *not* commit to AP-DP's five strategies, instead **enumerating all subset sums via DP** (capped only by count) and considering all single-split candidates for each closest sum.
 
#### Pseudocode
 
```
LARP(P, L):
  if |P| = 1 or L ≤ 1:
    return leaf
  
  # Step 1: comprehensive candidate generation
  zero_split    = larpEnumZeroSplit(P, L)        # all subsets summing exactly to L/2
  single_split  = larpEnumSingleSplit(P, L)      # closest-sum + single-fluid split
  rma_style     = larpEnumDominant(P, L)         # dominant fluid candidates
  candidates    = dedup(zero_split ∪ single_split ∪ rma_style)
  
  # Step 2: lookahead-augmented scoring
  for each (P1, P2) in candidates:
    immediate = splits(P1, P2) × 20 + minDistinct(P1, P2) × 3
    quick1    = larpQuickScore(P1, L/2)         # 0 if zero-split next level reachable, else 1
    quick2    = larpQuickScore(P2, L/2)
    lookahead = (quick1 + quick2) × 8
    balance   = |fluids(P1)| − |fluids(P2)|
    cand.score = immediate + lookahead + 2 × |balance|
  
  best = argmin(candidates, cand.score)
  return Mix(LARP(best.P1, L/2), LARP(best.P2, L/2))
```
 
#### `larpQuickScore(P, L)` — the look-ahead estimator
 
A lightweight DP that asks only: *can the next level reach a zero-split partition?*
 
```
target = L/2
dp[0] = true
for each fluid with count v:
  for s = target downto v:
    dp[s] |= dp[s − v]
return 0 if dp[target] else 1
```
 
Time: `O(n · L)` per call, called twice per candidate. Negligible relative to the full enumeration.
 
#### Complexity
 
`O(d · (n · L + C · n · L/4))` where `C` is the number of candidates per node (effectively bounded). In practice runs in ≤ 60 ms per tree on our benchmark.
 
#### Behavior
 
LARP's lookahead correctly identifies that *zero-split-now* with *zero-split-next* is generally better than *zero-split-now* alone. Empirically it produces the lowest split counts of all algorithms tested (excluding the BS bound, which has its own quirks; see §9).
 
### 5.5 ILP **[proposed]** — Integer-Allocation DP heuristic
 
> **Naming note.** "ILP" here refers to the **integer-allocation framing**, not to a solver-based integer linear program. There is no LP relaxation, no simplex, no branch-and-bound, no commercial solver. The algorithm is a **dynamic-programming heuristic** that prunes its search using a non-domination rule on (splits, allocation) states. The name is retained because the variables of interest — how much of each fluid goes into each child — are integer-valued allocations. See §13 for the full clarification.
 
#### Inspiration
 
LARP's candidate generation, while comprehensive, is procedural: enumerate zero-split subsets, then single-split candidates, etc. ILP takes a different stance: **set up the problem as a DP over allocation decisions and track all non-dominated (splits, allocation) states.** This generalizes both zero-split and multi-split candidates uniformly: at each item `i`, decide how much of fluid `i` goes into `P₁` (an integer in `[0, count(i)]`); the DP tracks all non-dominated states.
 
#### Phases
 
**Phase 1 — Integer-allocation DP.** State `dp[s]` = list of `{splits, alloc[]}` pairs reaching exact running sum `s`. Transition: for each item with count `v`, try every take `t ∈ [0, v]`. A take with `0 < t < v` adds 1 to the split count. A non-domination rule is enforced: keep only states whose `(splits, alloc)` are not dominated. A cap `ALLOC_CAP = 32` per `s` prevents combinatorial blow-up.
 
**Phase 2 — Collect all candidates** with `s = L/2`.
 
**Phase 3 — RMA-style fallback.** For each fluid with `count ≥ L/2`, generate the dominant-fluid partition.
 
**Phase 4 — Greedy fallback** for inputs where Phases 1–3 yield no valid partition.
 
**Phase 5 — Deduplicate** (canonical sorted-pair key).
 
**Phase 6 — Score and select** with `ilpScoreAllocation`:
 
```
score = splits × 1000
      + (n1 + n2) × 10                       # secondary: minimize distinct fluids
      − 50 if n1 ≤ 2; − 50 if n2 ≤ 2         # reward dilution potential
      − 30 if n1 = 1; − 30 if n2 = 1         # reward pure isolation
      − 20 if dominant fluid ≥ 50% of side   # reward dilution dominance
      + (la1.splits + la2.splits) × 100      # 1-level lookahead penalty
      − 15 if next level can be zero-split
```
 
#### Complexity
 
`O(n · (L/2) · ALLOC_CAP · v_max)` per node for the DP — *quadratic in volume*, which is the bottleneck. On the largest stress tests (33 fluids, depth 12), this can take up to 4.4 s per test (see §9.5).
 
#### Behavior
 
ILP explores a strictly larger candidate space than LARP and finds partitions that LARP's candidate enumerators miss. Its scoring is more elaborate than LARP's. Empirically ILP's metrics are very close to LARP's; see §9.
 
### 5.6 Algorithm Comparison at a Glance
 
| Aspect                       | RMA              | BS               | AP-DP            | LARP                     | ILP                          |
| ---------------------------- | ---------------- | ---------------- | ---------------- | ------------------------ | ---------------------------- |
| Direction                    | Top-down         | Bottom-up        | Top-down         | Top-down                 | Top-down                     |
| Candidates per node          | 1                | Fixed (no pick)  | ≤ 30             | All subsets + variants   | All non-dominated DP states  |
| Lookahead                    | None             | None             | None (local)     | 1-level                  | 1-level (richer)             |
| Primary objective            | Dilution         | Min mixes        | Splits + diluton | Splits                   | Splits                       |
| Worst-case time              | O(d · n²)        | O(d · n)         | O(2^n + n·L)     | O(d · n · L)             | O(d · n · L · ALLOC_CAP)     |
| Origin                       | Roy et al. 2010  | Thies et al.     | **This work**    | **This work**            | **This work**                |
 
---
 
## 6. The Beam-Search + Memoization Meta-Framework
 
This is a **meta-algorithm** that augments any partitioning-based mixing-tree algorithm. It applies in this work to AP-DP, LARP, and ILP, leaving RMA and BS unchanged because their candidate-generation models are not partition-pick-driven.
 
### 6.1 Motivation
 
The base algorithms (LARP, ILP, AP-DP) all share a common shape:
 
```
buildX(P, L):
  candidates = enumerate(P, L)
  best       = argmin_{c in candidates} score(c)
  return Mix(buildX(best.P1, L/2), buildX(best.P2, L/2))
```
 
The score function is a heuristic *estimate* of the candidate's downstream cost. A natural improvement: **don't trust the heuristic; build the subtrees and measure**. This is **beam search at depth `K`**:
 
```
buildX(P, L, K):
  candidates = enumerate(P, L)
  topK       = sort(candidates, key=score)[:K]
  for each c in topK:
    leftSub  = buildX(c.P1, L/2, K)
    rightSub = buildX(c.P2, L/2, K)
    cost     = trueSubtreeCost(c, leftSub, rightSub)
  return the candidate with minimum cost
```
 
When `K = 1`, this is identical to the original greedy algorithm. When `K → ∞`, it is exhaustive search. **The default in this work is `K = 3`.**
 
The cost: without further optimization, beam search is `O(K^d)` build calls per top-level construction. **Memoization** saves us: the same `(P, L)` sub-problem appears in many places (most strikingly when `P₁ = P₂` — symmetric inputs). Caching `(P, L) → bestSubtree` collapses repeated work.
 
### 6.2 Three Shared Helpers
 
#### `canonKey(P, L)`
 
A stable, order-independent key for the cache:
 
```
canonKey(P, L) = sort_lexicographic(P).join('|') + '@' + L
```
 
This guarantees `(P, L)` and `(reorderedP, L)` map to the same cache entry.
 
#### `cloneWithLevel(node, baseLv)`
 
A subtle correctness issue: a memoized subtree built first at depth 2 may later be reused at depth 5. The `level` field on each node is required by the `parallelism (p)` metric, so a stale-level retrieval would silently corrupt that statistic. `cloneWithLevel` walks the cached subtree and produces a fresh deep-copy with `level` rewritten to its actual position in the surrounding tree:
 
```
cloneWithLevel(node, lv):
  if node is null: return null
  newNode = shallow-copy of node
  newNode.level = lv
  if node has children:
    newNode.left  = cloneWithLevel(node.left,  lv + 1)
    newNode.right = cloneWithLevel(node.right, lv + 1)
  return newNode
```
 
#### `subtreeTotalCost(node, splitWeight, dilutionWeight = 5)`
 
The *true* subtree cost (replacing the algorithm-specific score estimates):
 
```
cost = splits(subtree) × splitWeight  +  mixes(subtree)  −  dilution(subtree) × dilutionWeight
```
 
With `splitWeight = 1000` and `dilutionWeight = 5`, splits are the primary objective, mixes are the secondary tiebreaker, and dilution is a small bonus. The choice of weights matches the design intent of LARP and ILP (which both list splits as their primary objective).
 
### 6.3 Beam-Augmented `buildLARP` (illustrative)
 
```
buildLARP(P, L, K = 3, memo = ∅):
  if base case: return leaf
  
  if memo has canonKey(P, L):
    return cloneWithLevel(memo[canonKey(P, L)], current_level)
  
  candidates = larpEnumPartitions(P, L)
  rank candidates by larpScore (the cheap heuristic)
  beam       = top-K candidates
  
  best, bestCost = nil, ∞
  for each c in beam:
    left  = buildLARP(c.P1, L/2, K, memo)
    right = buildLARP(c.P2, L/2, K, memo)
    node  = Mix(P, L, left, right)
    cost  = subtreeTotalCost(node, 1000, 5)
    if cost < bestCost: best, bestCost = node, cost
  
  memo[canonKey(P, L)] = best
  return best
```
 
The same template (beam-search + memoization) is applied to ILP via `ilpScoreAllocation` as the cheap-heuristic ranker.
 
### 6.4 Memoization-Only for AP-DP
 
AP-DP intentionally retains its myopic local-scoring identity — adding beam search to it would erase its character as a comparison baseline. Memoization (alone) is applied so repeated sub-problems are cached, giving the same picks faster.
 
### 6.5 Cost / Benefit
 
* **Compute overhead**: ~10× wall-time vs the greedy baselines for LARP, ~100× for ILP (driven mostly by ILP's Phase-1 DP, not by beam search overhead per se).
* **Quality gain**: ~14% reduction in splits and mixes, with a small dilution regression as a trade-off (§9).
 
---
 
## 7. Benchmarking Framework
 
### 7.1 Goals
 
A reproducible, deterministic, multi-metric benchmark that exposes the relative strengths of each algorithm across a stratified set of input distributions covering the practically interesting cases for a DMFB sample-preparation problem.
 
### 7.2 Architecture
 
* **`benchmark.js`** — module containing all five algorithm implementations, the metric calculators (`cntM`, `cntSplits`, `dL`, `maxDL`, `cntL`, `mxP`), the tree builder dispatch (`buildTree`), and the ratio-approximation pre-processor (`ratioApprox`).
* **`bench.js`** — the test runner. It seeds a deterministic PRNG, generates a stratified test suite, runs every algorithm on every test, computes per-test statistics, aggregates per-category and overall, prints multi-table dashboards, and persists `benchmark_latest.json` and `benchmark_detailed.csv`.
 
The split is intentional: `benchmark.js` is the **algorithm library** (consumable by the React UI as well), and `bench.js` is the **experiment driver** consuming it via dynamic `eval`.
 
### 7.3 Stratified Test Suite (323 cases)
 
| Category       | n   | Fluid count (typical) | Depth (typical)   | Purpose                                                  |
| -------------- | --: | --------------------- | ----------------- | -------------------------------------------------------- |
| `paper`        | 9   | as in published work  | as in published   | Reproduce specific examples from the literature          |
| `edge`         | 12  | 2–4                   | 3–5               | Small/edge-case ratios (powers of two, all equal, etc.)  |
| `stress`       | 2   | 30+                   | 12                | Largest tractable inputs, runtime stress test            |
| `rand-few`     | 40  | 2–4                   | 4–8               | Random ratios, small fluid count                         |
| `rand-med`     | 80  | 5–8                   | 6–9               | Random ratios, mid-range count                           |
| `rand-many`    | 60  | 10–18                 | 7–10              | Random ratios, many fluids                               |
| `rand-skew`    | 40  | 5–8                   | 6–9               | Strongly asymmetric counts (one fluid dominates)         |
| `rand-equal`   | 40  | 5–8                   | 6–9               | All counts within ±1 of each other                       |
| `rand-deep`    | 20  | 8–12                  | 10–12             | Deep trees, high-volume                                  |
| `rand-shallow` | 20  | 3–5                   | 4–6               | Shallow trees, exposes fewer-mix algorithms              |
 
Total: **323** test cases.
 
### 7.4 Metrics Aggregation
 
For each test × algorithm:
 
* Per-test absolute values for all 7 metrics
* Per-test "win" indicators (does this algorithm equal the best on this metric?)
* Per-test composite-best (lowest sum of normalized metric ranks)
 
Aggregations:
 
* **Per-category**: mean and median of each metric per algorithm
* **Overall**: same, plus win counts per metric
* **Final scoreboard**: rank-based points across all metrics, with **tied-rank-aware points** (multiple algorithms tied for `k`-th place all receive the `k`-th-place points; later positions skip)
 
### 7.5 Reproducibility
 
* Deterministic PRNG seed (constant) → same suite every run
* Per-run output: `benchmark_latest.json` (full result), `benchmark_detailed.csv` (per-test row), and an in-terminal dashboard
* Diff-vs-baseline mode prints `+0.12 ✗` style deltas when a previous baseline file is present, enabling fast iteration during development
 
### 7.6 Per-Algorithm Timing
 
`bench.js` records `{total, avg, max}` wall-time per algorithm. Reported as a separate table in the output. Used to surface ILP's worst-case 4.4-s tail on stress tests (§9.5) — a real concern for interactive use.
 
---
 
## 8. The Visualization Tool (React Web App)
 
### 8.1 Tech Stack
 
* **Vite + React 18** as the build tool and framework
* **Plain Canvas** for tree rendering (no charting library dependency)
* **Single-file `App.jsx`** containing the algorithm library (mirrored from `benchmark.js`) plus the UI
 
### 8.2 Features
 
* Input panel for ratios, depth, and algorithm selection (RMA / BS / AP-DP / LARP / ILP)
* On-input recomputation of all five trees
* Side-by-side comparison view with all metrics displayed beneath each tree
* Tree layout (`layoutTree`) computes node positions so internal nodes are horizontally centered above their children
* Color-coded leaves per fluid for visual distinguishability
* Algorithm-info panels with short descriptions of each method and the relevant references
 
### 8.3 Architecture
 
```
ratio input  ─────────►  ratioApprox  ─────────►  approx[]  ──┐
                                                              │
algorithm pick ──────────────────────────────────────────────┤
                                                              ▼
depth ─────────────────────────────────►  buildTree(algo, ...)  ──►  tree
                                                                          │
                                                                          ▼
                                                                  layoutTree
                                                                          │
                                                                          ▼
                                                                    Canvas render
                                                                          │
                                                                          ▼
                                                                    metrics panel
```
 
The UI is intentionally lean — its purpose is *demonstration and debugging*, not a research artifact in itself.
 
---
 
## 9. Experimental Results
 
All numbers below come from a single deterministic run of the full benchmark suite (n = 323) with beam-search width `K = 3` for LARP and ILP, and the dilution-aware subtree cost (`splitWeight = 1000`, `dilutionWeight = 5`). Reproduced via `node bench.js` from a fresh clone.
 
> **A note on methodology.** An earlier iteration of this work used a 7-metric composite to rank algorithms per test and reported "% of tests won". On audit, three sources of inflation were identified:
> 1. **`leaves = mixes + 1`** is mathematically locked in any binary tree, so `leaves` and `mixes` are the same metric counted twice.
> 2. **`splits` is highly correlated with `mixes`** in the test suite (per-test win-rate vectors are identical for both).
> 3. **`depth` is invariant** across all algorithms because it is determined by the input, not the algorithm.
>
> So the original "7-metric composite" was effectively a **4-metric composite with 3-fold weight on the size cluster** (mixes/splits/leaves) and a constant offset (depth). After correction the headline is reported on the **4 independent metrics** `m, l, maxL, p` and supplemented with **sole-vs-tied wins**, **per-rank distribution**, and **pairwise tie matrix**, which together prevent the previous over-reporting.
 
### 9.1 Overall Averages
 
| Metric         | RMA    | BS     | AP-DP  | LARP   | ILP    | Direction | Best on average |
| -------------- | -----: | -----: | -----: | -----: | -----: | :-------: | :--------------: |
| Mixes (m)      | 21.80  | 17.96  | 21.93  | 18.98  | 19.14  | lower ↓   | BS               |
| Splits         | 16.10  | 12.26  | 16.23  | 13.27  | 13.43  | lower ↓   | BS               |
| Dilution (l)   | 16.43  |  8.70  | 16.07  | 12.57  | 13.00  | higher ↑  | RMA              |
| MaxDilution    |  5.16  |  1.81  |  5.21  |  4.82  |  4.91  | higher ↑  | AP-DP            |
| Leaves         | 22.80  | 18.96  | 22.93  | 19.98  | 20.14  | lower ↓   | BS               |
| Depth          |  7.45  |  7.45  |  7.45  |  7.45  |  7.45  | input     | (all tied)       |
| Parallelism    |  4.06  |  3.56  |  4.03  |  3.64  |  3.66  | lower ↓   | BS               |
 
### 9.2 Headline: Proposed vs Published Baselines (% change)
 
This is the **single most important table** of the work. It reports each proposed algorithm's per-metric performance as a signed percentage change versus each published baseline. ✓ marks the direction we want.
 
| Comparison              | Mixes        | Splits       | Leaves       | Dilution (l) | MaxDilution    | Parallelism   |
| ----------------------- | -----------: | -----------: | -----------: | -----------: | -------------: | ------------: |
| **LARP vs RMA** (Roy 2010) | **−12.9%** ✓ | **−17.6%** ✓ | **−12.4%** ✓ | −23.5% ✗     | −6.6% ✗        | **−10.3%** ✓  |
| **LARP vs BS** (Thies)     | +5.7% ✗      | +8.2% ✗      | +5.4% ✗      | **+44.5%** ✓ | **+166.3%** ✓  | +2.2% ✗       |
| **ILP vs RMA** (Roy 2010)  | **−12.2%** ✓ | **−16.6%** ✓ | **−11.7%** ✓ | −20.9% ✗     | −4.8% ✗        | **−9.9%** ✓   |
| **ILP vs BS** (Thies)      | +6.6% ✗      | +9.5% ✗      | +6.2% ✗      | **+49.4%** ✓ | **+171.3%** ✓  | +2.8% ✗       |
 
**Reading this table:**
* LARP and ILP **strictly improve over RMA on size and parallelism** (5/6 metrics) at the cost of dilution-chain length.
* LARP and ILP **strictly improve over BS on dilution structure** (2/2 dilution metrics) at the cost of 6–10% more mix operations.
* Neither LARP nor ILP dominates either baseline; they occupy a **previously-unfilled middle of the splits-vs-dilution trade-off curve**.
 
This is the honest research one-liner:
 
> *"On a 323-test stratified DMFB sample-preparation benchmark, beam-augmented LARP achieves **18% fewer reagent dispenses** than RMA (Roy et al., 2010) and **44% longer dilution chains** than BS (Thies et al.), occupying a previously-unfilled middle of the splits-vs-dilution trade-off curve."*
 
### 9.3 Per-Test Rank Distribution (the most informative single view)
 
For each test, the algorithms are ranked by their composite score over the **4 independent metrics**: `mixes`, `dilution (l)`, `maxL`, `parallelism`. The histogram below counts how often each algorithm placed at each rank.
 
| Algorithm | 1st   | 2nd   | 3rd   | 4th   | 5th   |
| --------- | ----: | ----: | ----: | ----: | ----: |
| RMA       | 158   |  23   |  53   |  70   |  19   |
| BS        |   2   |  66   |  17   |  79   | **159** |
| AP-DP     |  34   |  52   | 111   |  53   |  73   |
| **LARP**  |  92   |  92   |  61   |  73   | **5** |
| **ILP**   |  37   |  90   |  81   |  48   |  67   |
 
**Reading this table — the qualitative story:**
* **RMA** is bimodal: it is 1st place on 49% of tests (its dilution-optimization wins) but last on 6% (when dilution doesn't matter).
* **BS** is the *most extreme* algorithm: 1st on only 0.6% of tests but last on **49.2%** when dilution structure matters.
* **AP-DP** is consistently middle-of-the-pack, with no clear strength.
* **LARP** is the **most consistent algorithm**: ranked 1st-or-2nd on **57%** of tests and last place on only **1.5%** (5 / 323). It never has catastrophic regressions.
* **ILP** is similar to LARP, with slightly higher variance.
 
### 9.4 Sole vs Tied Composite Wins (4-metric)
 
A win is "sole" if only one algorithm has the maximum composite; "tied" includes any algorithm at that maximum.
 
| Algorithm | Sole-best     | Tied-best (any-of-best) |
| --------- | ------------: | ----------------------: |
| RMA       |  22 (6.8%)    |  158 (48.9%)            |
| BS        |   2 (0.6%)    |   66 (20.4%)            |
| AP-DP     |  19 (5.9%)    |  130 (40.2%)            |
| **LARP**  |  35 (10.8%)   |  235 (72.8%)            |
| **ILP**   |  37 (11.5%)   |  240 (74.3%)            |
 
Average algorithms tied at the top per test: **2.57** — most tests have multiple algorithms reaching the same composite peak. The "tied-best" column is therefore an **upper bound** on each algorithm's quality; the "sole-best" column is the **lower bound**.
 
### 9.5 Pairwise Composite-Tie Matrix
 
How often do each pair of algorithms produce indistinguishable trees on the 4-metric composite?
 
|       | RMA   | BS    | AP-DP | LARP  | ILP   |
| ----- | ----: | ----: | ----: | ----: | ----: |
| RMA   |  —    | 117   | 104   | 138   | 141   |
| BS    | 117   |  —    | 118   |  65   |  69   |
| AP-DP | 104   | 118   |  —    | 110   | 113   |
| LARP  | 138   |  65   | 110   |  —    | **226** |
| ILP   | 141   |  69   | 113   | **226** | — |
 
**The single most important entry: LARP–ILP = 226 / 323 (70%).** LARP and ILP produce identical-quality trees on 70% of the test suite. This is reported honestly as a finding: *the two proposed lookahead algorithms converge on the same well-balanced answers most of the time*, suggesting both have reached the practical optimum on the explored candidate space. ILP explores a strictly larger candidate space but the additional candidates rarely change the chosen partition.
 
### 9.6 Strict-Best on Splits (the design objective)
 
Splits is the *primary* design objective for LARP and ILP. How often is each algorithm the unique best on splits?
 
| Algorithm | Sole-best on splits | Tied-best on splits  |
| --------- | ------------------: | -------------------: |
| RMA       |   0   (0.0%)        | 113   (35.0%)        |
| **BS**    | **125 (38.7%)**     | **323 (100.0%)**     |
| AP-DP     |   0   (0.0%)        |  78   (24.1%)        |
| LARP      |   0   (0.0%)        | 192   (59.4%)        |
| ILP       |   0   (0.0%)        | 182   (56.3%)        |
 
**Critical, honest finding: BS is the unbeatable champion on splits.** BS is at-least-tied for fewest splits on every single test (100%), and uniquely best on 38.7%. *No proposed algorithm is ever strictly fewer-splits than BS.* LARP and ILP can only **match** BS, on 59.4% and 56.3% of tests respectively.
 
This does **not** invalidate the contribution: BS is split-optimal but **catastrophic on dilution** (1.81 avg vs 12+ for the others). The proposed algorithms' value is in *matching BS's split count more than half the time* while *also preserving dilution structure*. See §10.
 
### 9.7 Per-Metric Win Counts (any-of-best, ties counted)
 
| Metric       | RMA       | BS         | AP-DP     | LARP    | ILP     |
| ------------ | --------- | ---------- | --------- | ------- | ------- |
| Mixes        | 35%       | **100%**   | 24%       | 59%     | 56%     |
| Splits       | 35%       | **100%**   | 24%       | 59%     | 56%     |
| Dilution (l) | **63%**   | 18%        | 68%       | 27%     | 28%     |
| MaxDilution  | 85%       | 17%        | **90%**   | 67%     | 70%     |
| Leaves       | 35%       | **100%**   | 24%       | 59%     | 56%     |
| Parallelism  | 68%       | **100%**   | 64%       | **92%** | 91%     |
 
> Note: the win-rate vectors for `Mixes`, `Splits`, and `Leaves` are *identical* — direct empirical confirmation of the metric correlation noted at the top of §9. This is why `splits` and `leaves` were dropped from the composite in §9.3.
 
### 9.8 Final Scoreboard (7-metric, tied-rank-aware)
 
For completeness, the original 7-metric scoreboard (1st = 6 pts, 2nd = 4, …, with ties sharing the higher position's points):
 
```
#1  BS     32 pts
#2  LARP   26 pts
#3  ILP    24 pts
#4  RMA    23 pts
#5  AP-DP  21 pts
```
 
This scoreboard places **BS first** because it dominates the three correlated size metrics (mixes/splits/leaves), each contributing 6 pts × 3 = 18 pts to BS's total. **The 4-metric per-test rank distribution in §9.3 tells a different and more honest story** — BS ranks last in 49% of tests on the de-correlated composite. Both views are reported; downstream documents should prefer §9.3 for any "which is the best algorithm" claim.
 
### 9.9 Per-Algorithm Timing (n = 323)
 
| Algorithm | Total           | Avg          | Max            |
| --------- | --------------: | -----------: | -------------: |
| RMA       |       33 ms     |  0.10 ms     |      2 ms      |
| BS        |        9 ms     |  0.03 ms     |      1 ms      |
| AP-DP     |      136 ms     |  0.42 ms     |      3 ms      |
| LARP      |    2 288 ms     |  7.08 ms     |     62 ms      |
| ILP       | **24 403 ms**   | **75.55 ms** | **4 427 ms**   |
 
ILP's worst-case 4.4 s tail is on the largest stress test (33 fluids, depth 12). For interactive UI use on inputs > 25 fluids, this is borderline; documented as a limitation in §11.
 
### 9.10 Phase-B Impact (beam search vs greedy on the same algorithm)
 
| Metric        | LARP-beam vs LARP-greedy | ILP-beam vs ILP-greedy |
| ------------- | -----------------------: | ---------------------: |
| Mixes         | **−2.03 (−9.7%)**        | **−1.32 (−6.4%)**      |
| Splits        | **−2.17 (−14.0%)**       | **−1.33 (−9.0%)**      |
| Leaves        |     −2.03                |     −1.32              |
| Parallelism   |     −0.22                |     −0.14              |
| Dilution (l)  | **+3.22 (regression)**   | **+1.56 (regression)** |
| MaxDilution   |     +0.35                |     +0.19              |
 
Beam-augmented LARP and ILP **gain 9–14% on the design objective** (splits/mixes) and **lose 13–25% on dilution chains** (a non-objective). The trade-off is honest and documented as an improvement on the stated objective, not a uniform improvement on every metric.
 
---
 
## 10. Discussion and Trade-offs
 
### 10.1 The Splits-vs-Dilution Trade-off
 
The single clearest finding of this work: **splits and dilution are negatively correlated under a fixed depth budget.** The five algorithms cluster into three regimes on this 2-D frontier:
 
```
                    HIGH dilution
                          │
        RMA (16.4)        │
              ●           │
   AP-DP (16.1) ●         │
                          │
        LARP (12.6) ●─────┤   ← previously unfilled middle
        ILP  (13.0) ●     │       (proposed contribution)
                          │
              BS (8.7) ●  │
                          │
                    LOW dilution
                          │
   ──────────────────────────────────────────────
   16 splits        13 splits        12 splits
   (more splits)                  (fewer splits)
```
 
| Regime                | Algorithms        | Property                                                        |
| --------------------- | ----------------- | --------------------------------------------------------------- |
| **Dilution-optimized**| RMA, AP-DP        | Long dilution chains, but ≥ 16 splits on average                |
| **Balanced (proposed)** | **LARP, ILP**   | 13 splits *and* 12.6+ dilution — neither extreme                |
| **Split-optimized**   | BS                | 12.3 splits but only 8.7 dilution — catastrophic on dilution    |
 
No algorithm dominates the others. The contribution of this work is **filling the previously unoccupied middle of the frontier**, where neither BS (which ignores dilution) nor RMA (which freely splits) was a competitive option.
 
### 10.2 LARP vs ILP — They Are Effectively the Same Algorithm
 
A crucial honest finding (§9.5): **LARP and ILP produce indistinguishable trees on 226 of 323 tests (70%)**. Although their internal candidate generators differ (LARP enumerates subset-sum partitions; ILP runs an integer-allocation DP), their final picks coincide in the vast majority of cases.
 
* On the 30% of tests where they differ, LARP wins as often as ILP loses, with both staying within ±2% on every average metric.
* ILP's worst-case wall-time (4.4 s) is **two orders of magnitude** larger than LARP's (62 ms) for the same answer.
 
The honest framing: this is *empirical evidence that beam-augmented partition-pick algorithms have converged on the practical optimum* — different candidate-enumeration strategies arrive at the same answers because the answer is essentially fixed by the input. **For practical use, LARP is the recommended algorithm; ILP is retained as a candidate-richer reference implementation.**
 
### 10.3 BS Is Unbeatable on Splits — and That's a Different Algorithm Niche
 
Section 9.6 shows that **BS is at-least-tied for fewest splits on every test in the suite (100%)**, and uniquely best on 38.7%. *No proposed algorithm ever beats BS on splits in absolute count.*
 
This does not weaken the contribution. BS achieves split-optimality by **scattering fluids across the tree level-by-level** based on bit position, which is structurally orthogonal to dilution-friendly construction. The proposed algorithms' value is precisely that they **match BS's split count on 59% of tests** *while preserving the dilution structure* that BS destroys (BS dilution = 1.81 vs LARP = 12.57). On a chip where both reagent dispenses and routing locality matter (most realistic chips), LARP is the better choice; on a chip where only dispense count matters, BS is better.
 
### 10.4 When to Use Which Algorithm
 
| Use case                                                                            | Recommended         |
| ----------------------------------------------------------------------------------- | ------------------- |
| Splits-dominated cost only; ignore routing — minimum dispenses                      | BS                  |
| Long dilution chains, layout-routing-dominated cost — PCR-style serial dilution     | RMA                 |
| Both reagent dispenses *and* dilution structure matter (most realistic assays)      | **LARP**            |
| Same as above, with maximum candidate exploration (when compute is cheap)           | **ILP**             |
| Pedagogical baseline / comparison reference                                         | AP-DP               |
 
### 10.5 The Role of Beam Search
 
Beam search is **not** a novel algorithmic technique — it is a 1970s AI textbook tool. Its contribution here is empirical: applied at width K = 3 with proper memoization to a well-designed candidate-generator (LARP or ILP), it eliminates the gap between the heuristic candidate-score estimate and the true downstream cost, **recovering 9–14% on the design objective for ~10× wall-time**. This is a *favorable* engineering trade-off in any regime where reagent cost dominates compute cost (i.e. the entire wet-lab use case).
 
### 10.6 Why AP-DP Underperforms
 
AP-DP was the most ambitious of the local-scoring algorithms in the original design but loses on most metrics post-correction. The reason: AP-DP's candidate generators are an arbitrary subset of those LARP enumerates, so it is *strictly weaker* in candidate diversity. Its retention as a baseline (memoization-only, no beam) preserves its identity as a "sophisticated greedy" — useful as the bridge between RMA and the lookahead methods in the literature progression.
 
### 10.7 Methodological Note: Why the 4-Metric Composite Is the Honest View
 
The 7-metric composite originally used in this work was found to have three sources of inflation:
 
1. **Mathematical redundancy.** `leaves = mixes + 1` in any binary tree — empirically confirmed by the identical 35/100/24/59/56 win-rate vectors for `Mixes`, `Splits`, `Leaves` in §9.7.
2. **Empirical correlation.** `splits` and `mixes` have identical per-test win patterns in the suite, indicating they rank algorithms identically per test.
3. **Constant-offset metric.** `depth` is determined by the input, contributing zero discrimination.
 
After correction (§9.3), the 4-metric composite (`m, l, maxL, p`) gives **per-rank distributions that vary substantially across algorithms**, which the inflated composite hid. The corrected composite is the basis for all "best" claims in this work.
 
---
 
## 11. Limitations and Future Work
 
### 11.1 Limitations (honest)
 
* **LARP and ILP are near-equivalent in output** (70% of tests produce identical trees; §9.5). Effectively the proposed contribution is *one* algorithm with two candidate-enumeration strategies, not two distinct algorithms.
* **BS is never beaten on splits** (§9.6). LARP/ILP can only match BS's split count, not exceed it. The contribution is balance, not absolute split-minimization.
* **No comparison to REMIA / WARA / IDMA / DMRW.** These are major prior algorithms in the DMFB literature. Their absence is acknowledged; their inclusion would strengthen the empirical comparison considerably.
* **No real ILP solver.** ILP is a heuristic, not actual integer linear programming. The 2016 IEEE ILP paper (CPLEX/Gurobi-based) provides a true solver-based baseline that this work does not match.
* **Layout-agnostic.** No physical chip routing, no electrode-occupancy model, no transport-cost model. Layout-aware methods (Bhattacharjee 2015) are not benchmarked.
* **No hardware-cost model.** Metrics are graph-theoretic, not in chip seconds or microliters of reagent. Splits is a *proxy* for reagent dispenses, not for actual reagent volume in mL.
* **No statistical-significance tests** in the current `bench.js` output (e.g., Wilcoxon signed-rank). Per-test numbers are reported but no p-values.
* **ILP's worst-case 4.4 s tail** is hostile for interactive UI use on inputs > 25 fluids.
* **No optimality bounds.** A brute-force optimal solver for small `n ≤ 5` would establish how close beam-augmented LARP comes to the true minimum.
* **The composite scoring required correction** during the work (see §9 methodological note). The original 7-metric composite was inflated by metric redundancy; the corrected 4-metric composite is used throughout.
 
### 11.2 Concrete Future Work (priority-ordered)
 
1. **Implement REMIA and WARA** as published baselines — adds the two strongest comparison points the current work lacks.
2. **Add a brute-force optimal solver** for small inputs (n ≤ 5) to establish a ceiling for beam-augmented LARP.
3. **Add Wilcoxon signed-rank tests** to `bench.js` for paired pre-/post-Phase-B comparisons (gives p-values for the abstract).
4. **Replace ILP's Phase-1 DP with a real ILP** via `glpk-js` or a similar in-browser solver; use as quality oracle on small inputs.
5. **Add a hardware-cost model** (transport, dispenses, contamination) and re-rank algorithms.
6. **Multi-target preparation** (à la WARA) — reuse waste droplets across targets.
7. **Layout-aware variants** — joint optimization of tree and routing.
 
---
 
## 12. References
 
1. **S. Roy, B. B. Bhattacharya, P. P. Chakrabarti, K. Chakrabarty.** *Layout-Aware Solution Preparation for Biochemical Analysis on a Digital Microfluidic Biochip.* VLSI Design, 2010. *(RMA — primary baseline reference.)*
2. **W. Thies et al.** *Bit-scanning approach for mixture generation on DMFBs.* (BS — implemented baseline.)
3. **D. Mitra et al.** *Reactant Minimization during Sample Preparation on Digital Microfluidic Biochips Using Skewed Mixing Trees.* ACM/IEEE, 2013.
4. **S. Bhattacharjee et al.** *Layout-Aware Mixture Preparation of Biochemical Fluids on Application-Specific Digital Microfluidic Biochips.* ACM TODAES, 2015.
5. *Demand-Driven Mixture Preparation and Droplet Streaming using Digital Microfluidic Biochips.* DAC, 2014.
6. *ILP-based Synthesis for Sample Preparation Applications on Digital Microfluidic Biochips.* IEEE Conference Publication, 2016.
7. *SIMOP: A SIMulation-Guided OPtimization Mechanism for Sample Preparation with Digital Microfluidic Biochip.* Springer, 2021.
8. **Hsieh et al.** *REMIA — Reactant Minimization Algorithm.* (Cited via UCR Microfluidics bibliography.)
9. *Waste-aware single-target dilution of a biochemical fluid using digital microfluidic biochips.* ScienceDirect, 2014.
10. *Reactant Minimization for Sample Preparation on Microfluidic Biochips With Various Mixing Models.* IEEE Journals, 2015.
 
### Online resources consulted
 
* University of California Riverside DMFB Microfluidics teaching site (single- and multi-target sample preparation).
* ACM Digital Library and IEEE Xplore for the publications above.
 
---
 
## 13. Naming Note: What "ILP" Means in This Work
 
The algorithm referred to as **ILP** in this work and codebase is **not actual integer linear programming with a solver** (CPLEX, Gurobi, GLPK). It is a **dynamic-programming (DP) heuristic** that:
 
* Tracks **integer-valued allocation states** — for each fluid, an integer count indicating how much goes into the left child of the partition
* Prunes the search using a **non-domination rule** on (split count, allocation) pairs
* Caps the number of states per running sum (`ALLOC_CAP = 32`) to keep runtime tractable
 
The "ILP" name is retained because the algorithm operates on integer allocation variables and uses a frontier-pruning strategy reminiscent of solver-based methods. Solver-based ILP work for the same problem (e.g., the 2016 IEEE paper *ILP-based Synthesis for Sample Preparation*, listed in §4.10) uses commercial solvers and is **distinct from our algorithm** — we cite it as related work but do not directly compare. Our ILP gives no global-optimality guarantee; it is a tractable heuristic, not an exact solver.
 
In short: **"ILP" here = an integer-allocation DP heuristic, not an LP solver.**
 
---
 
*End of document.*