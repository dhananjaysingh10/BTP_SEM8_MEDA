# Mixing Tree Construction for Sample Preparation on Digital Microfluidic Biochips

> **B.Tech Project (BTP) — Semester 8**  
> A comprehensive study of five partitioning algorithms with beam-search augmentation for constructing optimal mixing trees on DMFBs.

---

## 📋 Table of Contents

1. [About the Project](#about-the-project)
2. [Key Features](#key-features)
3. [Project Structure](#project-structure)
4. [Prerequisites](#prerequisites)
5. [Installation & Setup](#installation--setup)
6. [Running the Visualization Tool](#running-the-visualization-tool)
7. [Running the Benchmark Suite](#running-the-benchmark-suite)
8. [Compiling the LaTeX Report](#compiling-the-latex-report)
9. [Algorithms Implemented](#algorithms-implemented)
10. [Metrics Evaluated](#metrics-evaluated)
11. [Visualization Tool Features](#visualization-tool-features)
12. [Benchmark Framework Features](#benchmark-framework-features)
13. [Key Results](#key-results)
14. [References](#references)
15. [Authors](#authors)

---

## About the Project

Digital Microfluidic Biochips (DMFBs) manipulate discrete droplets on a 2D electrode grid using electrowetting-on-dielectric (EWOD) to perform biochemical protocols such as PCR preparation, immunoassays, and clinical diagnostics. A critical step is **sample preparation** — combining reagents in a precise volumetric ratio by constructing a **mixing tree** (a binary tree of 50:50 mix-merge operations).

This project implements **five mixing-tree construction algorithms**, a **beam-search meta-framework**, a **323-test benchmark suite**, and an **interactive web-based visualization tool** for comparing algorithm performance across 7 metrics.

### Research Claim

> On a 323-test stratified DMFB sample-preparation benchmark, beam-augmented **LARP** achieves **17.6% fewer fluid splits** than RMA and **44.5% longer dilution chains** than BS, occupying a previously unfilled middle of the splits-vs-dilution trade-off curve.

---

## Key Features

### 🔬 Five Mixing-Tree Construction Algorithms
- **RMA** (Ratioed Mixing Algorithm) — Baseline, top-down greedy
- **BS** (Bit-Scanning) — Baseline, bottom-up bitwise
- **AP-DP** (Adaptive Partitioning via DP) — *Proposed*, multi-strategy
- **LARP** (Lookahead-Augmented Recursive Partitioning) — *Proposed*, with 2-level lookahead
- **ILP** (Integer-Allocation DP Heuristic) — *Proposed*, Pareto-optimal state tracking

### 🔍 Beam-Search + Memoization Meta-Framework
- Augments AP-DP, LARP, and ILP with beam-width `K=3`
- Achieves **9–14% improvement** on the design objective (splits and mixes)
- Memoization with canonical keys collapses repeated subproblems

### 📊 323-Test Stratified Benchmark
- 10 test categories (paper examples, edge cases, stress tests, random variants)
- 7 metrics per algorithm per test
- Deterministic PRNG seed for full reproducibility
- Baseline diff mode for rapid iteration

### 🖥️ Interactive Web Visualization
- **Smart Optimizer** — Weighted multi-metric recommendation engine
- **Five-Way Comparison** — Side-by-side metric comparison with bar charts
- **Per-Algorithm Pages** — Interactive tree visualization with zoom, pan, and color-coded fluid leaves
- Built with **Vite + React 19**

### 📄 Complete LaTeX Thesis Report
- 7-chapter BTP report following NITT thesis guidelines
- Comprehensive literature review of 10+ algorithm families
- All experimental data, tables, and analysis

---

## Project Structure

```
BTP_SEM8_MEDA/
├── src/                        # React visualization application
│   ├── App.jsx                 # Main app — all 5 algorithms + UI (1,237 lines)
│   ├── Optimizer.jsx           # Smart Optimizer page with weight sliders
│   ├── App.css                 # Application styles
│   ├── index.css               # Global styles
│   └── main.jsx                # React entry point
│
├── benchmark.js                # Algorithm library (all 5 algorithms + metrics)
├── bench.js                    # Benchmark runner (323 tests, stats, scoreboard)
│
├── index.html                  # Vite entry point
├── vite.config.js              # Vite configuration
├── package.json                # Node.js dependencies
│
├── BTPReport.tex               # Main LaTeX report (master file)
├── Introduction.tex            # Chapter 1: Problem domain & formal statement
├── LitReview.tex               # Chapter 2: Literature review (10+ algorithms)
├── Methodology.tex             # Chapter 3: Proposed algorithms & beam search
├── Benchmarking.tex            # Chapter 4: Framework & visualization tool
├── Results.tex                 # Chapter 5: Experimental results (323 tests)
├── Discussion.tex              # Chapter 6: Trade-offs & recommendations
├── Conclusion.tex              # Chapter 7: Conclusion & future work
├── abstract.tex                # Abstract
├── titlepage.tex               # Title page
├── bonafide.tex                # Certificate
├── declaration.tex             # Declaration
├── acknowledgement.tex         # Acknowledgement
├── appendix.tex                # Code attachments
├── references.bib              # BibTeX bibliography (10 references)
│
├── thesis.md                   # Comprehensive knowledge base (65KB)
├── benchmark_latest.json       # Latest benchmark results
├── benchmark_baseline.json     # Baseline for delta comparison
├── benchmark_detailed.csv      # Per-test detailed results
└── benchmark_results.csv       # Summary results
```

---

## Prerequisites

| Tool | Version | Purpose |
|------|---------|---------|
| **Node.js** | v18+ | Run benchmark & visualization tool |
| **npm** | v9+ | Package management |
| **LaTeX** (optional) | TeX Live / MiKTeX | Compile the thesis report |

---

## Installation & Setup

### 1. Clone the Repository

```bash
git clone <repository-url>
cd BTP_SEM8_MEDA
```

### 2. Install Dependencies

```bash
npm install
```

This installs:
- `react` and `react-dom` (v19) — UI framework
- `vite` (v7) — Build tool and dev server
- `@vitejs/plugin-react` — React support for Vite
- ESLint and related plugins — Code linting

---

## Running the Visualization Tool

Start the development server:

```bash
npm run dev
```

Open your browser at **http://localhost:5173**

### Using the Visualization Tool

1. **Enter fluid ratios** in the input field (comma-separated, e.g., `2,3,5,7,11,13,87`)
2. **Set the depth** parameter `d` (controls volume: `2^d` units)
3. **Click Generate** or switch between algorithm tabs

### Available Pages

| Tab | Description |
|-----|-------------|
| **Optimizer** | Smart recommendation engine — set importance weights for each metric and get the best algorithm |
| **Compare** | Side-by-side 5-way comparison with bar charts and winner highlights |
| **RMA** | Interactive RMA tree with metrics panel |
| **BS** | Interactive BS tree with metrics panel |
| **AP-DP** | Interactive AP-DP tree with metrics panel |
| **LARP** | Interactive LARP tree with metrics panel |
| **ILP** | Interactive ILP tree with metrics panel |

### Build for Production

```bash
npm run build
npm run preview
```

---

## Running the Benchmark Suite

Run the full 323-test benchmark:

```bash
node bench.js
```

### What Happens

1. Generates **323 deterministic test cases** across 10 categories
2. Runs all **6 algorithms** on every test case
3. Computes **7 metrics** per algorithm per test
4. Prints multi-table dashboards:
   - Overall averages per algorithm
   - Median values
   - Win counts per metric
   - Composite Pareto wins
   - Per-category breakdowns
   - Final scoreboard with points and medals
5. Saves output files:
   - `benchmark_latest.json` — Full JSON results
   - `benchmark_detailed.csv` — Per-test CSV with all metrics

### Baseline Comparison

On first run, a baseline is saved as `benchmark_baseline.json`. Subsequent runs will show delta indicators (`✓` for improvements, `✗` for regressions) against the baseline.

To reset the baseline:
```bash
cp benchmark_latest.json benchmark_baseline.json
```

### Sample Output

```
╔══════════════════════════════════════════════════╗
║        MIXING TREE BENCHMARK FRAMEWORK          ║
╚══════════════════════════════════════════════════╝

Tests: 323 total
  paper          9 cases
  edge           12 cases
  stress         2 cases
  rand-few       40 cases
  rand-med       80 cases
  ...

╔══════════════════════════════════════════════════
║  FINAL SCOREBOARD (7 metrics × [5,3,2,1,0] pts)
╠══════════════════════════════════════════════════
║  🥇 #1  BS     32 pts  ████████████████
║  🥈 #2  LARP   26 pts  █████████████
║  🥉 #3  ILP    24 pts  ████████████
║     #4  RMA    23 pts  ███████████
║     #5  AP-DP  21 pts  ██████████
╚══════════════════════════════════════════════════
```

---

## Algorithms Implemented

### Baselines (from literature)

| Algorithm | Direction | Strategy | Complexity |
|-----------|-----------|----------|------------|
| **RMA** | Top-down | Greedy dominant-fluid partitioning | O(d·n²) |
| **BS** | Bottom-up | Binary bit-scanning | O(d·n) |

### Proposed (this work)

| Algorithm | Strategy | Lookahead | Complexity |
|-----------|----------|-----------|------------|
| **AP-DP** | 5-strategy candidate generation + local scoring | None | O(2ⁿ + n·L) |
| **LARP** | Comprehensive subset-sum DP + lookahead scoring | 2-level | O(d·n·L) |
| **ILP** | Pareto-optimal integer-allocation DP | 1-level | O(d·n·L·CAP) |

### Recommended: **LARP with beam search K=3**
- Within ±10% of every algorithm's per-metric optimum
- Ranks last on only 5/323 tests (1.5%)
- 7ms average, 62ms worst case per construction

---

## Metrics Evaluated

| Metric | Symbol | Direction | Description |
|--------|--------|-----------|-------------|
| **Mixes** | `m` | ↓ Lower better | Number of 50:50 mix-merge operations |
| **Splits** | `splits` | ↓ Lower better | Fluids appearing in both subtrees |
| **Dilution** | `l` | ↑ Higher better | Total dilution-subtree depth |
| **Max Dilution** | `maxL` | ↑ Higher better | Longest single dilution chain |
| **Leaves** | `leaves` | ↓ Lower better | Number of dispense events |
| **Parallelism** | `p` | ↓ Lower better | Max concurrent mix operations |
| **Depth** | `d` | — | Tree depth (input-determined) |

---

## Visualization Tool Features

### 🎯 Smart Optimizer
- **7 weighted sliders** — Set importance for each metric (0–100)
- **3 presets** — Cost-optimized, Quality-optimized, Balanced
- **Auto-recommendation** — Ranks all 5 algorithms and shows the best
- **Winner tree display** — Interactive visualization of the recommended tree

### 📊 Five-Way Comparison
- **Bar chart view** — Visual comparison of all metrics across all algorithms
- **Summary table** — Per-metric winners with tie detection
- **Algorithm descriptions** — Brief info about each method

### 🌳 Per-Algorithm Tree Viewer
- **Canvas-based rendering** — Efficient rendering of large trees
- **Scroll to zoom** — Smooth zoom in/out
- **Drag to pan** — Navigate large trees
- **Color-coded leaves** — Each fluid gets a distinct color
- **Fluid legend** — Shows fluid names and approximated counts
- **7-metric dashboard** — All metrics displayed below the tree

### 🔄 Shared Input Panel
- Fluid ratios persist across all tabs
- Depth parameter shared globally
- Live info bar shows: `d=7 · 7 fluids · 2^7=128 units`

---

## Benchmark Framework Features

### Test Suite Categories (323 total)

| Category | Tests | Fluids | Depth | Purpose |
|----------|-------|--------|-------|---------|
| `paper` | 9 | as published | as published | Reproduce literature examples |
| `edge` | 12 | 2–4 | 3–5 | Small/edge-case ratios |
| `stress` | 2 | 30+ | 12 | Largest tractable inputs |
| `rand-few` | 40 | 2–4 | 4–8 | Random, small fluid count |
| `rand-med` | 80 | 5–8 | 6–9 | Random, mid-range count |
| `rand-many` | 60 | 10–18 | 7–10 | Random, many fluids |
| `rand-skew` | 40 | 5–8 | 6–9 | Strongly asymmetric counts |
| `rand-equal` | 40 | 5–8 | 6–9 | All counts within ±1 |
| `rand-deep` | 20 | 8–12 | 10–12 | Deep trees, high-volume |
| `rand-shallow` | 20 | 3–5 | 4–6 | Shallow trees |

### Output Tables

1. **Overall Averages** — Mean of each metric per algorithm
2. **Medians** — Robust central tendency
3. **Win Counts** — Per-test best counts (ties counted)
4. **Composite Pareto Wins** — Best all-around per test
5. **Per-Category Breakdown** — Averages for each test category
6. **Final Scoreboard** — Rank-based points with medals

### Reproducibility
- **Deterministic PRNG seed** (`42`) ensures identical test cases every run
- **JSON + CSV output** for external analysis
- **Baseline diff mode** for tracking improvements across iterations

---

## Key Results

| Comparison | Splits | Dilution | Mixes | Parallelism |
|------------|--------|----------|-------|-------------|
| **LARP vs RMA** | **−17.6% ✓** | −23.5% | **−12.9% ✓** | **−10.3% ✓** |
| **LARP vs BS** | +8.2% | **+44.5% ✓** | +5.7% | +2.2% |
| **Beam vs Greedy (LARP)** | **−14.0% ✓** | −20% | **−9.7% ✓** | −0.22 |

### Key Findings

1. **LARP and ILP fill the middle** of the splits-vs-dilution trade-off curve
2. **BS is unbeatable on splits** (100% tied-best) but catastrophic on dilution
3. **LARP and ILP produce identical trees** on 70% of tests (converged on practical optimum)
4. **Beam search yields 9–14% improvement** on the design objective at 10× wall-time

---

## References

1. S. Roy et al. — *Layout-Aware Solution Preparation for Biochemical Analysis on a DMFB* — VLSI Design, 2010
2. W. Thies et al. — *Abstraction Layers for Scalable Microfluidic Biocomputing* — DNA Computing, 2008
3. D. Mitra et al. — *Reactant Minimization using Skewed Mixing Trees* — ACM/IEEE, 2013
4. S. Bhattacharjee et al. — *Layout-Aware Mixture Preparation on DMFBs* — ACM TODAES, 2015
5. S. Bhattacharjee et al. — *Demand-Driven Mixture Preparation and Droplet Streaming* — DAC, 2014
6. S. Bhattacharjee et al. — *ILP-based Synthesis for Sample Preparation* — IEEE ICCD, 2016
7. A. Banerjee et al. — *SIMOP: Simulation-Guided Optimization for Sample Preparation* — Springer, 2021
8. T. Hsieh et al. — *REMIA — Reactant Minimization Algorithm* — IEEE TCAD, 2012
9. S. Roy et al. — *Waste-Aware Single-Target Dilution on DMFBs* — ScienceDirect, 2014
10. C.Y. Huang et al. — *Reactant Minimization with Various Mixing Models* — IEEE TCAD, 2015

---

## Authors

**Dhananjay Singh, Kalpit Sancheti, Krishankant Sharma**  
B.Tech Project — Semester 8  
Netaji Subhas University of Technology, Delhi

---

## License

This project is part of an academic B.Tech thesis. Please cite appropriately if using any algorithms or benchmark data from this work.
