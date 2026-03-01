import fs from 'fs';

const appContent = fs.readFileSync('src/App.jsx', 'utf-8');

const startIndex = appContent.indexOf('// ========== RATIO APPROXIMATION ==========');
const endIndex = appContent.indexOf('// ========== UI ==========');

if (startIndex === -1 || endIndex === -1) {
    console.error('Could not find boundaries in App.jsx');
    process.exit(1);
}

const algorithmsCode = appContent.substring(startIndex, endIndex);

let benchmarkContent = fs.readFileSync('benchmark.js', 'utf-8');

// Replace everything between RATIO APPROXIMATION and BENCHMARK SCRIPT
const benchStart = benchmarkContent.indexOf('// ========== RATIO APPROXIMATION ==========');
const benchEnd = benchmarkContent.indexOf('// ========== BENCHMARK SCRIPT ==========');

if (benchStart === -1 || benchEnd === -1) {
    console.error('Could not find boundaries in benchmark.js');
    process.exit(1);
}

const newBenchmarkContent =
    benchmarkContent.substring(0, benchStart) +
    algorithmsCode +
    benchmarkContent.substring(benchEnd);

// Now patch the BENCHMARK SCRIPT to include l and maxL
let patchedContent = newBenchmarkContent.replace(
    /const resultsMap = \{([\s\S]*?)\};/,
    `const resultsMap = {
    rma: { m: 0, d: 0, p: 0, leaves: 0, splits: 0, l: 0, maxL: 0 },
    bs: { m: 0, d: 0, p: 0, leaves: 0, splits: 0, l: 0, maxL: 0 },
    apdp: { m: 0, d: 0, p: 0, leaves: 0, splits: 0, l: 0, maxL: 0 },
    larp: { m: 0, d: 0, p: 0, leaves: 0, splits: 0, l: 0, maxL: 0 },
    ilp: { m: 0, d: 0, p: 0, leaves: 0, splits: 0, l: 0, maxL: 0 }
};`
);

patchedContent = patchedContent.replace(
    'csvDetails += `${i + 1},${numFluids},${rawRatiosStr},${algo},${stats.m},${stats.d},${stats.p},${stats.leaves},${stats.splits}\\n`;',
    'csvDetails += `${i + 1},${numFluids},${rawRatiosStr},${algo},${stats.m},${stats.d},${stats.p},${stats.leaves},${stats.splits},${stats.l},${stats.maxL}\\n`;'
);

patchedContent = patchedContent.replace(
    'let csvDetails = "TestId,NumFluids,Ratios,Algo,M_Mixes,D_Depth,P_Parallelism,Leaves,Splits\\n";',
    'let csvDetails = "TestId,NumFluids,Ratios,Algo,M_Mixes,D_Depth,P_Parallelism,Leaves,Splits,L_Dilution,MaxL_Longest\\n";'
);

patchedContent = patchedContent.replace(
    /resultsMap\[algo\]\.splits \+= stats\.splits;/,
    `resultsMap[algo].splits += stats.splits;
        resultsMap[algo].l += stats.l || 0;
        resultsMap[algo].maxL += stats.maxL || 0;`
);

patchedContent = patchedContent.replace(
    /splits: \(resultsMap\[algo\]\.splits \/ NUM_TESTS\)\.toFixed\(2\)/g,
    `splits: (resultsMap[algo].splits / NUM_TESTS).toFixed(2),
        l: (resultsMap[algo].l / NUM_TESTS).toFixed(2),
        maxL: (resultsMap[algo].maxL / NUM_TESTS).toFixed(2)`
);

patchedContent = patchedContent.replace(
    /\$\{formatRow\('splits', 'Splits'\)\}/,
    `\${formatRow('splits', 'Splits')}
\${formatRow('l', 'l(Sum)')}
\${formatRow('maxL', 'maxL')}`
);

fs.writeFileSync('benchmark.js', patchedContent);
console.log('Successfully synced benchmark.js with App.jsx improvements!');
