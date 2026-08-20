#!/usr/bin/env python3
"""Confirmed-interaction (GoldD3R) network analysis — regenerates every
confirmed-network value reported in the manuscript from the raw file.

Input : GoldD3R.txt  (lines: CUI_a,CUI_b[,\[MechanismLabel[/Label2...]\]])
Output: confirmed_network_results.json
"""
import collections, json, re, sys
from pathlib import Path

import networkx as nx
import numpy as np

ROOT = Path(__file__).resolve().parent
GOLD = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "GoldD3R.txt"
OUT = ROOT / "out"
OUT.mkdir(exist_ok=True)


def gini(values):
    v = np.sort(np.asarray(values, dtype=float))
    n = len(v)
    return float((2 * np.sum(np.arange(1, n + 1) * v)) / (n * v.sum()) - (n + 1) / n)


edges = set()
mech_pairs = 0
mech_counts = collections.Counter()
raw_lines = 0
for ln in GOLD.read_text().strip().split("\n"):
    ln = ln.strip()
    if not ln:
        continue
    raw_lines += 1
    m = re.match(r"^(C\d+),(C\d+)(?:,\[([^\]]*)\])?", ln)
    if not m:
        continue
    a, b = m.group(1), m.group(2)
    if a == b:
        continue
    edges.add((min(a, b), max(a, b)))
    if m.group(3):
        mech_pairs += 1
        for lab in re.split(r"[,/;|]", m.group(3)):
            lab = lab.strip()
            if lab:
                mech_counts[lab] += 1

G = nx.Graph()
G.add_edges_from(edges)
deg = dict(G.degree())
degv = np.array(sorted(deg.values(), reverse=True))
N, E = G.number_of_nodes(), G.number_of_edges()
tot_endpoints = 2 * E

conc = {}
for p in (0.01, 0.05, 0.10, 0.20):
    k = int(round(p * N))
    conc[f"{int(p*100)}%"] = [k, round(100 * degv[:k].sum() / tot_endpoints, 1)]

comps = sorted(nx.connected_components(G), key=len, reverse=True)

res = {
    "raw_lines": raw_lines,
    "nodes": N,
    "edges": E,
    "mean_degree": round(float(np.mean(degv)), 2),
    "median_degree": float(np.median(degv)),
    "max_degree": int(degv.max()),
    "deg1_count": int((degv == 1).sum()),
    "deg1_pct": round(100 * (degv == 1).sum() / N, 1),
    "concentration": conc,
    "gini": round(gini(degv), 3),
    "density": round(nx.density(G), 4),
    "avg_clustering": round(nx.average_clustering(G), 3),
    "n_components": len(comps),
    "giant_size": len(comps[0]),
    "giant_pct": round(100 * len(comps[0]) / N, 1),
    "component_sizes_head": [len(c) for c in comps[:8]],
    "mechanism_labelled_pairs": mech_pairs,
    "mechanism_label_occurrences": int(sum(mech_counts.values())),
    "mechanism_labels": dict(mech_counts.most_common()),
}

try:
    import powerlaw
    fit = powerlaw.Fit(degv, discrete=True, verbose=False)
    res["powerlaw_alpha"] = round(float(fit.alpha), 2)
    res["powerlaw_xmin"] = float(fit.xmin)
    res["powerlaw_tail_n"] = int((degv >= fit.xmin).sum())
    for alt in ("lognormal", "exponential"):
        R, p = fit.distribution_compare("power_law", alt, normalized_ratio=True)
        res[f"pl_vs_{alt}_R_normalized"] = round(float(R), 2)
        res[f"pl_vs_{alt}_p"] = float(p)
except Exception as e:  # pragma: no cover
    res["powerlaw_error"] = str(e)

(OUT / "confirmed_network_results.json").write_text(json.dumps(res, indent=2))
print(json.dumps(res, indent=2))
