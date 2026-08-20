#!/usr/bin/env python3
"""Candidate DDI network — full analysis from the raw candidate-pair CSV.

Regenerates every candidate-network quantity reported in the manuscript
(Tables 1-3, S1, S3; concentration; power-law battery; connectivity;
centralities; Louvain; robustness) from Dataset_d3_candidate_ddi_pairs.csv.

Usage: python candidate_network_analysis.py <candidate_csv>
Output: out/candidate_network_results.json
Deterministic: Louvain seed 42; robustness RNG seed 42.
"""
import collections, json, sys, warnings
from pathlib import Path

import numpy as np
import pandas as pd
import networkx as nx

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parent
OUT = ROOT / "out"
OUT.mkdir(exist_ok=True)
SEED = 42

NAMES = {"DB00252": "Phenytoin", "DB00682": "Warfarin", "DB00489": "Sotalol",
         "DB00501": "Cimetidine", "DB00564": "Carbamazepine", "DB00199": "Erythromycin",
         "DB00814": "Meloxicam", "DB00945": "Acetylsalicylic acid",
         "DB01234": "Dexamethasone", "DB00717": "Norethisterone",
         "DB00363": "Clozapine", "DB00091": "Cyclosporine"}


def gini(values):
    v = np.sort(np.asarray(values, dtype=float))
    n = len(v)
    return float((2 * np.sum(np.arange(1, n + 1) * v)) / (n * v.sum()) - (n + 1) / n)


csv_path = Path(sys.argv[1])
df = pd.read_csv(csv_path)
id_cols = [c for c in df.columns if df[c].astype(str).str.match(r"^DB\d{5}$").mean() > 0.9]
if len(id_cols) < 2:
    sys.exit(f"could not find two DrugBank-ID columns in {list(df.columns)}")
a, b = id_cols[:2]

n_self = int((df[a] == df[b]).sum())
edges_sorted = sorted({tuple(sorted((x, y))) for x, y in zip(df[a], df[b]) if x != y})
G = nx.Graph()
G.add_edges_from(edges_sorted)
deg = dict(G.degree())
N, E = G.number_of_nodes(), G.number_of_edges()
degv = np.array(sorted(deg.values(), reverse=True))
order = sorted(deg.items(), key=lambda kv: -kv[1])
top10 = [k for k, _ in order[:10]]
tot_end = 2 * E

res = {"input_file": csv_path.name, "input_rows": int(len(df)),
       "self_loops_removed": n_self,
       "duplicate_pairs_collapsed": int(len(df) - n_self - E),
       "nodes": N, "edges": E,
       "mean_degree": round(float(degv.mean()), 2),
       "median_degree": float(np.median(degv)),
       "max_degree": int(degv.max()),
       "deg1_count": int((degv == 1).sum()),
       "deg1_pct": round(100 * (degv == 1).sum() / N, 1),
       "total_endpoints": tot_end,
       "gini": round(gini(degv), 3),
       "density": round(nx.density(G), 4),
       "avg_clustering": round(nx.average_clustering(G), 3)}

res["concentration"] = {}
for p in (0.01, 0.05, 0.10, 0.20):
    k = int(round(p * N))
    res["concentration"][f"{int(p*100)}%"] = [k, round(100 * degv[:k].sum() / tot_end, 1)]

# ---- connectivity
comps = sorted(nx.connected_components(G), key=len, reverse=True)
res["n_components"] = len(comps)
res["giant_size"] = len(comps[0])
res["giant_pct"] = round(100 * len(comps[0]) / N, 1)
res["component_sizes"] = [len(c) for c in comps]
res["isolated_pairs"] = [sorted(c) for c in comps if len(c) == 2]

# ---- power-law battery (pinned powerlaw package, automatic xmin)
import powerlaw
fit = powerlaw.Fit(degv, discrete=True, verbose=False)
res["powerlaw"] = {"alpha": round(float(fit.alpha), 2), "xmin": float(fit.xmin),
                   "tail_n": int((degv >= fit.xmin).sum())}
for alt in ("lognormal", "exponential"):
    R, p = fit.distribution_compare("power_law", alt, normalized_ratio=True)
    res["powerlaw"][f"vs_{alt}"] = {"R_normalized": round(float(R), 2), "p": float(p)}

# ---- threshold-sensitivity scan: does the model comparison depend on xmin?
scan = {}
for xm in (10, 20, 32, 50, 81):
    tail_n = int((degv >= xm).sum())
    if tail_n < 50:
        continue
    f2 = powerlaw.Fit(degv, discrete=True, xmin=xm, verbose=False)
    row = {"alpha": round(float(f2.power_law.alpha), 2), "tail_n": tail_n}
    for alt in ("lognormal", "exponential"):
        R, p = f2.distribution_compare("power_law", alt, normalized_ratio=True)
        row[f"vs_{alt}"] = {"R_normalized": round(float(R), 2), "p": float(p)}
    scan[f"xmin={xm}"] = row
res["powerlaw_threshold_sensitivity"] = scan

# ---- hubs
res["top12_degree"] = [[d, int(deg[d]), NAMES.get(d, "")] for d, _ in order[:12]]

# ---- centralities (all computed exactly on the giant component, matching Table S1)
Gc0 = G.subgraph(sorted(comps[0])).copy()
bet = nx.betweenness_centrality(Gc0, normalized=True)
bet_order = sorted(bet.items(), key=lambda kv: -kv[1])
deg_rank = {d: r + 1 for r, (d, _) in enumerate(order)}
res["top10_betweenness"] = [[d, round(v, 4), deg_rank[d], NAMES.get(d, "")]
                            for d, v in bet_order[:10]]
res["betweenness_overlap_with_degree_top10"] = len(set(top10) & {d for d, _ in bet_order[:10]})

Gc = G.subgraph(comps[0]).copy()
eig = nx.eigenvector_centrality(Gc, max_iter=1000)
clo = nx.closeness_centrality(Gc)
for nm, cen in (("eigenvector", eig), ("closeness", clo)):
    top = [d for d, _ in sorted(cen.items(), key=lambda kv: -kv[1])[:10]]
    res[f"overlap_{nm}"] = [len(set(top10) & set(top)),
                            [[d, deg_rank[d], NAMES.get(d, "")] for d in top if d not in top10]]

# ---- Louvain (python-louvain, seed 42)
import community as community_louvain
part = community_louvain.best_partition(Gc, random_state=SEED)
Q = community_louvain.modularity(part, Gc)
cc = collections.Counter(part.values())
ordered_comm = cc.most_common()
size_rank = {c: i + 1 for i, (c, _) in enumerate(ordered_comm)}
res["louvain"] = {"modularity": round(Q, 3),
                  "n_comm_giant": len(ordered_comm),
                  "sizes_giant": [s for _, s in ordered_comm],
                  "n_comm_full": len(ordered_comm) + sum(1 for c in comps if len(c) == 2),
                  "hub_communities": {NAMES.get(h, h): size_rank[part[h]] for h in top10},
                  "hub_n_distinct_communities": len({size_rank[part[h]] for h in top10})}

# ---- robustness: random edge removal
rng = np.random.default_rng(SEED)
edges_list = list(G.edges())


def top10_after_removal(frac, rng):
    keep = rng.choice(len(edges_list), size=int(round((1 - frac) * len(edges_list))),
                      replace=False)
    H = nx.Graph()
    H.add_nodes_from(G.nodes())
    H.add_edges_from(edges_list[i] for i in keep)
    d2 = dict(H.degree())
    return [k for k, _ in sorted(d2.items(), key=lambda kv: -kv[1])[:10]]


rob = {}
for frac, ntrials in ((0.10, 20), (0.20, 20)):
    jac = []
    for _ in range(ntrials):
        t = set(top10_after_removal(frac, rng))
        jac.append(len(t & set(top10)) / len(t | set(top10)))
    rob[f"jaccard_{int(frac*100)}pct"] = {"mean": round(float(np.mean(jac)), 3),
                                          "min": round(float(np.min(jac)), 3),
                                          "n_trials": ntrials}
ret = collections.Counter()
NT = 50
for _ in range(NT):
    t = set(top10_after_removal(0.20, rng))
    for h in top10:
        if h in t:
            ret[h] += 1
rob["per_drug_retention_20pct"] = {NAMES.get(h, h): round(100 * ret[h] / NT)
                                   for h in top10}
res["robustness"] = rob

(OUT / "candidate_network_results.json").write_text(json.dumps(res, indent=2))
print(json.dumps(res, indent=2))
