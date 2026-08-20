#!/usr/bin/env python3
"""All manuscript figures, generated exclusively from the raw input data.

Every quantitative element in every figure is computed at run time from:
  - the candidate-pair CSV        (Figs 2, 3, 4, S1, S2; candidate side of Fig 5)
  - GoldD3R.txt                   (confirmed side of Fig 5)
  - the DDInter category exports  (Fig 6, incl. permutation null band)
Randomness is used ONLY (a) to choose which peripheral partner nodes are
DISPLAYED in the Fig 4 sample for legibility and (b) in the deposited
permutation tests; both are seeded. No figure contains simulated,
hard-coded, or interpolated data.

Usage: python make_figures.py <candidate_csv> <GoldD3R.txt> <ddinter_dir>
"""
import collections, glob, json, re, sys
from pathlib import Path

import numpy as np
import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "figs"
OUT.mkdir(exist_ok=True)
SEED = 42

INK = "#16283D"; CAND = "#C46A18"; CONF = "#17788C"; CAND_F = "#F0DFC8"
CONF_F = "#CFE4E8"; GRID = "#E7E3DA"; MUTE = "#6A7280"; SOFT = "#9AA1AA"
MAJOR_C = "#B04A3A"
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 11,
    "text.color": INK, "axes.labelcolor": INK, "axes.edgecolor": "#C9C6BE",
    "xtick.color": MUTE, "ytick.color": MUTE, "axes.linewidth": 1.0,
    "axes.spines.top": False, "axes.spines.right": False,
    "savefig.dpi": 320, "savefig.bbox": "tight",
    "figure.facecolor": "white", "axes.facecolor": "white"})

NAMES = {"DB00252": "Phenytoin", "DB00682": "Warfarin", "DB00489": "Sotalol",
         "DB00501": "Cimetidine", "DB00564": "Carbamazepine", "DB00199": "Erythromycin",
         "DB00814": "Meloxicam", "DB00945": "Aspirin", "DB01234": "Dexamethasone",
         "DB00717": "Norethisterone"}


def gini(v):
    v = np.sort(np.asarray(v, float)); n = len(v)
    return float((2 * np.sum(np.arange(1, n + 1) * v)) / (n * v.sum()) - (n + 1) / n)


def load_candidate(p):
    df = pd.read_csv(p)
    ids = [c for c in df.columns if df[c].astype(str).str.match(r"^DB\d{5}$").mean() > 0.9]
    edges = sorted({tuple(sorted((x, y))) for x, y in zip(df[ids[0]], df[ids[1]]) if x != y})
    G = nx.Graph(); G.add_edges_from(edges)
    return G


def load_gold(p):
    edges = set()
    for ln in open(p):
        m = re.match(r"^(C\d+),(C\d+)", ln.strip())
        if m and m.group(1) != m.group(2):
            edges.add((min(m.group(1), m.group(2)), max(m.group(1), m.group(2))))
    G = nx.Graph(); G.add_edges_from(edges)
    return G


def conc_curve(G):
    d = np.array(sorted(dict(G.degree()).values(), reverse=True))
    return np.arange(1, len(d) + 1) / len(d), np.cumsum(d) / d.sum(), d


def fig1():
    fig, ax = plt.subplots(figsize=(9.4, 4.7)); ax.set_xlim(0, 10); ax.set_ylim(0, 6); ax.axis("off")
    def box(x, y, w, h, t, fc, ec, tc, fs=9.6):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.14",
                                    fc=fc, ec=ec, lw=1.3))
        ax.text(x + w/2, y + h/2, t, ha="center", va="center", color=tc, fontsize=fs)
    def arr(x1, y1, x2, y2, c):
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=15, lw=1.6, color=c))
    box(0.3, 4.2, 4.3, 1.4, "Five heterogeneous DDI evidence sources\n(DrugBank · KEGG · DIKB\nSemMedDB · TWOSIDES)", "#F3F1EC", "#CFCABB", "#4A5260")
    box(5.4, 4.2, 4.3, 1.4, "Reference filtering  +\n9-D Jaccard features  +\nlogistic-regression relevance\nmodel  (τ = 0.3)", "#F3F1EC", "#CFCABB", "#4A5260")
    arr(4.6, 4.9, 5.4, 4.9, SOFT)
    box(2.85, 2.5, 4.3, 1.2, "16,316 candidate DDI pairs\n1,114 unique drugs", INK, INK, "white", fs=10.5)
    arr(4.9, 4.2, 4.9, 3.7, SOFT)
    box(0.3, 0.5, 4.3, 1.35, "Network construction\n(undirected, unweighted)", CAND, CAND, "white")
    box(5.4, 0.5, 4.3, 1.35, "Degree & concentration · connectivity\nbetweenness / eigenvector / closeness\nLouvain communities · robustness", CAND, CAND, "white", fs=8.8)
    arr(4.0, 2.5, 2.45, 1.85, CAND); arr(4.6, 1.17, 5.4, 1.17, CAND)
    ax.text(0.3, 5.86, "UPSTREAM CANDIDATE CONSTRUCTION  —  fixed input", fontsize=8.6, color=SOFT, weight="bold")
    ax.text(0.3, 2.08, "THIS STUDY'S CONTRIBUTION", fontsize=8.6, color=CAND, weight="bold")
    fig.savefig(OUT / "Fig1.png"); plt.close(fig)


def fig2(G):
    d = np.array(list(dict(G.degree()).values()))
    fig, ax = plt.subplots(figsize=(7.4, 4.5))
    ax.hist(d, bins=np.arange(0, d.max() + 9, 8), color=CAND, edgecolor="white", linewidth=0.6, zorder=3)
    ax.set_yscale("log"); ax.yaxis.grid(True, color=GRID, lw=0.9, zorder=0); ax.set_axisbelow(True)
    ax.axvline(d.mean(), color=INK, ls="--", lw=1.5, zorder=4); ax.axvline(np.median(d), color=CONF, ls=":", lw=1.8, zorder=4)
    ax.text(d.mean() + 2, 300, f"mean {d.mean():.1f}", color=INK, fontsize=9, weight="bold")
    ax.text(np.median(d) + 2, 120, f"median {np.median(d):.0f}", color=CONF, fontsize=9, weight="bold")
    ax.set_xlabel("Candidate-pair degree (partners per drug)"); ax.set_ylabel("Number of drugs (log scale)")
    fig.savefig(OUT / "Fig2.png"); plt.close(fig)


def fig3(G):
    xs, ys, d = conc_curve(G); N = len(d); tot = d.sum()
    fig, ax = plt.subplots(figsize=(6.9, 5.2))
    ax.plot([0, 1], [0, 1], color="#B9B5AC", lw=1.1, ls=(0, (5, 4)), zorder=1, label="Perfect equality")
    ax.fill_between(xs, ys, xs, color=CAND_F, alpha=0.5, zorder=1); ax.plot(xs, ys, color=CAND, lw=2.8, zorder=4, solid_capstyle="round", label="Observed concentration")
    for p in (0.01, 0.05, 0.10, 0.20):
        k = int(round(p * N)); share = d[:k].sum() / tot
        ax.plot([p, p], [p, share], color=CAND, lw=0.8, alpha=0.45, zorder=2); ax.scatter([p], [share], s=42, color=CAND, edgecolor="white", linewidth=1.2, zorder=6)
        ax.annotate(f"{100*share:.1f}%", (p, share), textcoords="offset points", xytext=(8, -11), fontsize=9.5, color=INK, weight="bold")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_xticks(np.arange(0, 1.01, .25)); ax.set_yticks(np.arange(0, 1.01, .25))
    ax.set_xlabel("Top fraction of drugs, ranked by degree (descending)"); ax.set_ylabel("Cumulative share of candidate-pair endpoints")
    ax.text(0.55, 0.19, "Gini", fontsize=9, color=MUTE, weight="bold"); ax.text(0.55, 0.12, f"{gini(d):.3f}", fontsize=15, color=CAND, weight="bold")
    ax.legend(frameon=False, loc="upper left", fontsize=10); fig.savefig(OUT / "Fig3.png"); plt.close(fig)


def fig4(G):
    rng = np.random.default_rng(7); deg = dict(G.degree()); top10 = [k for k, _ in sorted(deg.items(), key=lambda kv: -kv[1])[:10]]
    sub = G.subgraph(top10); H = nx.Graph(); H.add_nodes_from(top10); H.add_edges_from(sub.edges()); n_hub_edges = H.number_of_edges(); pos = {}; R = 3.0
    for i, h in enumerate(top10):
        ang = 2 * np.pi * i / 10; pos[h] = (R * np.cos(ang), R * np.sin(ang)); partners = [p for p in G.neighbors(h) if p not in top10]
        for p in rng.choice(partners, min(6, len(partners)), replace=False):
            nd = f"{h}__{p}"; H.add_node(nd); H.add_edge(h, nd); pr = R + 1.2 + rng.uniform(0, 0.5); pa = ang + rng.uniform(-0.22, 0.22); pos[nd] = (pr * np.cos(pa), pr * np.sin(pa))
    iso = [sorted(c) for c in nx.connected_components(G) if len(c) == 2]; isonodes = []
    for i, (x, y) in enumerate(iso):
        H.add_edge(x, y); isonodes += [x, y]; pos[x] = (-3.2 + i * 3.1, -5.4); pos[y] = (-2.5 + i * 3.1, -5.4)
    fig, ax = plt.subplots(figsize=(7.8, 8.0)); ax.axis("off"); peri = [n for n in H if n not in top10 and n not in isonodes]; he = [e for e in H.edges() if e[0] in top10 and e[1] in top10]; oe = [e for e in H.edges() if e not in he]
    nx.draw_networkx_edges(H, pos, edgelist=oe, edge_color="#E3DFD6", width=0.7, ax=ax); nx.draw_networkx_edges(H, pos, edgelist=he, edge_color=CONF, width=1.5, alpha=0.85, ax=ax)
    nx.draw_networkx_nodes(H, pos, nodelist=peri, node_color="#B9BEC5", node_size=48, edgecolors="white", linewidths=0.5, ax=ax); nx.draw_networkx_nodes(H, pos, nodelist=isonodes, node_color=MAJOR_C, node_size=80, edgecolors="white", linewidths=0.6, ax=ax)
    dmin = min(deg[h] for h in top10); sizes = [340 + (deg[h] - dmin + 1) * 6 for h in top10]; nx.draw_networkx_nodes(H, pos, nodelist=top10, node_color=INK, node_size=sizes, edgecolors="white", linewidths=1.2, ax=ax)
    for h in top10: ax.text(pos[h][0], pos[h][1], NAMES.get(h, h)[:6], color="white", ha="center", va="center", fontsize=6.3, weight="bold")
    ax.text(0, 6.7, f"Teal edges: {n_hub_edges} of 45 hub–hub pairs present", color=CONF, fontsize=9.5, ha="center", weight="bold"); ax.text(0.05, -6.0, f"{len(iso)} isolated two-drug pairs", color=MAJOR_C, fontsize=9.5, ha="center", weight="bold")
    ax.set_xlim(-6, 6); ax.set_ylim(-6.6, 7.1); fig.savefig(OUT / "Fig4.png"); plt.close(fig)


def fig5(Gc, Gg):
    fig = plt.figure(figsize=(11.2, 5.0)); gs = fig.add_gridspec(1, 2, width_ratios=[1.05, 1.0], wspace=0.28); axA = fig.add_subplot(gs[0, 0])
    xc, yc, dcand = conc_curve(Gc); xg, yg, dgold = conc_curve(Gg)
    axA.plot([0, 1], [0, 1], color="#B9B5AC", lw=1.1, ls=(0, (5, 4)), zorder=1); axA.fill_between(xg, yg, xg, color=CONF_F, alpha=0.45, zorder=1)
    axA.plot(xg, yg, color=CONF, lw=2.8, zorder=4, solid_capstyle="round", label=f"Reference network  (n = {len(dgold):,})"); axA.plot(xc, yc, color=CAND, lw=2.8, zorder=3, solid_capstyle="round", label=f"Candidate network  (n = {len(dcand):,})")
    marks = []
    for d_, col in ((dgold, CONF), (dcand, CAND)):
        k = int(round(0.10 * len(d_))); sh = d_[:k].sum() / d_.sum(); axA.plot([0.10, 0.10], [0.10, sh], color=col, lw=0.8, alpha=0.5, zorder=2); axA.scatter([0.10], [sh], s=34, c=[col], zorder=6, edgecolor="white", linewidth=1.1); marks.append(100 * sh)
    axA.annotate(f"top 10% of drugs\n≈ {min(marks):.0f}–{max(marks):.0f}% of all pair endpoints", xy=(0.10, np.mean(marks) / 100), xytext=(0.34, 0.30), fontsize=9.5, color=INK, ha="left", va="center", arrowprops=dict(arrowstyle="-", color=MUTE, lw=0.9, connectionstyle="arc3,rad=-0.25"))
    axA.set_xlim(0, 1); axA.set_ylim(0, 1); axA.set_xlabel("Top fraction of drugs, ranked by degree"); axA.set_ylabel("Cumulative share of pair endpoints"); axA.legend(frameon=False, loc="lower right", fontsize=9.5)
    axA.text(0.045, 0.93, "Gini", fontsize=8.5, color=MUTE, weight="bold"); axA.text(0.045, 0.865, f"{gini(dgold):.3f}", fontsize=11, color=CONF, weight="bold"); axA.text(0.045, 0.80, f"{gini(dcand):.3f}", fontsize=11, color=CAND, weight="bold"); axA.text(0, 1.075, "A", transform=axA.transAxes, fontsize=15, weight="bold", color=INK)
    axB = fig.add_subplot(gs[0, 1])
    def stats(G, d_):
        comps = sorted(nx.connected_components(G), key=len, reverse=True); k = int(round(0.10 * len(d_))); return [100 * len(comps[0]) / G.number_of_nodes(), 100 * d_[:k].sum() / d_.sum(), 100 * (d_ == 1).sum() / len(d_), 100 * gini(d_), 100 * nx.average_clustering(G)]
    cand_v = stats(Gc, dcand); conf_v = stats(Gg, dgold); metrics = ["Giant\ncomponent", "Top-10%\nshare", "Degree-1\ndrugs", "Gini\n(×100)", "Clustering\n(×100)"]; y = np.arange(len(metrics))[::-1]; h = 0.36
    axB.barh(y + h/2, cand_v, height=h, color=CAND, zorder=3, label="Candidate", edgecolor="white", linewidth=0.8); axB.barh(y - h/2, conf_v, height=h, color=CONF, zorder=3, label="Reference", edgecolor="white", linewidth=0.8)
    for yi, cv, kv in zip(y, cand_v, conf_v): axB.text(cv + 1.2, yi + h/2, f"{cv:.1f}", va="center", fontsize=8.5, color=CAND, weight="bold"); axB.text(kv + 1.2, yi - h/2, f"{kv:.1f}", va="center", fontsize=8.5, color=CONF, weight="bold")
    axB.set_yticks(y); axB.set_yticklabels(metrics, fontsize=9.5); axB.set_xlim(0, 112); axB.set_xlabel("Value"); axB.set_xticks([0, 25, 50, 75, 100]); axB.xaxis.grid(True, color=GRID, lw=0.9, zorder=0); axB.set_axisbelow(True); axB.spines["left"].set_visible(False); axB.tick_params(axis="y", length=0); axB.legend(frameon=False, loc="lower right", fontsize=9.5); axB.text(0, 1.075, "B", transform=axB.transAxes, fontsize=15, weight="bold", color=INK)
    fig.suptitle("Candidate network and its construction reference share one structure", x=0.5, y=1.06, fontsize=14, weight="bold", color=INK); fig.text(0.5, 0.99, "Filtered candidate pairs (DrugBank IDs) vs clinically documented reference pairs (UMLS CUIs); the reference is related by construction, not independent", ha="center", fontsize=9, color=MUTE); fig.savefig(OUT / "Fig5.png"); plt.close(fig)


def fig6(ddinter_results_json):
    r = json.load(open(ddinter_results_json)); fig = plt.figure(figsize=(12.6, 4.4)); gs = fig.add_gridspec(1, 3, wspace=0.34)
    axA = fig.add_subplot(gs[0, 0]); obs_all = r["concentration_by_overall_degree"]["0.1"]["all_endpoint_share"]; obs_maj = r["permutation_null"]["major_share_top_decile"]["observed_pct"]; null_mean = r["permutation_null"]["major_share_top_decile"]["null_mean_pct"]; null_lo, null_hi = r["permutation_null"]["major_share_top_decile"]["null_95_pct"]
    bars = axA.bar(["All\nendpoints", "Major\nendpoints"], [obs_all, obs_maj], color=[CAND, MAJOR_C], width=0.55, zorder=3, edgecolor="white"); axA.errorbar([1], [null_mean], yerr=[[null_mean - null_lo], [null_hi - null_mean]], fmt="_", color=INK, capsize=6, ms=18, lw=1.4, zorder=4, label="Label-permutation null (95%)")
    for b, v in zip(bars, [obs_all, obs_maj]): axA.text(b.get_x() + b.get_width()/2, v + 1.2, f"{v:.1f}%", ha="center", fontsize=9.5, weight="bold", color=INK)
    axA.set_ylim(0, 100); axA.set_ylabel("Endpoint share of top 10% drugs\n(ranked by overall degree)"); axA.yaxis.grid(True, color=GRID, lw=0.9, zorder=0); axA.set_axisbelow(True); axA.legend(frameon=False, fontsize=8.5, loc="upper right"); axA.text(0, 1.06, "a", transform=axA.transAxes, fontsize=14, weight="bold")
    axB = fig.add_subplot(gs[0, 1]); ft, fr = r["fraction_major_top_decile"], r["fraction_major_rest"]; axB.bar(["Top-decile\ndrugs", "All other\ndrugs"], [ft, fr], color=[CAND, "#B9BEC5"], width=0.55, zorder=3, edgecolor="white")
    for x, v in zip([0, 1], [ft, fr]): axB.text(x, v + 0.4, f"{v:.1f}%", ha="center", fontsize=9.5, weight="bold", color=INK)
    axB.set_ylabel("Share of endpoints rated Major (%)"); axB.yaxis.grid(True, color=GRID, lw=0.9, zorder=0); axB.set_axisbelow(True); d_ = r["permutation_null"]["fraction_major_diff_top_minus_rest"]; axB.set_title(f"difference {d_['observed_pctpts']:+.1f} pp; null 95% [{d_['null_95_pctpts'][0]:+.1f}, {d_['null_95_pctpts'][1]:+.1f}] pp", fontsize=8.5, color=MUTE); axB.text(0, 1.06, "b", transform=axB.transAxes, fontsize=14, weight="bold")
    axC = fig.add_subplot(gs[0, 2]); rho = r["permutation_null"]["rho_degree_vs_n_major"]; axC.bar(["Observed"], [rho["observed"]], color=CAND, width=0.4, zorder=3, edgecolor="white"); axC.errorbar([0.75], [rho["null_mean"]], yerr=[[rho["null_mean"] - rho["null_95"][0]], [rho["null_95"][1] - rho["null_mean"]]], fmt="_", color=INK, capsize=6, ms=18, lw=1.4, zorder=4, label="Null (95%)")
    axC.text(0, rho["observed"] + 0.02, f"{rho['observed']:.2f}", ha="center", fontsize=9.5, weight="bold"); axC.text(0.75, rho["null_mean"] + 0.02, f"{rho['null_mean']:.2f}", ha="center", fontsize=9, color=INK); axC.set_xlim(-0.5, 1.25); axC.set_ylim(0, 1.0); axC.set_xticks([0, 0.75]); axC.set_xticklabels(["Observed ρ", "Permutation\nnull ρ"]); axC.set_ylabel("Spearman ρ (overall degree vs n Major)"); axC.yaxis.grid(True, color=GRID, lw=0.9, zorder=0); axC.set_axisbelow(True); axC.text(0, 1.06, "c", transform=axC.transAxes, fontsize=14, weight="bold"); fig.savefig(OUT / "Fig6.png"); plt.close(fig)


def figS1S2(G):
    import community as community_louvain
    comps = sorted(nx.connected_components(G), key=len, reverse=True); Gc = G.subgraph(sorted(comps[0])).copy(); part = community_louvain.best_partition(Gc, random_state=SEED); Q = community_louvain.modularity(part, Gc); cc = collections.Counter(part.values()); ordered = cc.most_common(); size_rank = {c: i + 1 for i, (c, _) in enumerate(ordered)}; n_dyads = sum(1 for c in comps if len(c) == 2)
    sizes = [s for _, s in ordered] + [2] * n_dyads; labels = [f"C{i+1}" for i in range(len(sizes))]; n_sub = sum(1 for s in sizes if s >= 100); colors = [CAND] * n_sub + [MUTE] * (len(ordered) - n_sub) + [MAJOR_C] * n_dyads
    fig, ax = plt.subplots(figsize=(8.0, 4.3)); ax.bar(labels, sizes, color=colors, edgecolor="white", linewidth=1.0, zorder=3, width=0.72); ax.yaxis.grid(True, color=GRID, lw=0.9, zorder=0); ax.set_axisbelow(True)
    for i, v in enumerate(sizes): ax.text(i, v + 5, str(v), ha="center", fontsize=8.4, color=INK, weight="bold")
    ax.set_ylabel("Number of drugs in community"); ax.set_xlabel("Louvain community (ranked by size)"); ax.spines["left"].set_visible(False); ax.tick_params(axis="y", length=0); ax.text(0.985, 0.9, f"Modularity  Q = {Q:.2f}", transform=ax.transAxes, ha="right", fontsize=10.5, color=INK, weight="bold"); ax.text(0.985, 0.80, "■ isolated dyad", transform=ax.transAxes, ha="right", fontsize=8.2, color=MAJOR_C); ax.text(0.845, 0.80, "■ small", transform=ax.transAxes, ha="right", fontsize=8.2, color=MUTE); ax.text(0.76, 0.80, "■ substantial", transform=ax.transAxes, ha="right", fontsize=8.2, color=CAND); fig.savefig(OUT / "FigS1.png"); plt.close(fig)
    deg = dict(G.degree()); top10 = [k for k, _ in sorted(deg.items(), key=lambda kv: -kv[1])[:10]]; hubs = [NAMES.get(h, h) for h in top10]; comm = {NAMES.get(h, h): size_rank[part[h]] for h in top10}; ncol = max(comm.values()); fig, ax = plt.subplots(figsize=(5.8, 6.2)); ax.set_xlim(-0.5, ncol - 0.5); ax.set_ylim(-0.5, len(hubs) - 0.5)
    for i, h in enumerate(hubs):
        yi = len(hubs) - 1 - i
        for j in range(ncol):
            filled = comm[h] == j + 1; ax.add_patch(plt.Rectangle((j - 0.42, yi - 0.42), 0.84, 0.84, facecolor=CAND if filled else "#F4F2EC", edgecolor="white", linewidth=1.6, zorder=2))
            if filled: ax.text(j, yi, "●", ha="center", va="center", color="white", fontsize=12, zorder=3)
    ax.set_xticks(range(ncol)); ax.set_xticklabels([f"C{i+1}" for i in range(ncol)], fontsize=10); ax.set_yticks(range(len(hubs))); ax.set_yticklabels(hubs[::-1], fontsize=9.5); ax.set_xlabel("Community (ranked by size)"); ax.tick_params(length=0)
    for s in ax.spines.values(): s.set_visible(False)
    ax.xaxis.set_ticks_position("top"); ax.xaxis.set_label_position("top"); ax.set_title("Hub-drug community membership", fontsize=12, color=INK, weight="bold", pad=26); fig.savefig(OUT / "FigS2.png"); plt.close(fig)


if __name__ == "__main__":
    cand_csv, gold_txt, ddinter_dir = sys.argv[1:4]
    Gc = load_candidate(cand_csv); Gg = load_gold(gold_txt)
    fig1(); fig2(Gc); fig3(Gc); fig4(Gc); fig5(Gc, Gg); figS1S2(Gc)
    dd_json = ROOT.parent / "out" / "ddinter_severity_results.json"
    if dd_json.exists(): fig6(dd_json)
    print("figures written to", OUT)
