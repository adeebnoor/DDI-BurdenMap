#!/usr/bin/env python3
"""
ddinter_severity_analysis.py
Non-circular severity analysis of the complete DDInter export for
"Candidate drug-drug interactions concentrate on hub drugs".

Reproduces, from the raw export alone, every DDInter number in the manuscript
AND regenerates Figure 6 AND the candidate hub-drug rank lookup:
  - concentration (top-k endpoint share) with a random-drug-set null
  - top-decile Major-endpoint capture with a label-permutation null
  - fraction-Major (top decile vs rest) with a label-permutation null
  - Spearman(overall degree, per-drug Major count) vs part-whole null
  - DDInter apex drugs (top 10 by degree)
  - location of the candidate network's 10 named hub drugs in DDInter
  - Figure 6 (3 panels) written to Fig6_ddinter_severity_corrected.png

Ranking is ALWAYS by overall degree (never Major-degree): the design that
avoids the circularity identified in review. All randomness is seeded (42).

Usage:  python ddinter_severity_analysis.py ddinter2_unique.csv out.json fig6.png
Input columns: DDInterID_A, Drug_A, DDInterID_B, Drug_B, Level[, atc_cat]
  Level in {Major, Moderate, Minor, Unknown}; export covers ATC A,B,D,H,L,P,R,V.
"""
import sys, json
import numpy as np, pandas as pd, networkx as nx
from scipy.stats import spearmanr

SEED = 42
B = 10000

PAPER_HUBS = ["Phenytoin", "Warfarin", "Sotalol", "Cimetidine", "Carbamazepine",
              "Erythromycin", "Meloxicam", "Acetylsalicylic acid", "Dexamethasone",
              "Norethisterone"]
SYNONYMS = {"Acetylsalicylic acid": ["acetylsalicylic acid", "aspirin"],
            "Norethisterone": ["norethisterone", "norethindrone"]}

def gini(x):
    x = np.sort(np.asarray(x, float)); n = len(x); c = np.cumsum(x)
    return (n + 1 - 2 * np.sum(c) / c[-1]) / n

def analyse(inp):
    raw = pd.read_csv(inp)
    raw = raw[raw.DDInterID_A != raw.DDInterID_B].copy()
    order = {"Major": 3, "Moderate": 2, "Minor": 1, "Unknown": 0}
    raw["key"] = raw.apply(lambda r: tuple(sorted((r.DDInterID_A, r.DDInterID_B))), axis=1)
    best, name = {}, {}
    for _, r in raw.iterrows():
        name[r.DDInterID_A] = r.Drug_A; name[r.DDInterID_B] = r.Drug_B
        k = r["key"]
        if k not in best or order[r.Level] > order[best[k]]:
            best[k] = r.Level
    G = nx.Graph(); G.add_edges_from(best.keys())

    deg = dict(G.degree()); nodes = list(deg)
    degv = np.array([deg[n] for n in nodes]); N = len(nodes)
    ni = {n: i for i, n in enumerate(nodes)}

    def topk(kf):
        m = round(kf * N); idx = np.argsort(degv)[::-1][:m]
        return m, degv[idx].sum() / degv.sum() * 100
    conc = {f"{int(k*100)}%": [topk(k)[0], round(topk(k)[1], 1)] for k in (.01,.05,.10,.20)}
    comps = sorted((len(c) for c in nx.connected_components(G)), reverse=True)

    ai = np.array([ni[a] for a, b in best]); bi = np.array([ni[b] for a, b in best])
    is_major = np.array([1 if best[(a, b)] == "Major" else 0 for a, b in best])
    def major_ep(flags):
        mc = np.zeros(N); np.add.at(mc, ai, flags); np.add.at(mc, bi, flags); return mc
    mc_obs = major_ep(is_major); tmaj = is_major.sum() * 2

    m10 = round(0.10 * N)
    od = np.argsort(degv)[::-1]; top, rest = od[:m10], od[m10:]
    cap_obs = mc_obs[top].sum() / tmaj * 100
    frac_top = mc_obs[top].sum() / degv[top].sum() * 100
    frac_rest = mc_obs[rest].sum() / degv[rest].sum() * 100
    rho_obs = spearmanr(degv, mc_obs).statistic

    rng = np.random.default_rng(SEED)
    obs10 = topk(.10)[1]; rs = np.empty(B); ds = degv.sum()
    for i in range(B):
        rs[i] = degv[rng.choice(N, m10, replace=False)].sum() / ds * 100
    E = len(is_major); nmaj = int(is_major.sum())
    tt, rt = degv[top].sum(), degv[rest].sum()
    cap_n = np.empty(B); diff_n = np.empty(B); rho_n = np.empty(B)
    for i in range(B):
        p = np.zeros(E); p[rng.choice(E, nmaj, replace=False)] = 1
        mc = major_ep(p)
        cap_n[i] = mc[top].sum() / (nmaj*2) * 100
        diff_n[i] = mc[top].sum()/tt*100 - mc[rest].sum()/rt*100
        rho_n[i] = spearmanr(degv, mc).statistic
    ep = lambda null, o: (np.sum(null >= o) + 1) / (B + 1)

    deg_by_name = {name[n].lower(): deg[n] for n in nodes}
    ranked = sorted(deg_by_name.items(), key=lambda x: -x[1])
    rank_of = {nm: i + 1 for i, (nm, _) in enumerate(ranked)}
    apex = {name[n]: int(deg[n]) for n in sorted(nodes, key=lambda x: -deg[x])[:10]}
    hub_rank, hubs_in_decile = {}, 0
    for h in PAPER_HUBS:
        cands = SYNONYMS.get(h, [h.lower()])
        found = [(deg_by_name[c], rank_of[c]) for c in cands if c in deg_by_name]
        if found:
            dg, rk = max(found); pct = round(rk / N * 100, 1)
            hub_rank[h] = {"degree": int(dg), "rank": int(rk), "percentile": pct}
            if pct <= 10: hubs_in_decile += 1

    res = {
      "pairs": G.number_of_edges(), "drugs": N,
      "severity_counts": {s: int(sum(1 for v in best.values() if v == s))
                          for s in ("Major","Moderate","Minor","Unknown")},
      "unknown_pct": round(sum(1 for v in best.values() if v=="Unknown")/G.number_of_edges()*100, 1),
      "mean_degree": round(float(degv.mean()), 1), "median_degree": float(np.median(degv)),
      "max_degree": int(degv.max()),
      "gini": round(gini(degv), 3), "density": round(nx.density(G), 4),
      "avg_clustering": round(nx.average_clustering(G), 3),
      "concentration": conc, "top10_endpoint_share": round(obs10, 1),
      "top20_endpoint_share": round(conc["20%"][1], 1),
      "n_components": len(comps), "giant_pct": round(comps[0]/N*100, 1),
      "random_set_mean": round(rs.mean(), 1), "conc_p": ep(rs, obs10),
      "top10_major_capture": round(cap_obs, 1),
      "major_capture_null_mean": round(cap_n.mean(), 1),
      "major_capture_null_ci": [round(np.percentile(cap_n, 2.5), 1),
                                round(np.percentile(cap_n, 97.5), 1)],
      "major_capture_p": ep(cap_n, cap_obs),
      "frac_major_top": round(frac_top, 1), "frac_major_rest": round(frac_rest, 1),
      "frac_diff_p": ep(diff_n, frac_top - frac_rest),
      "spearman_obs": round(rho_obs, 3), "spearman_null_mean": round(rho_n.mean(), 3),
      "spearman_null_ci": [round(np.percentile(rho_n, 2.5), 3),
                           round(np.percentile(rho_n, 97.5), 3)],
      "n_permutations": B, "seed": SEED,
      "apex_drugs_top10": apex,
      "paper_hub_ddinter_rank": hub_rank,
      "paper_hubs_in_top_decile": hubs_in_decile,
    }
    nulls = {"conc_obs": obs10, "cap_obs": cap_obs, "cap_null": cap_n,
             "rho_obs": rho_obs, "rho_null": rho_n, "conc_curve": conc,
             "rand_mean": rs.mean()}
    return res, nulls

def make_figure(res, nulls, outpng):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    GREY = "#8a8f94"
    fig, ax = plt.subplots(1, 3, figsize=(11, 3.5))
    ks = np.array([1, 5, 10, 20]); obs = [nulls["conc_curve"][f"{k}%"][1] for k in ks]
    ax[0].plot(ks, obs, "o-", color="#1b4965", lw=2, ms=6, zorder=3, label="DDInter (observed)")
    ax[0].plot(ks, ks, "--", color=GREY, lw=1.5, label="uniform expectation")
    ax[0].set_xlabel("Top X% of drugs by degree"); ax[0].set_ylabel("Share of interaction endpoints (%)")
    ax[0].set_title("Burden concentrates on hub drugs", loc="left")
    ax[0].annotate(f"{obs[2]:.0f}% on top 10%\n(random {nulls['rand_mean']:.0f}%, P<10$^{{-4}}$)",
                   xy=(10, obs[2]), xytext=(11, obs[2]-16),
                   arrowprops=dict(arrowstyle="->", color="#1b4965", lw=1), fontsize=6)
    ax[0].legend(frameon=False, fontsize=6, loc="upper left"); ax[0].margins(0.06)

    cap_n = nulls["cap_null"]; cap_obs = nulls["cap_obs"]
    ax[1].hist(cap_n, bins=40, color="#cfd8dc", zorder=1, label="label-permutation null")
    lo, hi = np.percentile(cap_n, [2.5, 97.5])
    ax[1].axvspan(lo, hi, color="#9aa7ad", alpha=0.3, zorder=0, label="null 95% interval")
    ax[1].axvline(cap_obs, color="#c1440e", lw=2, zorder=3, label="observed")
    p_cap = (np.sum(cap_n >= cap_obs)+1)/(len(cap_n)+1)
    ax[1].set_xlabel("Major-endpoint capture by top-decile drugs (%)"); ax[1].set_ylabel("Permutations")
    ax[1].set_title("Severe capture barely exceeds null", loc="left")
    ax[1].annotate(f"obs {cap_obs:.1f}%\nnull {cap_n.mean():.1f}%\nP={p_cap:.1e}",
                   xy=(cap_obs, ax[1].get_ylim()[1]*0.6),
                   xytext=(cap_obs-1.4, ax[1].get_ylim()[1]*0.55), ha="right",
                   fontsize=6, color="#c1440e")
    ax[1].legend(frameon=False, fontsize=6, loc="upper left")

    rho_n = nulls["rho_null"]; rho_obs = nulls["rho_obs"]
    ax[2].hist(rho_n, bins=40, color="#cfd8dc", zorder=1, label="null (part-whole)")
    ax[2].axvline(rho_obs, color="#c1440e", lw=2, zorder=3, label="observed")
    ax[2].set_xlabel(r"Spearman $\rho$ (degree vs # Major)"); ax[2].set_ylabel("Permutations")
    ax[2].set_title(r"Observed $\rho$ below part-whole null", loc="left")
    ax[2].annotate(f"obs {rho_obs:.2f}\nnull {rho_n.mean():.2f}",
                   xy=(rho_obs, ax[2].get_ylim()[1]*0.7),
                   xytext=(rho_obs+0.02, ax[2].get_ylim()[1]*0.65), fontsize=6, color="#c1440e")
    ax[2].legend(frameon=False, fontsize=6, loc="upper left")
    for a, L in zip(ax, "abc"):
        a.text(-0.08, 1.06, L, transform=a.transAxes, fontweight="bold", fontsize=11)
    fig.tight_layout(); fig.savefig(outpng, dpi=300, bbox_inches="tight")

def main(inp, outp, outpng):
    from pathlib import Path
    if not Path(inp).exists():
        raise FileNotFoundError(
            f"{inp} not found. Run prepare_ddinter_from_public.py --download first; "
            "the DDInter-derived pair file is intentionally not redistributed."
        )
    res, nulls = analyse(inp)
    json.dump(res, open(outp, "w"), indent=2)
    make_figure(res, nulls, outpng)
    print(json.dumps({k: res[k] for k in
        ("pairs","drugs","gini","top10_endpoint_share","top10_major_capture",
         "spearman_obs","paper_hubs_in_top_decile")}, indent=2))
    print("wrote", outp, "and", outpng)

if __name__ == "__main__":
    inp = sys.argv[1] if len(sys.argv) > 1 else "ddinter2_unique.csv"
    outp = sys.argv[2] if len(sys.argv) > 2 else "ddinter_corrected_results.json"
    outpng = sys.argv[3] if len(sys.argv) > 3 else "Fig6_ddinter_severity_corrected.png"
    main(inp, outp, outpng)
