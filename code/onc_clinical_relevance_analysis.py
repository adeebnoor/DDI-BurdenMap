#!/usr/bin/env python3
"""Clinical-relevance discrimination test against the ONC expert-consensus DDI lists.

Compares a fixed degree-ranked watchlist against ONC High-Priority (should interrupt)
and Non-Interruptive (should not interrupt) reference pairs. Coverage is the share
of list pairs with at least one endpoint in the watchlist. Nulls use random equal-size
drug sets and a pooled expert-label permutation with the watchlist held fixed.

The public formatted Non-Interruptive file contains a documented Salicylates-class
mapping error; this analysis re-expands that class to true salicylates and retains
the uncorrected result as a sensitivity analysis.

Usage:
  python onc_clinical_relevance_analysis.py <candidate_csv> <pddi_datasets_dir> [out.json]
"""
import csv, json, re, sys
from pathlib import Path
import numpy as np
import pandas as pd
import networkx as nx

SEED = 42
N_RAND = 10_000
N_PERM = 10_000
SALICYLATES = {"DB00945", "DB00861", "DB01399", "DB00936", "DB09216", "DB13743"}
BAD_SALICYLATE_IDS = {"DB06709", "DB06803"}


def load_candidate(path):
    df = pd.read_csv(path)
    ids = [c for c in df.columns if df[c].astype(str).str.match(r"^DB\d{5}$").mean() > 0.9]
    edges = sorted({tuple(sorted((a, b))) for a, b in zip(df[ids[0]], df[ids[1]]) if a != b})
    G = nx.Graph(); G.add_edges_from(edges); return G


def load_high_priority(path):
    pairs = set()
    for line in open(path, encoding="utf-8-sig"):
        ids = re.findall(r"DB\d{5}", line)
        if len(ids) >= 2 and ids[0] != ids[1]: pairs.add(tuple(sorted(ids[:2])))
    return pairs


def load_non_interruptive(path, correct_salicylates=True):
    raw = []
    for row in csv.reader(open(path, encoding="utf-8-sig")):
        if len(row) >= 9: raw.append((row[1].strip(), row[3].strip(), row[6].strip(), row[8].strip()))
    pairs = set()
    for c1, d1, c2, d2 in raw:
        set1 = SALICYLATES if (correct_salicylates and c1.lower().startswith("salicylate")) else {d1}
        set2 = SALICYLATES if (correct_salicylates and c2.lower().startswith("salicylate")) else {d2}
        if correct_salicylates:
            set1 = {x for x in set1 if x not in BAD_SALICYLATE_IDS} or set1
            set2 = {x for x in set2 if x not in BAD_SALICYLATE_IDS} or set2
        for a in set1:
            for b in set2:
                if re.match(r"^DB\d{5}$", a) and re.match(r"^DB\d{5}$", b) and a != b:
                    pairs.add(tuple(sorted((a, b))))
    return pairs


def coverage(pairs, drug_set):
    return 100.0 * sum(1 for a, b in pairs if a in drug_set or b in drug_set) / len(pairs) if pairs else float("nan")


def main(cand_csv, pddi_dir, out_path="onc_results.json"):
    pddi = Path(pddi_dir); G = load_candidate(cand_csv); deg = dict(G.degree()); N = len(deg)
    ranked = [d for d, _ in sorted(deg.items(), key=lambda kv: (-kv[1], kv[0]))]
    hp_all = load_high_priority(pddi / "ONC-High-Priority" / "ONC_High_Priority_Mapped.csv")
    ni_all = load_non_interruptive(pddi / "ONC-Non-Interuptive" / "ONC_Non_Interuptive_List_Mapped_Formatted.csv", True)
    ni_raw = load_non_interruptive(pddi / "ONC-Non-Interuptive" / "ONC_Non_Interuptive_List_Mapped_Formatted.csv", False)
    inn = set(deg)
    hp = {p for p in hp_all if p[0] in inn and p[1] in inn}
    ni = {p for p in ni_all if p[0] in inn and p[1] in inn}
    ni_s = {p for p in ni_raw if p[0] in inn and p[1] in inn}
    res = {"candidate_network":{"drugs":N,"pairs":G.number_of_edges()},
           "onc_high_priority":{"pairs_total":len(hp_all),"pairs_mappable":len(hp)},
           "onc_non_interruptive":{"pairs_total":len(ni_all),"pairs_mappable":len(ni),"pairs_mappable_uncorrected":len(ni_s)},
           "salicylate_correction":"class 'Salicylates' re-expanded to true salicylates; source file mapped it to methacholine/niclosamide",
           "coverage_by_watchlist_size":{},"seed":SEED}
    rng=np.random.default_rng(SEED); drugs=np.array(ranked,dtype=object)
    for frac in (0.01,0.05,0.10,0.20):
        k=int(round(frac*N)); top=set(ranked[:k]); obs_hp,obs_ni=coverage(hp,top),coverage(ni,top)
        null_hp=np.empty(N_RAND); null_ni=np.empty(N_RAND)
        for i in range(N_RAND):
            S=set(rng.choice(drugs,k,replace=False)); null_hp[i]=coverage(hp,S); null_ni[i]=coverage(ni,S)
        res["coverage_by_watchlist_size"][f"{int(frac*100)}%"]={
            "k_drugs":k,"high_priority_pct":round(obs_hp,1),"high_priority_null_mean_pct":round(float(null_hp.mean()),1),
            "high_priority_null_95":[round(float(q),1) for q in np.quantile(null_hp,[.025,.975])],
            "high_priority_p":float((np.count_nonzero(null_hp>=obs_hp)+1)/(N_RAND+1)),
            "non_interruptive_pct":round(obs_ni,1),"non_interruptive_null_mean_pct":round(float(null_ni.mean()),1),
            "non_interruptive_null_95":[round(float(q),1) for q in np.quantile(null_ni,[.025,.975])],
            "non_interruptive_p":float((np.count_nonzero(null_ni>=obs_ni)+1)/(N_RAND+1)),
            "contrast_pctpts":round(obs_ni-obs_hp,1)}
    k=int(round(0.10*N)); top=set(ranked[:k]); pooled=sorted(hp|ni); n_hp=len({p for p in pooled if p in hp})
    hit=np.array([1 if (a in top or b in top) else 0 for a,b in pooled]); obs=res["coverage_by_watchlist_size"]["10%"]["contrast_pctpts"]
    perm=np.empty(N_PERM); idx=np.arange(len(pooled))
    for i in range(N_PERM):
        rng.shuffle(idx); a,b=idx[:n_hp],idx[n_hp:]; perm[i]=100*hit[b].mean()-100*hit[a].mean()
    res["pooled_label_permutation"]={"pooled_pairs":len(pooled),"n_assigned_high_priority":n_hp,"observed_contrast_pctpts":obs,
        "null_mean_pctpts":round(float(perm.mean()),2),"null_95_pctpts":[round(float(q),1) for q in np.quantile(perm,[.025,.975])],
        "empirical_p":float((np.count_nonzero(np.abs(perm)>=abs(obs))+1)/(N_PERM+1)),"n_permutations":N_PERM}
    res["sensitivity_uncorrected_salicylates"]={"non_interruptive_pct_top10":round(coverage(ni_s,top),1)}
    Path(out_path).write_text(json.dumps(res,indent=2)); print(json.dumps(res,indent=2))


if __name__ == "__main__":
    main(sys.argv[1],sys.argv[2],sys.argv[3] if len(sys.argv)>3 else "onc_results.json")
