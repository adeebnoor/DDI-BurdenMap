#!/usr/bin/env python3
"""Patient-cohort validation in NHANES: does the fixed hub watchlist reach real
concurrent medication use?

The watchlist is derived from the candidate interaction network before patient
data are introduced. This script asks what share of actually co-taken medication
pairs, including the subset also present in that prespecified candidate network,
involve a watchlist drug. Candidate-network membership is a structural target;
it does not imply a confirmed DDI or an alert that actually fired.

Inputs (NHANES public files, .XPT):
  RXQ_RX_I.XPT   Prescription medications, 2015-2016
  RXQ_RX_J.XPT   Prescription medications, 2017-2018
  DEMO_I / DEMO_J (optional)  Demographics, for cohort description

Also required:
  name_to_drugbank.json       drug-name -> DrugBank ID lookup
  d3_candidate_ddi_pairs.csv candidate network (defines degree/watchlist)
Optional:
  ddinter2_unique.csv         independent curated-pair sensitivity target

All randomness is seeded (42).
"""
import argparse, collections, itertools, json, re
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd

SEED = 42
N_RAND = 10_000
INVALID = {"", "nan", "none", "refused", "don t know", "dont know", "unknown",
           "55555", "77777", "99999"}


def norm(s):
    s = re.sub(r"\(.*?\)", " ", str(s)).lower()
    s = re.sub(r"[^a-z0-9 ]", " ", s)
    return re.sub(r"\s+", " ", s).strip()


def read_xpt(path):
    df = pd.read_sas(path, format="xport")
    for c in df.columns:
        if df[c].dtype == object:
            df[c] = df[c].astype(str).str.replace(r"^b'|'$", "", regex=True).str.strip()
    return df


def load_candidate(path):
    df = pd.read_csv(path)
    ids = [c for c in df.columns if df[c].astype(str).str.match(r"^DB\d{5}$").mean() > 0.9]
    edges = sorted({tuple(sorted((a, b))) for a, b in zip(df[ids[0]], df[ids[1]]) if a != b})
    G = nx.Graph(); G.add_edges_from(edges)
    return G, set(edges)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rx", nargs="+", required=True)
    ap.add_argument("--demo", nargs="*", default=[])
    ap.add_argument("--lookup", required=True)
    ap.add_argument("--candidate", required=True)
    ap.add_argument("--ddinter", default=None)
    ap.add_argument("--out", default="nhanes_results.json")
    a = ap.parse_args()

    lookup = json.load(open(a.lookup))
    G, cand_edges = load_candidate(a.candidate)
    deg = dict(G.degree()); N = len(deg)
    ranked = [d for d, _ in sorted(deg.items(), key=lambda kv: (-kv[1], kv[0]))]

    per_person = collections.defaultdict(set)
    rows_total = rows_named = rows_mapped = 0
    unmapped = collections.Counter(); cycles = {}
    for f in a.rx:
        df = read_xpt(f); cyc = Path(f).stem; seqns = set()
        for seqn, drug in zip(df["SEQN"], df.get("RXDDRUG", pd.Series(dtype=str))):
            rows_total += 1; k = norm(drug)
            if not k or k in INVALID: continue
            rows_named += 1
            db = lookup.get(k) or lookup.get(k.split(" ")[0])
            if db is None:
                unmapped[k] += 1; continue
            rows_mapped += 1
            per_person[(cyc, int(seqn))].add(db); seqns.add(int(seqn))
        cycles[cyc] = {"records": int(len(df)), "participants_with_rx": len(seqns)}

    pair_counter = collections.Counter(); n_multi = 0
    for drugs in per_person.values():
        if len(drugs) >= 2:
            n_multi += 1
            for p in itertools.combinations(sorted(drugs), 2): pair_counter[p] += 1
    pairs = set(pair_counter); total_events = sum(pair_counter.values())
    candidate_pairs = {p for p in pairs if p in cand_edges}

    ddinter_pairs = set()
    if a.ddinter:
        dd = pd.read_csv(a.ddinter, dtype=str); dd_edges = set()
        for n1, n2 in zip(dd.get("Drug_A", []), dd.get("Drug_B", [])):
            d1, d2 = lookup.get(norm(n1)), lookup.get(norm(n2))
            if d1 and d2 and d1 != d2: dd_edges.add(tuple(sorted((d1, d2))))
        ddinter_pairs = {p for p in pairs if p in dd_edges}

    rng = np.random.default_rng(SEED); drugs_arr = np.array(ranked, dtype=object)

    def cover_pairs(P, S, weighted=False):
        if not P: return float("nan")
        if weighted:
            den = sum(c for p, c in pair_counter.items() if p in P)
            num = sum(c for p, c in pair_counter.items() if p in P and (p[0] in S or p[1] in S))
            return 100 * num / den if den else float("nan")
        return 100 * sum(1 for x, y in P if x in S or y in S) / len(P)

    res = {
        "cycles": cycles, "records_total": rows_total,
        "records_with_drug_name": rows_named, "records_mapped_to_drugbank": rows_mapped,
        "mapping_rate_pct": round(100 * rows_mapped / max(rows_named, 1), 1),
        "top_unmapped_names": unmapped.most_common(25),
        "participants_with_any_rx": len(per_person), "participants_on_2plus_drugs": n_multi,
        "unique_concurrent_pairs": len(pairs), "concurrent_pair_events": total_events,
        "pairs_in_candidate_network": len(candidate_pairs), "pairs_in_ddinter": len(ddinter_pairs),
        "seed": SEED, "n_random_sets": N_RAND, "coverage": {}
    }

    targets = [("all_concurrent_pairs", pairs), ("candidate_network_pairs", candidate_pairs)]
    if ddinter_pairs: targets.append(("ddinter_pairs", ddinter_pairs))

    for frac in (0.05, 0.10, 0.20):
        k = int(round(frac * N)); top = set(ranked[:k]); block = {"k_drugs": k}
        for label, P in targets:
            obs = cover_pairs(P, top); obs_w = cover_pairs(P, top, weighted=True)
            null = np.empty(N_RAND)
            for i in range(N_RAND): null[i] = cover_pairs(P, set(rng.choice(drugs_arr, k, replace=False)))
            block[label] = {
                "n_pairs": len(P), "covered_pct": round(obs, 1),
                "covered_pct_event_weighted": round(obs_w, 1),
                "null_mean_pct": round(float(null.mean()), 1),
                "null_95_pct": [round(float(q), 1) for q in np.quantile(null, [.025, .975])],
                "enrichment_vs_null": round(obs / null.mean(), 2) if null.mean() else None,
                "empirical_p": float((np.count_nonzero(null >= obs) + 1) / (N_RAND + 1)),
            }
        exposed = sum(1 for d in per_person.values() if d & top)
        block["participants_exposed_to_watchlist_drug_pct"] = round(100 * exposed / len(per_person), 1)
        res["coverage"][f"{int(frac*100)}%"] = block

    if a.demo:
        with_rx = {seqn for _, seqn in per_person}
        multi = {seqn for (_, seqn), drugs in per_person.items() if len(drugs) >= 2}
        frames = []
        for f in a.demo:
            d = read_xpt(f); frames.append(pd.DataFrame({c: d[c] for c in ("SEQN", "RIDAGEYR", "RIAGENDR") if c in d}))
        demo = pd.concat(frames, ignore_index=True)
        demo["SEQN"] = pd.to_numeric(demo["SEQN"], errors="coerce").astype("Int64")
        demo["RIDAGEYR"] = pd.to_numeric(demo["RIDAGEYR"], errors="coerce")
        demo["RIAGENDR"] = pd.to_numeric(demo["RIAGENDR"], errors="coerce")
        def describe(sub, label):
            ages = sub["RIDAGEYR"].dropna(); sexes = sub["RIAGENDR"].dropna()
            return {"population": label, "n": int(len(sub)),
                    "age_median": float(ages.median()) if len(ages) else None,
                    "age_iqr": [float(ages.quantile(.25)), float(ages.quantile(.75))] if len(ages) else None,
                    "age_65plus_pct": round(100 * float((ages >= 65).mean()), 1) if len(ages) else None,
                    "pct_female": round(100 * float((sexes == 2).mean()), 1) if len(sexes) else None}
        res["cohort"] = {
            "all_participants": describe(demo, "all NHANES participants (both cycles)"),
            "with_any_prescription": describe(demo[demo["SEQN"].isin(with_rx)], "participants with >=1 mapped prescription"),
            "on_2plus_drugs": describe(demo[demo["SEQN"].isin(multi)], "participants on >=2 mapped drugs (contribute pairs)")
        }

    Path(a.out).write_text(json.dumps(res, indent=2))
    print(json.dumps({k: v for k, v in res.items() if k != "top_unmapped_names"}, indent=2))


if __name__ == "__main__": main()
