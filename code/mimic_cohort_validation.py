#!/usr/bin/env python3
"""Inpatient-cohort validation in MIMIC-IV: does the hub watchlist reach real
co-administered medication pairs, and the co-administrations that would fire a
DDI alert?

Complements the NHANES analysis on two axes that NHANES cannot address:
1. Setting: inpatient prescriber orders, where DDI alerts fire in CPOE.
2. Concurrency: true temporal overlap of prescription windows within admission.

Two prespecified rules are applied before analysis: exclusion of non-drug orders
and a documented inpatient formulary synonym bridge. Nulls are 10,000 random
equal-size drug sets from the candidate network, seed 42.
"""
import argparse, collections, itertools, json, re
from pathlib import Path
import numpy as np
import pandas as pd
import networkx as nx

SEED = 42
N_RAND = 10_000
NON_DRUG_PATTERNS = [
    r"^bag$", r"^vial$", r"^syringe", r"^soln$", r"^sw$", r"^ns$", r"^d5w$",
    r"sterile water", r"sodium chloride.*flush", r"^0?\.9% sodium chloride",
    r"^sodium chloride 0\.9%", r"iso-osmotic", r"^lactated ringers", r"^dextrose",
    r"free water", r"heparin flush", r"flush$", r"contrast", r"omnipaque",
    r"^nutrition", r"^tpn", r"^parenteral nutrition", r"^dialysate", r"vaccine",
]
NON_DRUG_RE = re.compile("|".join(NON_DRUG_PATTERNS))
INPATIENT_SYNONYMS = {
    "insulin": "insulin human", "insulin human regular": "insulin regular",
    "docusate sodium": "docusate", "docusate": "dioctyl sulfosuccinate",
    "senna": "senna glycoside", "ipratropium": "ipratropium bromide",
    "albuterol sulfate": "albuterol", "milk of magnesia": "magnesium hydroxide",
    "vitamin d3": "cholecalciferol", "vitamin d2": "ergocalciferol",
    "vitamin b12": "cyanocobalamin", "metoprolol tartrate": "metoprolol",
    "amiodarone hcl": "amiodarone", "heparin sodium": "heparin",
    "enoxaparin sodium": "enoxaparin", "warfarin sodium": "warfarin",
    "levothyroxine sodium": "levothyroxine", "pantoprazole sodium": "pantoprazole",
    "furosemide sodium": "furosemide", "morphine sulfate": "morphine",
    "hydromorphone hcl": "hydromorphone", "aspirin ec": "aspirin",
}
SPLIT_RE = re.compile(r"\s*[-/;]\s*|\s+and\s+")
FORM_NOISE = re.compile(r"\b(neb|nebu|nebulizer|inhaler|inhalation|hfa|mdi|po|iv|im|sc|subq|pr|sl|oral|liquid|solution|soln|susp|tablet|tab|capsule|cap|cream|patch|gel|spray|ec|er|xl|xr|sr|cr|dr|extended|release|immediate|ir|premix|bolus|drip|infusion|inj|injection|human|sodium|hcl|hydrochloride|sulfate|succinate|tartrate|citrate|maleate|besylate|mesylate|fumarate|acetate|phosphate|bitartrate)\b")

def norm(s):
    s = re.sub(r"\(.*?\)", " ", str(s)).lower()
    s = re.sub(r"[^a-z0-9 /;-]", " ", s)
    return re.sub(r"\s+", " ", s).strip()

def resolve(raw, lookup, unmapped=None):
    k = norm(raw)
    if not k or NON_DRUG_RE.search(k): return set(), "nondrug"
    out = set(); parts = [p.strip() for p in SPLIT_RE.split(k) if p.strip()] or [k]
    for p in parts:
        cands = [INPATIENT_SYNONYMS.get(p), p]
        stripped = re.sub(r"\s+", " ", FORM_NOISE.sub(" ", p)).strip()
        cands += [INPATIENT_SYNONYMS.get(stripped), stripped, p.split(" ")[0]]
        for c in [x for x in cands if x]:
            db = lookup.get(c)
            if db: out.add(db); break
    if not out and unmapped is not None: unmapped[k] += 1
    return out, ("mapped" if out else "unmapped")

def load_candidate(path):
    df = pd.read_csv(path)
    ids = [c for c in df.columns if df[c].astype(str).str.match(r"^DB\d{5}$").mean() > 0.9]
    edges = sorted({tuple(sorted((a,b))) for a,b in zip(df[ids[0]],df[ids[1]]) if a != b})
    G=nx.Graph(); G.add_edges_from(edges); return G,set(edges)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--prescriptions",required=True); ap.add_argument("--patients"); ap.add_argument("--admissions"); ap.add_argument("--lookup",required=True); ap.add_argument("--candidate",required=True); ap.add_argument("--out",default="mimic_results.json"); a=ap.parse_args()
    lookup=json.load(open(a.lookup)); G,cand_edges=load_candidate(a.candidate); deg=dict(G.degree()); N=len(deg); ranked=[d for d,_ in sorted(deg.items(),key=lambda kv:(-kv[1],kv[0]))]; cand_nodes=set(deg)
    rx=pd.read_csv(a.prescriptions,low_memory=False); rx["starttime"]=pd.to_datetime(rx["starttime"],errors="coerce"); rx["stoptime"]=pd.to_datetime(rx["stoptime"],errors="coerce")
    unmapped=collections.Counter(); tally=collections.Counter(); cache={}; dbs=[]
    for raw in rx["drug"]:
        if raw not in cache: cache[raw]=resolve(raw,lookup,unmapped)
        ids,status=cache[raw]; tally[status]+=1; dbs.append(ids)
    rx["db_ids"]=dbs; rx_drug=rx[rx["db_ids"].map(len)>0].copy()
    overlap_pairs=collections.Counter(); admission_pairs=collections.Counter(); pt_overlap=collections.defaultdict(set); n_adm=0
    for hadm,g in rx_drug.groupby("hadm_id"):
        n_adm+=1; recs=[]
        for _,r in g.iterrows():
            for d in r["db_ids"]: recs.append((d,r["starttime"],r["stoptime"]))
        subj=g["subject_id"].iloc[0]; drugs_here={d for d,_,_ in recs}
        for p in itertools.combinations(sorted(drugs_here),2): admission_pairs[p]+=1
        by=collections.defaultdict(list)
        for d,st,sp in recs: by[d].append((st,sp))
        ds=sorted(by)
        for i in range(len(ds)):
            for j in range(i+1,len(ds)):
                found=False
                for s1,e1 in by[ds[i]]:
                    for s2,e2 in by[ds[j]]:
                        lo=max(s1,s2) if pd.notna(s1) and pd.notna(s2) else (s1 if pd.notna(s1) else s2)
                        hi=min(e1 if pd.notna(e1) else pd.Timestamp.max,e2 if pd.notna(e2) else pd.Timestamp.max)
                        if pd.isna(lo) or lo<=hi: found=True; break
                    if found: break
                if found:
                    pair=(ds[i],ds[j]); overlap_pairs[pair]+=1; pt_overlap[subj].update(pair)
    def restrict(counter): return {p:c for p,c in counter.items() if p[0] in cand_nodes and p[1] in cand_nodes}
    ov_in=restrict(overlap_pairs); ad_in=restrict(admission_pairs); ov_alert={p:c for p,c in ov_in.items() if p in cand_edges}; ad_alert={p:c for p,c in ad_in.items() if p in cand_edges}
    rng=np.random.default_rng(SEED); drugs_arr=np.array(ranked,dtype=object)
    def cover(counter,S,weighted=False):
        if not counter:return float("nan")
        if weighted:
            den=sum(counter.values()); return 100*sum(c for (x,y),c in counter.items() if x in S or y in S)/den
        return 100*sum(1 for x,y in counter if x in S or y in S)/len(counter)
    targets=[("overlapping_pairs",ov_in),("overlapping_alertable_pairs",ov_alert),("admission_pairs",ad_in),("admission_alertable_pairs",ad_alert)]
    res={"source":"MIMIC-IV Clinical Database Demo v2.2 (PhysioNet; ODbL v1.0; doi:10.13026/dp1f-ex47)","prescription_rows":int(len(rx)),"rows_nondrug_excluded":int(tally["nondrug"]),"rows_unmapped":int(tally["unmapped"]),"rows_mapped":int(tally["mapped"]),"mapping_rate_pct_of_drug_rows":round(100*tally["mapped"]/max(tally["mapped"]+tally["unmapped"],1),1),"admissions":n_adm,"subjects":int(rx_drug["subject_id"].nunique()),"unique_overlapping_pairs":len(overlap_pairs),"unique_admission_pairs":len(admission_pairs),"overlapping_pairs_in_candidate_network":len(ov_in),"admission_pairs_in_candidate_network":len(ad_in),"overlapping_pairs_alertable":len(ov_alert),"admission_pairs_alertable":len(ad_alert),"seed":SEED,"n_random_sets":N_RAND,"coverage":{}}
    for frac in (0.05,0.10,0.20):
        k=int(round(frac*N)); top=set(ranked[:k]); block={"k_drugs":k}
        for label,P in targets:
            obs=cover(P,top); obsw=cover(P,top,True); null=np.empty(N_RAND)
            for z in range(N_RAND): null[z]=cover(P,set(rng.choice(drugs_arr,k,replace=False)))
            block[label]={"n_pairs":len(P),"covered_pct":round(obs,1),"covered_pct_event_weighted":round(obsw,1),"null_mean_pct":round(float(np.nanmean(null)),1),"null_95_pct":[round(float(q),1) for q in np.nanquantile(null,[.025,.975])],"enrichment_vs_null":round(obs/np.nanmean(null),2) if np.nanmean(null) else None,"empirical_p":float((np.count_nonzero(null>=obs)+1)/(N_RAND+1))}
        res["coverage"][f"{int(frac*100)}%"] = block
    Path(a.out).write_text(json.dumps(res,indent=2)); print(json.dumps(res,indent=2))
if __name__=="__main__": main()
