#!/usr/bin/env python3
"""Patient-level extension of the ONC expert-consensus analysis.

Counts ONC High-Priority and Non-Interruptive pairs actually co-prescribed in
NHANES and MIMIC-IV Demo, including co-exposure-event counts, watchlist coverage,
random-watchlist nulls and a label-permutation contrast on realized pairs.

Important: this patient-realization analysis uses the upstream fully expanded
`ONC_Non_Interuptive_Mapped.csv` representation. The knowledge-base Fig. 7
analysis uses the class-annotated formatted representation. The two denominators
are analysis-specific and are not compared directly.
"""
import argparse, collections, itertools, json, re
from pathlib import Path
import numpy as np
import pandas as pd
import networkx as nx

SEED=42; N_RAND=10_000; N_PERM=10_000
SALICYLATES={"DB00945","DB00861","DB01399","DB00936","DB09216","DB13743"}
BAD_SALICYLATE_IDS={"DB06709","DB06803"}
INVALID={"","nan","none","refused","don t know","dont know","unknown","55555","77777","99999"}

def norm(s):
    s=re.sub(r"\(.*?\)"," ",str(s)).lower(); s=re.sub(r"[^a-z0-9 ]"," ",s)
    return re.sub(r"\s+"," ",s).strip()

def load_onc(path):
    pairs=set()
    for line in open(path,encoding="utf-8-sig"):
        ids=re.findall(r"DB\d{5}",line)
        if len(ids)>=2 and ids[0]!=ids[1]: pairs.add(tuple(sorted(ids[:2])))
    return pairs

def correct_salicylates(pairs):
    out=set(); touched=0
    for a,b in pairs:
        bad={a,b}&BAD_SALICYLATE_IDS
        if not bad: out.add((a,b)); continue
        touched+=1; other=((({a,b}-bad) or {a}).pop())
        for s in SALICYLATES:
            if s!=other: out.add(tuple(sorted((s,other))))
    return out,touched

def load_candidate(path):
    df=pd.read_csv(path); ids=[c for c in df.columns if df[c].astype(str).str.match(r"^DB\d{5}$").mean()>0.9]
    edges=sorted({tuple(sorted((a,b))) for a,b in zip(df[ids[0]],df[ids[1]]) if a!=b})
    G=nx.Graph(); G.add_edges_from(edges); return G

def nhanes_pairs(files,lookup):
    per=collections.defaultdict(set)
    for f in files:
        df=pd.read_sas(f,format="xport")
        for c in df.columns:
            if df[c].dtype==object: df[c]=df[c].astype(str).str.replace(r"^b'|'$","",regex=True).str.strip()
        cyc=Path(f).stem
        for seqn,drug in zip(df["SEQN"],df.get("RXDDRUG",pd.Series(dtype=str))):
            k=norm(drug)
            if not k or k in INVALID: continue
            db=lookup.get(k) or lookup.get(k.split(" ")[0])
            if db: per[(cyc,int(seqn))].add(db)
    counter=collections.Counter()
    for drugs in per.values():
        if len(drugs)>=2:
            for p in itertools.combinations(sorted(drugs),2): counter[p]+=1
    return counter,len(per)

def mimic_pairs(path,lookup):
    import importlib.util
    spec=importlib.util.spec_from_file_location("mimicmod",str(Path(__file__).with_name("mimic_cohort_validation.py")))
    m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    rx=pd.read_csv(path,low_memory=False); rx["starttime"]=pd.to_datetime(rx["starttime"],errors="coerce"); rx["stoptime"]=pd.to_datetime(rx["stoptime"],errors="coerce")
    cache={}; ids=[]
    for raw in rx["drug"]:
        if raw not in cache: cache[raw]=m.resolve(raw,lookup)[0]
        ids.append(cache[raw])
    rx["db_ids"]=ids; rx=rx[rx["db_ids"].map(len)>0]
    counter=collections.Counter(); n_adm=0
    for _,g in rx.groupby("hadm_id"):
        n_adm+=1; by=collections.defaultdict(list)
        for _,r in g.iterrows():
            for d in r["db_ids"]: by[d].append((r["starttime"],r["stoptime"]))
        ds=sorted(by)
        for i in range(len(ds)):
            for j in range(i+1,len(ds)):
                found=False
                for s1,e1 in by[ds[i]]:
                    for s2,e2 in by[ds[j]]:
                        lo=max(s1,s2) if pd.notna(s1) and pd.notna(s2) else (s1 if pd.notna(s1) else s2)
                        e1x=e1 if pd.notna(e1) else pd.Timestamp.max; e2x=e2 if pd.notna(e2) else pd.Timestamp.max
                        if pd.isna(lo) or lo<=min(e1x,e2x): found=True; break
                    if found: break
                if found: counter[(ds[i],ds[j])]+=1
    return counter,n_adm

def analyse(name,pair_counter,n_units,hp,ni,ranked,N,rng):
    out={"cohort":name,"units":n_units,"unique_co_exposure_pairs":len(pair_counter),"co_exposure_events":int(sum(pair_counter.values()))}
    realised={}
    for label,L in (("high_priority",hp),("non_interruptive",ni)):
        pr={p:pair_counter[p] for p in L if p in pair_counter}; realised[label]=pr
        out[label]={"pairs_in_list":len(L),"pairs_realised_in_cohort":len(pr),"pct_of_list_realised":round(100*len(pr)/len(L),1) if L else None,"co_exposure_events":int(sum(pr.values()))}
    drugs=np.array(ranked,dtype=object)
    def ev_cov(pr,S):
        den=sum(pr.values()); return 100*sum(c for (x,y),c in pr.items() if x in S or y in S)/den if den else float("nan")
    def pair_cov(pr,S): return 100*sum(1 for x,y in pr if x in S or y in S)/len(pr) if pr else float("nan")
    out["coverage"]={}
    for frac in (0.05,0.10,0.20):
        k=int(round(frac*N)); top=set(ranked[:k]); blk={"k_drugs":k}
        for label in ("high_priority","non_interruptive"):
            pr=realised[label]; obs_e,obs_p=ev_cov(pr,top),pair_cov(pr,top); null=np.empty(N_RAND)
            for i in range(N_RAND): null[i]=ev_cov(pr,set(rng.choice(drugs,k,replace=False)))
            blk[label]={"n_pairs_realised":len(pr),"events":int(sum(pr.values())),"event_covered_pct":round(obs_e,1),"pair_covered_pct":round(obs_p,1),"null_mean_event_pct":round(float(np.nanmean(null)),1),"null_95_event_pct":[round(float(q),1) for q in np.nanquantile(null,[.025,.975])],"empirical_p":float((np.count_nonzero(null>=obs_e)+1)/(N_RAND+1))}
        blk["event_contrast_pctpts"]=round(blk["non_interruptive"]["event_covered_pct"]-blk["high_priority"]["event_covered_pct"],1)
        blk["pair_contrast_pctpts"]=round(blk["non_interruptive"]["pair_covered_pct"]-blk["high_priority"]["pair_covered_pct"],1)
        out["coverage"][f"{int(frac*100)}%"] = blk
    k=int(round(.10*N)); top=set(ranked[:k]); pooled=list(realised["high_priority"].items())+list(realised["non_interruptive"].items()); n_hp=len(realised["high_priority"])
    if pooled and n_hp and len(pooled)>n_hp:
        obs=ev_cov(realised["non_interruptive"],top)-ev_cov(realised["high_priority"],top); idx=np.arange(len(pooled)); null=np.empty(N_PERM)
        for i in range(N_PERM):
            rng.shuffle(idx); a={pooled[j][0]:pooled[j][1] for j in idx[:n_hp]}; b={pooled[j][0]:pooled[j][1] for j in idx[n_hp:]}; null[i]=ev_cov(b,top)-ev_cov(a,top)
        out["pooled_label_permutation"]={"pooled_realised_pairs":len(pooled),"n_assigned_high_priority":n_hp,"watchlist_pct":10,"observed_event_contrast_pctpts":round(obs,1),"null_mean_pctpts":round(float(null.mean()),2),"null_95_pctpts":[round(float(q),1) for q in np.quantile(null,[.025,.975])],"empirical_p":float((np.count_nonzero(np.abs(null)>=abs(obs))+1)/(N_PERM+1)),"n_permutations":N_PERM}
    return out

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--pddi",required=True); ap.add_argument("--candidate",required=True); ap.add_argument("--lookup",required=True); ap.add_argument("--nhanes-rx",nargs="*",default=[]); ap.add_argument("--mimic-rx",default=None); ap.add_argument("--out",default="cohort_onc_results.json"); a=ap.parse_args()
    lookup=json.load(open(a.lookup)); G=load_candidate(a.candidate); deg=dict(G.degree()); N=len(deg); ranked=[d for d,_ in sorted(deg.items(),key=lambda kv:(-kv[1],kv[0]))]; nodes=set(deg); d=Path(a.pddi)
    hp_raw=load_onc(d/"ONC-High-Priority"/"ONC_High_Priority_Mapped.csv"); ni_raw=load_onc(d/"ONC-Non-Interuptive"/"ONC_Non_Interuptive_Mapped.csv"); ni_corr,touched=correct_salicylates(ni_raw)
    restrict=lambda P:{p for p in P if p[0] in nodes and p[1] in nodes}; hp,ni=restrict(hp_raw),restrict(ni_corr)
    res={"onc_lists":{"high_priority_pairs_raw":len(hp_raw),"non_interruptive_pairs_raw":len(ni_raw),"salicylate_pairs_corrected":touched,"high_priority_in_candidate_network":len(hp),"non_interruptive_in_candidate_network":len(ni)},"seed":SEED,"n_random_sets":N_RAND,"n_permutations":N_PERM,"cohorts":{}}
    rng=np.random.default_rng(SEED)
    if a.nhanes_rx:
        pc,n=nhanes_pairs(a.nhanes_rx,lookup); res["cohorts"]["NHANES_2015_2018"]=analyse("NHANES 2015-2018 (ambulatory, participants)",pc,n,hp,ni,ranked,N,rng)
    if a.mimic_rx:
        pc,n=mimic_pairs(a.mimic_rx,lookup); res["cohorts"]["MIMIC_IV_demo"]=analyse("MIMIC-IV demo (inpatient, admissions; temporal overlap)",pc,n,hp,ni,ranked,N,rng)
    Path(a.out).write_text(json.dumps(res,indent=2)); print(json.dumps(res,indent=2)[:4000])

if __name__=="__main__": main()
