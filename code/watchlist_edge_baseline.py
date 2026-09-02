#!/usr/bin/env python3
"""Compute unconditional candidate-edge coverage of degree-ranked drug watchlists."""
import argparse, json
import pandas as pd
import networkx as nx

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('candidate')
    ap.add_argument('--out', default='out/watchlist_edge_baseline.json')
    a=ap.parse_args()
    df=pd.read_csv(a.candidate)
    cols=[c for c in df.columns if df[c].astype(str).str.match(r'^DB\d{5}$').mean()>0.9][:2]
    edges=sorted({tuple(sorted((str(x),str(y)))) for x,y in zip(df[cols[0]],df[cols[1]]) if str(x)!=str(y)})
    G=nx.Graph(); G.add_edges_from(edges)
    ranked=[d for d,_ in sorted(G.degree(), key=lambda kv:(-kv[1],kv[0]))]
    res={'n_drugs':G.number_of_nodes(),'n_edges':G.number_of_edges(),'coverage':{}}
    for frac in (0.05,0.10,0.20):
        k=round(frac*G.number_of_nodes()); S=set(ranked[:k])
        covered=sum(1 for u,v in edges if u in S or v in S)
        res['coverage'][f'{int(frac*100)}%']={'k_drugs':k,'covered_edges':covered,'covered_pct':100*covered/len(edges)}
    with open(a.out,'w') as f: json.dump(res,f,indent=2)
    print(json.dumps(res,indent=2))
if __name__=='__main__': main()
