#!/usr/bin/env python3
"""Mechanism-stratified topology of the labelled GoldD3R reference network.

For each interaction-mechanism label in GoldD3R.txt, constructs the subnetwork
formed by pairs carrying that label and reports nodes, edges, Gini coefficient,
top-decile endpoint share, and giant-component fraction. A pair with multiple
labels contributes to each corresponding mechanism stratum.

Usage:
  python mechanism_stratified_analysis.py GoldD3R.txt mechanism_topology.csv mechanism_topology.json
"""
import collections, csv, json, re, sys
from pathlib import Path
import networkx as nx
import numpy as np


def gini(values):
    v = np.sort(np.asarray(values, dtype=float))
    n = len(v)
    if n == 0 or v.sum() == 0:
        return float('nan')
    return float((2 * np.sum(np.arange(1, n + 1) * v)) / (n * v.sum()) - (n + 1) / n)


def metrics(edges):
    G = nx.Graph(); G.add_edges_from(sorted(edges))
    d = np.array(sorted(dict(G.degree()).values(), reverse=True), dtype=float)
    n, e = G.number_of_nodes(), G.number_of_edges()
    k = max(1, int(round(0.10 * n)))
    comps = sorted(nx.connected_components(G), key=len, reverse=True)
    return {
        'nodes': n,
        'edges': e,
        'gini': round(gini(d), 3),
        'top10_endpoint_share_pct': round(float(100 * d[:k].sum() / d.sum()), 1),
        'giant_component_pct': round(float(100 * len(comps[0]) / n), 1) if n else 0.0,
    }


def main(inp, out_csv, out_json):
    strata = collections.defaultdict(set)
    for ln in Path(inp).read_text().splitlines():
        m = re.match(r'^(C\d+),(C\d+)(?:,\[([^\]]*)\])?', ln.strip())
        if not m or m.group(1) == m.group(2) or not m.group(3):
            continue
        edge = tuple(sorted((m.group(1), m.group(2))))
        for lab in re.split(r'[,/;|]', m.group(3)):
            lab = lab.strip()
            if lab:
                strata[lab].add(edge)
    rows = []
    for label, edges in strata.items():
        row = {'mechanism': label, **metrics(edges)}
        rows.append(row)
    rows.sort(key=lambda r: (-r['edges'], r['mechanism']))
    Path(out_json).write_text(json.dumps(rows, indent=2))
    with open(out_csv, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['mechanism','nodes','edges','gini','top10_endpoint_share_pct','giant_component_pct'])
        w.writeheader(); w.writerows(rows)
    print(json.dumps(rows, indent=2))

if __name__ == '__main__':
    inp = sys.argv[1] if len(sys.argv) > 1 else 'GoldD3R.txt'
    out_csv = sys.argv[2] if len(sys.argv) > 2 else 'mechanism_topology.csv'
    out_json = sys.argv[3] if len(sys.argv) > 3 else 'mechanism_topology.json'
    main(inp, out_csv, out_json)
