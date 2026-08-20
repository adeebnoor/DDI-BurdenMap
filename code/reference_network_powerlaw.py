#!/usr/bin/env python3
"""
reference_network_powerlaw.py
Power-law fit for the construction-reference (GoldD3R) network, using the same
automatic-xmin procedure as the candidate network, plus an explicit
fitting-threshold sensitivity scan.

The manuscript rests NO conclusion on the fitted exponent: the load-bearing
result is the likelihood-ratio model comparison, reported for both networks
under one pinned environment (see requirements.txt). This script prints the
automatic fit and the sensitivity scan that the manuscript discloses.

Input : GoldD3R.txt   (lines: CUI_A,CUI_B[,[MechanismLabels]])
Usage : python reference_network_powerlaw.py GoldD3R.txt
"""
import sys, warnings
import numpy as np
import networkx as nx
import powerlaw

warnings.filterwarnings("ignore")


def compare(fit, tag):
    for alt in ("lognormal", "exponential"):
        try:
            R, p = fit.distribution_compare("power_law", alt, normalized_ratio=True)
            print(f"   {tag} vs {alt}: R={R:+.2f} P={p:.2e}")
        except Exception as exc:                       # pragma: no cover
            print(f"   {tag} vs {alt}: comparison unavailable ({exc})")


def main(path):
    G = nx.Graph()
    with open(path) as fh:
        for line in fh:
            parts = line.strip().split(",", 2)
            if len(parts) >= 2 and parts[0] and parts[1] and parts[0] != parts[1]:
                G.add_edge(parts[0], parts[1])
    deg = np.array([d for _, d in G.degree()])
    print(f"nodes={G.number_of_nodes()} edges={G.number_of_edges()}")

    fit = powerlaw.Fit(deg, discrete=True, verbose=False)
    xmin_auto = float(fit.xmin)
    print(f"\nAUTOMATIC xmin={xmin_auto:g} -> alpha={fit.alpha:.3f} "
          f"(tail n={(deg >= xmin_auto).sum()})   <-- value reported in the manuscript")
    compare(fit, "auto")

    print("\nFITTING-THRESHOLD SENSITIVITY (disclosed in Methods; no conclusion rests on alpha):")
    for xm in (10, 32, 50, 93):
        if xm == xmin_auto:
            continue                                   # already reported above
        tail_n = int((deg >= xm).sum())
        if tail_n < 50:
            print(f"   xmin={xm}: tail n={tail_n} (too sparse; skipped)")
            continue
        f2 = powerlaw.Fit(deg, discrete=True, xmin=xm, verbose=False)
        try:
            alpha = float(f2.power_law.alpha)
        except Exception:                              # pragma: no cover
            print(f"   xmin={xm}: fit unavailable in this package version")
            continue
        print(f"   xmin={xm}: alpha={alpha:.3f} (tail n={tail_n})")
        compare(f2, f"   xmin={xm}")

    print("\nNote: automatic threshold selection is implementation-dependent. Under an earlier "
          "release of the powerlaw package the automatic threshold for this network is xmin=93 "
          "(alpha=4.56), where the 149-node tail leaves the likelihood-ratio tests underpowered. "
          "The manuscript reports the values produced by the pinned environment in "
          "requirements.txt and states this sensitivity explicitly.")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "GoldD3R.txt")
