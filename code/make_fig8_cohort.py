#!/usr/bin/env python3
"""Figure 8 — transport of a frozen drug-level DDI-burden ranking into
observed medication exposure.

Panel A shows the *primary* H3 transport test: coverage of observed co-exposed
pairs without conditioning on candidate-edge membership. The two cohorts use
different denominators, by construction of the two cohort scripts:
  NHANES   all unique co-taken pairs among DrugBank-mapped drugs (17,229);
           no rankability restriction is applied.
  MIMIC-IV overlapping pairs whose two drugs are both candidate-network
           nodes (9,270 of 14,677).
In each cohort the observed value and its random-set null are computed on that
cohort's own pair set, so the within-cohort contrast is valid; the absolute
percentages are not on a common denominator. Panel B shows the separate ONC
list-level realization analysis.

Usage: python make_fig8_cohort.py nhanes_results.json mimic_results.json \
           cohort_onc_results.json outdir/
"""
import json, os, sys
from pathlib import Path
import numpy as np, matplotlib as mpl, matplotlib.pyplot as plt

# Render resolution: default 300 reproduces the archived reference PNGs exactly.
RENDER_DPI = int(os.environ.get("FIGURE_DPI", "300"))

mpl.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8,
    "axes.titlesize": 8, "axes.labelsize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 6,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": RENDER_DPI, "savefig.dpi": RENDER_DPI,
})
META_GREY = "#8A8A8A"
C_NH, C_MI = "#C1652F", "#1F6F78"

nh_p, mi_p, co_p, outdir = sys.argv[1], sys.argv[2], sys.argv[3], Path(sys.argv[4])
outdir.mkdir(parents=True, exist_ok=True)
nh, mi, co = (json.load(open(p)) for p in (nh_p, mi_p, co_p))

sizes = ["5%", "10%", "20%"]
x = np.arange(len(sizes))
# Primary H3 targets (see docstring: NHANES unrestricted, MIMIC-IV rankable-restricted).
NH_TARGET = [nh["coverage"][s]["all_concurrent_pairs"] for s in sizes]
MI_TARGET = [mi["coverage"][s]["overlapping_pairs"] for s in sizes]

# Panel B inputs: list-level realization rates.
def rates(key):
    c = co["cohorts"][key]
    hp, ni = c["high_priority"], c["non_interruptive"]
    # Use the rounded list-level percentages stored in the audited output so
    # the figure exactly matches Table 5 and the manuscript text.
    r1 = float(hp["pct_of_list_realised"])
    r2 = float(ni["pct_of_list_realised"])
    rr = (ni["pairs_realised_in_cohort"] / ni["pairs_in_list"]) / (
          hp["pairs_realised_in_cohort"] / hp["pairs_in_list"])
    return r1, r2, round(rr, 1)

nh_hp, nh_ni, nh_rr = rates("NHANES_2015_2018")
mi_hp, mi_ni, mi_rr = rates("MIMIC_IV_demo")
hp_r, ni_r, rr = [nh_hp, mi_hp], [nh_ni, mi_ni], [nh_rr, mi_rr]


def panel_letter(ax, letter):
    ax.text(-0.16, 1.06, letter, transform=ax.transAxes,
            fontsize=11, fontweight="bold", va="bottom", ha="left")


def build_fig8():
    fig = plt.figure(figsize=(7.4, 3.5))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.34, 1], wspace=0.42)
    axA = fig.add_subplot(gs[0, 0]); axB = fig.add_subplot(gs[0, 1])

    for arr, col, lab, mk in [
        (NH_TARGET, C_NH, "NHANES (ambulatory)", "o"),
        (MI_TARGET, C_MI, "MIMIC-IV (inpatient)", "s"),
    ]:
        obs = [b["covered_pct"] for b in arr]
        null = [b["null_mean_pct"] for b in arr]
        lo = [b["null_95_pct"][0] for b in arr]
        hi = [b["null_95_pct"][1] for b in arr]
        axA.fill_between(x, lo, hi, color=col, alpha=0.10, lw=0, zorder=1)
        axA.plot(x, null, color=col, ls="--", lw=1.0, alpha=0.75, zorder=2)
        axA.plot(x, obs, color=col, marker=mk, lw=2.0, ms=6, zorder=3, label=lab)

    axA.set_xticks(x); axA.set_xticklabels(sizes)
    axA.set_xlabel("Frozen watchlist size (top X% by candidate-pair degree)")
    axA.set_ylabel("Observed medication pairs covered (%)")
    axA.set_title("Frozen drug ranking is enriched in observed\nmedication pairings in both settings", loc="left")
    axA.set_ylim(0, 65); axA.margins(x=0.09)
    axA.annotate("29.2%", (x[1], NH_TARGET[1]["covered_pct"]),
                 textcoords="offset points", xytext=(-43, 10),
                 color=C_NH, fontweight="bold", ha="left")
    axA.annotate("30.8%", (x[1], MI_TARGET[1]["covered_pct"]),
                 textcoords="offset points", xytext=(10, 3),
                 color=C_MI, fontweight="bold", ha="left")
    axA.text(0.02, 0.02,
             "NHANES: all mapped co-taken pairs; MIMIC-IV: candidate-universe overlaps\n"
             "Dashed / shaded: cohort-specific random-set mean / 95% interval",
             transform=axA.transAxes, fontsize=5.2, color=META_GREY, va="bottom")
    axA.legend(frameon=False, loc="upper left", fontsize=6)
    panel_letter(axA, "a")

    xb = np.arange(2); w = 0.33
    b1 = axB.bar(xb - w/2, hp_r, w, color=C_MI)
    b2 = axB.bar(xb + w/2, ni_r, w, color=C_NH)
    for xi, v in zip(xb - w/2, hp_r):
        axB.annotate(f"{v:.1f}%", (xi, v), ha="center", va="bottom",
                     textcoords="offset points", xytext=(0, 3), fontsize=6)
    for xi, v, rv in zip(xb + w/2, ni_r, rr):
        axB.annotate(f"{v:.1f}%  ({rv}×)", (xi, v), ha="center", va="bottom",
                     textcoords="offset points", xytext=(0, 3), fontsize=6, fontweight="bold")
    axB.set_xticks(xb); axB.set_xticklabels(["NHANES\n(ambulatory)", "MIMIC-IV\n(inpatient)"])
    axB.set_ylabel("ONC list pairs realized in cohort (%)")
    axB.set_title("Consensus DDI categories differ in\nreal-world realization frequency", loc="left")
    axB.set_ylim(0, 16); axB.margins(x=0.20)
    axB.legend([b1, b2], ["ONC high-priority", "ONC non-interruptive"],
               frameon=False, fontsize=5.8, loc="upper center",
               bbox_to_anchor=(0.5, -0.20), ncol=1, handlelength=1.1)
    panel_letter(axB, "b")
    return fig


fig = build_fig8()
fig.savefig(outdir / "Fig8_cohort_validation.png", dpi=RENDER_DPI, bbox_inches="tight")
fig.savefig(outdir / "Fig8_cohort_validation.tif", dpi=RENDER_DPI, bbox_inches="tight",
            pil_kwargs={"compression": "tiff_lzw"})
print("wrote", outdir / "Fig8_cohort_validation.png")
