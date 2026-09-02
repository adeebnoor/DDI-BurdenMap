# DDI-BurdenMap v1.1.0 — Pharmaceutics framing lock

Released: 2 September 2026

This release aligns the public repository with the final Pharmaceutics manuscript framing.

## What changed

- Reframed the primary patient analysis as **transport of a frozen drug-level burden structure into observed medication exposure**, rather than generic "patient validation".
- Promoted the non-conditioned H3 test to the primary patient result:
  - NHANES: 29.2% coverage vs 16.6% random-set mean (P = 3.0e-4).
  - MIMIC-IV Demo: 30.8% vs 19.0% (P = 0.0018).
- Explicitly documented the 65.0% unconditional candidate-edge baseline; 70.0%/68.0% candidate-edge cohort coverage is therefore labeled **secondary operational reach**, not independent validation.
- Clarified MIMIC-IV as temporal overlap of **prescription-order** windows, not medication-administration events.
- Retained the non-significant watchlist-specific ONC realized-pair contrasts as negative results.
- Updated citation/Zenodo metadata and landing-page language around mechanism-aware drug-level burden concentration.

The submission reproducibility archive remains the version of record for exact redistributable fixed inputs, scripts, environment, and machine-readable outputs.
