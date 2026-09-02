# Patient-cohort transport analysis

The Pharmaceutics submission includes a primary ambulatory medication-exposure transport test and a secondary inpatient transport/sensitivity analysis. Patient data are introduced only after the candidate-network hub ranking and watchlist have been fixed.

## NHANES 2015-2018 — primary ambulatory transport test

Public CDC/NCHS prescription-medication files:
- `RXQ_RX_I.XPT` (2015-2016)
- `RXQ_RX_J.XPT` (2017-2018)
- `DEMO_I.XPT` and `DEMO_J.XPT` for cohort description

The analytical cohort contains 7,669 participants with at least one mapped prescription; 5,301 used two or more mapped drugs. The **primary H3 test** uses all 17,229 observed drug pairs whose drugs are rankable in the candidate network; the observed pair is not required to be a candidate-network edge. At the frozen 10% watchlist, 29.2% of these pairs were covered versus a 16.6% random equal-size drug-set mean (95% interval 10.8-22.8%; 10,000 draws; seed 42; empirical P = 3.0e-4).

Restricting to the 1,133 observed pairs that are also candidate-network edges gives 70.0% coverage. Because the same top-decile watchlist already covers 65.0% of all candidate-network edges before patient data are introduced, this 70.0% value is reported only as **operational reach**, not independent validation.

## MIMIC-IV Clinical Database Demo v2.2 — secondary inpatient transport/sensitivity analysis

PhysioNet DOI: `10.13026/dp1f-ex47`. The open demo contains 100 de-identified patients. The analysis contains 250 admissions with at least one mapped prescription order. Co-exposure is defined by temporal overlap of **prescription-order** start/stop windows within one admission; a same-admission sensitivity definition is also computed. No administration-event claim is made.

Temporal overlap produced 14,677 unique mapped co-exposed pairs; 9,270 had both drugs rankable in the candidate network. At the frozen 10% watchlist, 30.8% of these rankable observed pairs were covered versus a 19.0% random equal-size drug-set mean (95% interval 11.7-26.9%; 10,000 draws; seed 42; empirical P = 0.0018).

Restricting to the 704 overlapping pairs that are also candidate-network edges gives 68.0% coverage; this is again interpreted as operational reach because of the unconditional 65.0% candidate-edge baseline.

Because the MIMIC-IV Demo is intentionally small, it is not presented as a second large definitive clinical cohort. It tests transport to an inpatient prescription-order stream under a stricter concurrency definition.

## ONC expert-consensus analyses

Two ONC analyses use two mapped representations from the same pinned upstream PDDI repository commit:

- **Knowledge-base coverage (Fig. 7):** the class-annotated formatted Non-Interruptive representation; after correction and candidate-network restriction, 1,895 Non-Interruptive pairs.
- **Patient realization (Fig. 8b):** the fully expanded mapped-pair representation; after the corresponding correction/restriction, 2,025 Non-Interruptive pairs.

These denominators are analysis-specific and are not compared directly. `prepare_onc_from_public.py` fetches both representations plus the High-Priority list from pinned commit `8199ee66b60bcb337f777889a210dd0d72a96e8f`.

Across patient cohorts, Non-Interruptive pairs are realized 6.9-fold more often than High-Priority pairs in NHANES and 3.7-fold more often in MIMIC-IV Demo. The **watchlist-specific** ONC coverage contrast among realized pairs is underpowered/non-significant (P = 0.26 and P = 0.80) and is not claimed as replicated.

## Reproducibility and data boundary

Raw NHANES and MIMIC-IV tables are not committed. Download them from their source custodians and run the scripts locally. NHANES public-use materials are federal public-use data (generally public domain, subject to NCHS requirements). MIMIC-IV Demo is distributed under ODbL v1.0. Aggregate results are provided under `out/`.

The fixed drug-name-to-DrugBank mapping is redistributable and is stored as `data/name_to_drugbank.json` in the repository/reproducibility archive.
