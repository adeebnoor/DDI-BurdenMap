# Patient-cohort validation

The Frontiers resubmission adds a **primary independent ambulatory patient-cohort validation** and a **secondary inpatient transport/sensitivity analysis**. Patient data are introduced only after the candidate-network hub ranking and watchlist have been fixed.

## NHANES 2015-2018 — primary validation

Public CDC/NCHS prescription-medication files:
- `RXQ_RX_I.XPT` (2015-2016)
- `RXQ_RX_J.XPT` (2017-2018)
- `DEMO_I.XPT` and `DEMO_J.XPT` for cohort description

The analytical cohort contains 7,669 participants with at least one mapped prescription; 5,301 used two or more mapped drugs. At the 10% watchlist, 70.0% of **co-taken candidate-network pairs** were covered versus an 18.9% random-watchlist null (10,000 draws; seed 42; empirical P < 1e-4).

Candidate-network membership is a prespecified structural target. It does **not** mean every pair is a confirmed interaction or that an alert actually fired.

## MIMIC-IV Clinical Database Demo v2.2 — secondary inpatient sensitivity analysis

PhysioNet DOI: `10.13026/dp1f-ex47`. The open demo contains 100 de-identified patients and 250 admissions. Prescriptions are evaluated using temporal overlap of start/stop windows within admission; a same-admission sensitivity definition is also computed.

At the 10% watchlist, 68.0% of temporally overlapping **co-exposed candidate-network pairs** were covered versus an 18.9% random-watchlist null (10,000 draws; seed 42; empirical P < 1e-4).

Because the MIMIC-IV Demo is intentionally small, it is not presented as a second large definitive clinical cohort. NHANES is the primary patient validation; MIMIC-IV Demo tests transport to an inpatient order stream under a stricter concurrency definition.

## ONC expert-consensus analyses

Two ONC analyses use two mapped representations from the same pinned upstream PDDI repository commit:

- **Knowledge-base coverage (Fig. 7):** the class-annotated formatted Non-Interruptive representation; after correction and candidate-network restriction, 1,895 Non-Interruptive pairs.
- **Patient realization (Fig. 8b):** the fully expanded mapped-pair representation; after the corresponding correction/restriction, 2,025 Non-Interruptive pairs.

These denominators are analysis-specific and are not compared directly. `prepare_onc_from_public.py` fetches both representations plus the High-Priority list from pinned commit `8199ee66b60bcb337f777889a210dd0d72a96e8f`.

Across patient cohorts, Non-Interruptive pairs are realized 6.9-fold more often than High-Priority pairs in NHANES and 3.7-fold more often in MIMIC-IV Demo. The **watchlist-specific** ONC coverage contrast among realized pairs is underpowered/non-significant (P=0.26 and P=0.80) and is not claimed as replicated.

## Reproducibility and data boundary

Raw NHANES and MIMIC-IV tables are not committed. Download them from their source custodians and run the scripts locally. NHANES public-use materials are federal public-use data (generally public domain, subject to NCHS requirements). MIMIC-IV Demo is distributed under ODbL v1.0. Aggregate results are provided under `out/`.

The fixed drug-name-to-DrugBank mapping is included in the redistribution-safe manuscript archive.
