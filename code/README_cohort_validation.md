# Patient-cohort validation

This repository now includes two independent real-world medication-use validations added in response to the Frontiers editorial request for patient/clinical validation.

## NHANES 2015-2018 (ambulatory)

Public CDC/NCHS prescription-medication files:
- RXQ_RX_I.XPT (2015-2016)
- RXQ_RX_J.XPT (2017-2018)
- DEMO_I.XPT and DEMO_J.XPT for cohort description

Run `nhanes_cohort_validation.py` with the fixed candidate network and the documented drug-name mapping. The analytical cohort contains 7,669 participants with at least one mapped prescription; 5,301 used two or more mapped drugs. At the 10% watchlist, 70.0% of candidate-alertable co-taken pairs were covered versus a 18.9% random-watchlist null (10,000 draws; seed 42; empirical P < 1e-4).

## MIMIC-IV Clinical Database Demo v2.2 (inpatient sensitivity cohort)

PhysioNet DOI: 10.13026/dp1f-ex47. The open demo contains 100 de-identified patients and 250 admissions. Prescriptions are evaluated using true temporal overlap of start/stop windows within admission; a same-admission sensitivity definition is also computed.

At the 10% watchlist, 68.0% of candidate-alertable temporally overlapping pairs were covered versus a 18.9% random-watchlist null (10,000 draws; seed 42; empirical P < 1e-4).

Because the MIMIC-IV demo is intentionally small, the manuscript treats it as an inpatient transport/sensitivity cohort rather than as a large definitive clinical-effectiveness cohort. The primary patient-level validation is NHANES.

## Reproducibility boundary

Raw NHANES and MIMIC-IV files are not committed. Download them from the source custodians and run the scripts locally. Aggregate results are provided under `out/`. The fixed mapping file used by the submitted analysis is included in the redistribution-safe manuscript archive.
