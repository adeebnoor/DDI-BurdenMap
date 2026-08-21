# Third-party data and redistribution notice

DDI-BurdenMap deliberately separates **open analysis code and aggregate results** from third-party source data and patient-level/pair-level derivatives.

## DDInter

DDInter states that its data are available under **CC BY-NC-SA 4.0**. Its source CSVs and the processed `ddinter2_unique.csv` pair-level derivative are not redistributed here. `code/prepare_ddinter_from_public.py` reconstructs the analytical edge list locally from DDInter's public category downloads. Users remain responsible for DDInter's terms.

- Terms: https://ddinter.scbdd.com/terms/
- Downloads: https://ddinter.scbdd.com/download/

## NHANES 2015-2018

NHANES files used for the human validation are **NCHS public-use files**, not material re-licensed by this repository. Users should obtain RXQ_RX_I, RXQ_RX_J and, if needed, DEMO_I/DEMO_J directly from CDC/NCHS and comply with the NCHS Data User Agreement:

- https://www.cdc.gov/nchs/policy/data-user-agreement.html
- https://wwwn.cdc.gov/nchs/nhanes/

The repository does not contain NHANES source records or participant-level/pair-level derived data. `out/nhanes_results.json` contains aggregate analytical results only.

## MIMIC-IV Clinical Database Demo v2.2

The MIMIC-IV Demo is distributed by PhysioNet under the **Open Data Commons Open Database License v1.0 (ODbL)**. Obtain the source tables directly from PhysioNet (DOI: 10.13026/dp1f-ex47). This repository does not redistribute MIMIC source tables or patient/pair-level derivatives. `out/mimic_results.json` contains aggregate results only.

## ONC / public-PDDI-analysis

The expert-consensus analysis reads the public ONC-derived mapped files from the `dbmi-pitt/public-PDDI-analysis` resource. Those source files are not copied into this repository. The manuscript reports aggregate coverage/permutation results produced by `code/onc_clinical_relevance_analysis.py`.

## No blanket relicensing

Nothing in this repository should be read as relicensing third-party databases. Only repository-authored code, documentation, and aggregate outputs are supplied here, subject to the notices and upstream terms applicable to each source.
