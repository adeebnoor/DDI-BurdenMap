# Third-party data and license notice

## DDInter

DDInter states that its data are available under **Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)**:

- Terms: https://ddinter.scbdd.com/terms/
- Public downloads: https://ddinter.scbdd.com/download/

This repository does **not** contain DDInter source CSVs or the processed pair-level `ddinter2_unique.csv` derivative. `code/prepare_ddinter_from_public.py` reconstructs the local analytical edge list from the public category downloads. Users remain responsible for DDInter's license and terms.

## ONC expert-consensus mapped reference files

The DrugBank-mapped ONC High-Priority and Non-Interruptive files are distributed in the public `dbmi-pitt/public-PDDI-analysis` repository, which declares no licence. They are therefore not re-hosted here.

`code/prepare_onc_from_public.py` retrieves all three mapped files required by the submitted analyses from pinned upstream commit `8199ee66b60bcb337f777889a210dd0d72a96e8f`:

- `ONC_High_Priority_Mapped.csv`
- `ONC_Non_Interuptive_List_Mapped_Formatted.csv` — class-annotated representation used for the Fig. 7 knowledge-base analysis
- `ONC_Non_Interuptive_Mapped.csv` — fully expanded representation used for the patient realization analysis

The pinned git commit fixes file identity; recorded SHA-256 checks additionally verify the High-Priority and class-annotated formatted Non-Interruptive files.

## NHANES

NHANES prescription and demographic public-use files are obtained from CDC/NCHS. Federal NHANES materials are generally in the public domain, subject to NCHS data-use and citation requirements. Raw source tables and patient-level derivatives are not committed here.

## MIMIC-IV Clinical Database Demo v2.2

The open MIMIC-IV Demo is distributed by PhysioNet under the **Open Data Commons Open Database License v1.0 (ODbL v1.0)**, DOI `10.13026/dp1f-ex47`. Raw source tables and patient-level derivatives are not committed here.

Only analysis code and aggregate outputs needed to audit the manuscript are public in this repository. The submission's redistribution-safe reproducibility archive contains the redistributable fixed network inputs and the fixed drug-name mapping used by the cohort scripts.
