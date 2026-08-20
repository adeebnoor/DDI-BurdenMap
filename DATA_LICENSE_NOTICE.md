# Third-party data license notice — DDInter

DDInter states that its data are made available under **Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)**:

- Terms: https://ddinter.scbdd.com/terms/
- Public downloads: https://ddinter.scbdd.com/download/

Frontiers publishes article content under **CC BY 4.0**. To avoid redistributing a DDInter-derived pair-level dataset under an incompatible article/supplement license, this repository **does not contain** the DDInter source CSVs or the processed `ddinter2_unique.csv` derivative.

Instead, `code/prepare_ddinter_from_public.py` downloads the eight public ATC-category files (A, B, D, H, L, P, R, V) directly from DDInter and reconstructs the local de-duplicated analytical edge list. The local file is then consumed by `code/ddinter_severity_analysis.py`. Users remain responsible for complying with DDInter's license and terms.

Only aggregate analytical outputs needed to audit the manuscript are included here; no new license over DDInter data is asserted.
