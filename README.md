<div align="center">

# DDI-BurdenMap

### Hub structure in drug-drug interaction burden — validated in real medication use

**Mechanism-, severity-, expert-consensus-, and patient-cohort-aware analysis.**

[![Python 3.12](https://img.shields.io/badge/Python-3.12-3776AB?logo=python&logoColor=white)](requirements.txt)
[![Primary patient validation](https://img.shields.io/badge/primary%20validation-NHANES-0E7490)](#patient-cohort-validation)
[![Inpatient sensitivity](https://img.shields.io/badge/inpatient%20sensitivity-MIMIC--IV%20Demo-6B7280)](#patient-cohort-validation)
[![DDInter boundary](https://img.shields.io/badge/DDInter-CC%20BY--NC--SA%204.0-important)](DATA_LICENSE_NOTICE.md)

*Reproducibility repository for:*  
**Hub drugs concentrate drug–drug interaction burden across curated networks and real prescribing: a mechanism-aware network analysis with patient-cohort validation**

</div>

---

## What this study tests

Most DDI work is pair-centered. This study asks a complementary systems question:

> **Is candidate DDI burden concentrated on a manageable minority of drugs, and does that fixed drug-level structure transport into real medication use?**

The hub ranking is created **before** patient data are introduced. Patient cohorts are validation/transport layers, not inputs to the ranking.

## Evidence stack

| Evidence layer | Main result | Role |
|---|---|---|
| Filtered candidate network | 1,114 drugs; 16,316 pairs; top 10% carry 42.3% of endpoints; Gini 0.637 | Discovery |
| Construction-reference network | top 10% carry 45.1%; Gini 0.654 | Internal consistency |
| Independent DDInter export | 160,235 pairs; top 10% carry 31.6%; Gini 0.503 | Independent structural replication |
| ONC expert-consensus lists | 77.4% Non-Interruptive vs 39.5% High-Priority coverage at top decile | Clinical-relevance benchmark |
| NHANES 2015–2018 | 7,669 medication users; 70.0% of co-taken candidate-network pairs covered vs 18.9% random null | **Primary independent patient validation** |
| MIMIC-IV Demo v2.2 | 100 patients / 250 admissions; 68.0% of temporally overlapping candidate-network pairs covered vs 18.9% null | Secondary inpatient transport/sensitivity |

**Interpretation boundary:** degree is a workload/connectivity and prioritization signal. It is **not** a confirmed-interaction label, patient-level risk score, stand-alone severity score, fired-alert measure, or validated suppression rule.

<p align="center">
  <img src="assets/concentration_curve.svg" width="820" alt="Concentration summary across candidate, reference and DDInter networks">
</p>

## Patient-cohort validation

### NHANES 2015–2018 — primary validation

`code/nhanes_cohort_validation.py` evaluates the fixed watchlist in CDC/NCHS prescription-medication data. The analytical cohort contains **7,669 participants with at least one mapped prescription**, including 5,301 using two or more mapped drugs.

At the 10% watchlist, **70.0% of co-taken candidate-network pairs** are covered versus an **18.9%** random-watchlist null (10,000 draws; seed 42; empirical P < 1e-4).

### MIMIC-IV Demo v2.2 — secondary inpatient sensitivity analysis

`code/mimic_cohort_validation.py` requires temporal overlap of prescription windows within admission. The open demo contains **100 de-identified patients and 250 admissions**.

At the 10% watchlist, **68.0% of temporally overlapping candidate-network pairs** are covered versus an **18.9%** random-watchlist null (empirical P < 1e-4).

Because the demo is intentionally small, the manuscript treats it as a **transport/sensitivity analysis**, not as a second large definitive clinical cohort.

See [`code/README_cohort_validation.md`](code/README_cohort_validation.md).

## ONC expert-consensus analyses

Two ONC analyses use two mapped representations from the same pinned upstream PDDI repository commit:

- **Fig. 7 knowledge-base coverage:** class-annotated formatted Non-Interruptive representation → 1,895 candidate-network pairs after correction/restriction.
- **Fig. 8b patient realization:** fully expanded mapped Non-Interruptive representation → 2,025 candidate-network pairs after the corresponding correction/restriction.

These denominators are analysis-specific and are not compared directly. `code/prepare_onc_from_public.py` retrieves both representations plus the High-Priority list from pinned upstream commit `8199ee66b60bcb337f777889a210dd0d72a96e8f`.

Across patient cohorts, Non-Interruptive pairs are co-prescribed **6.9×** more often than High-Priority pairs in NHANES and **3.7×** more often in MIMIC-IV Demo. The **watchlist-specific ONC contrast among realized pairs was underpowered/non-significant** (P=0.26 and P=0.80) and is not claimed as replicated.

## Independent DDInter replication

DDInter source CSVs and the processed pair-level derivative are **not committed** because DDInter states a **CC BY-NC-SA 4.0** license.

```bash
python code/prepare_ddinter_from_public.py --download \
  --raw-dir data/raw_ddinter --output data/ddinter2_unique.csv

python code/ddinter_severity_analysis.py \
  data/ddinter2_unique.csv out/ddinter_severity_results.json
```

See [`DATA_LICENSE_NOTICE.md`](DATA_LICENSE_NOTICE.md).

## Reproducibility

```bash
git clone https://github.com/adeebnoor/DDI-BurdenMap.git
cd DDI-BurdenMap
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

The submission reproducibility archive contains the two redistributable fixed network inputs and the fixed name-to-DrugBank mapping used by the patient analyses. Raw NHANES, MIMIC-IV, DDInter and ONC mapped source files are obtained from their custodians rather than re-hosted here.

Aggregate result JSONs are mirrored in `out/`; cohort and ONC scripts are under `code/`.

## Reproducibility design

- Candidate-network ranking is fixed before patient data are introduced.
- Random procedures use fixed seed 42 where applicable.
- Expert-label permutation and random-watchlist nulls are explicit.
- MIMIC-IV Demo is labeled as a sensitivity/transport analysis.
- The negative patient-level ONC watchlist contrast is retained rather than hidden.
- Candidate-network co-exposure is not equated with a confirmed DDI or an actual clinical alert.

## Citation

Until the journal article receives its final bibliographic citation, use [`CITATION.cff`](CITATION.cff). No Zenodo DOI is claimed unless one is actually minted.

## Author

**Adeeb Noor**  
Department of Information Technology, Faculty of Computing and Information Technology  
King Abdulaziz University, Jeddah, Saudi Arabia

---

<div align="center">

**DDI-BurdenMap — from pair lists to an auditable, patient-grounded map of interaction burden.**

</div>
