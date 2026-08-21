<div align="center">

# DDI-BurdenMap

### Mapping where drug-drug interaction burden concentrates

**A severity-, mechanism-, and patient-cohort-aware analysis of DDI burden concentration.**

[![Python 3.12](https://img.shields.io/badge/Python-3.12-3776AB?logo=python&logoColor=white)](requirements.txt)
[![Reproducibility](https://img.shields.io/badge/reproducibility-fixed%20seeds-2E8B57)](code/)
[![Patient validation](https://img.shields.io/badge/patient%20validation-NHANES%20%2B%20MIMIC--IV-0E7490)](#patient-cohort-validation)
[![DDInter boundary](https://img.shields.io/badge/DDInter-CC%20BY--NC--SA%204.0-important)](DATA_LICENSE_NOTICE.md)

*Companion reproducibility repository for:*  
**Hub drugs concentrate drug-drug interaction burden across curated networks and real prescribing**

</div>

---

## The question

Most DDI research asks whether a **pair** of drugs interacts. DDI-BurdenMap asks a complementary systems-level question:

> **Is interaction burden spread broadly across drugs, or does it concentrate around a relatively small set of hub drugs, and does that structure transport into real medication use?**

## Results at a glance

| Evidence layer | Main result |
|---|---|
| Candidate network | 1,114 drugs; 16,316 pairs; top 10% carry 42.3% of endpoints; Gini 0.637 |
| Construction-reference network | top 10% carry 45.1% of endpoints; Gini 0.654 |
| Independent DDInter export | 160,235 pairs; top 10% carry 31.6% of endpoints; Gini 0.503 |
| NHANES 2015-2018 | 7,669 medication users; top-decile watchlist covers 70.0% of candidate-alertable co-taken pairs vs 18.9% random null |
| MIMIC-IV demo v2.2 | 100 patients / 250 admissions; top-decile watchlist covers 68.0% of candidate-alertable temporally overlapping pairs vs 18.9% random null |

**Interpretation:** degree is a transparent workload/connectivity signal. It is **not** a patient-level risk score and is **not** a stand-alone severity score.

<p align="center">
  <img src="assets/concentration_curve.svg" width="820" alt="Concentration summary for DDI burden across candidate, reference and DDInter networks">
</p>

## Patient-cohort validation

The Frontiers resubmission adds patient-level transport validation that is fully separate from the databases used to build the hub ranking.

### Ambulatory: NHANES 2015-2018

`code/nhanes_cohort_validation.py` uses the CDC/NCHS prescription-medication files for 2015-2016 and 2017-2018. The analytical cohort contains **7,669 participants with at least one mapped prescription** and 5,301 with two or more mapped drugs. At the 10% watchlist, **70.0%** of candidate-alertable co-taken pairs are covered versus a **18.9%** random-watchlist null (10,000 draws; seed 42; empirical P < 1e-4).

### Inpatient sensitivity/transport cohort: MIMIC-IV demo v2.2

`code/mimic_cohort_validation.py` evaluates true temporal overlap of prescription start/stop windows within hospitalization. The public demo contains **100 de-identified patients and 250 admissions**. At the 10% watchlist, **68.0%** of candidate-alertable temporally overlapping pairs are covered versus a **18.9%** random-watchlist null (10,000 draws; seed 42; empirical P < 1e-4).

Because the MIMIC-IV demo is intentionally small, it is treated as an **inpatient transport/sensitivity cohort**, not as a large definitive clinical-effectiveness cohort. The primary patient-level validation is NHANES. See [`code/README_cohort_validation.md`](code/README_cohort_validation.md).

## Independent DDInter validation

DDInter source CSVs and the processed pair-level derivative are **not committed** because DDInter states a **CC BY-NC-SA 4.0** license. Reconstruct the analytical network locally from its public category downloads:

```bash
python code/prepare_ddinter_from_public.py --download --raw-dir data/raw_ddinter --output data/ddinter2_unique.csv
python code/ddinter_severity_analysis.py data/ddinter2_unique.csv out/ddinter_severity_results.json
```

See [`DATA_LICENSE_NOTICE.md`](DATA_LICENSE_NOTICE.md).

## Reproduce the analyses

```bash
git clone https://github.com/adeebnoor/DDI-BurdenMap.git
cd DDI-BurdenMap
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

The two redistributable fixed inputs used for the candidate/reference analyses are supplied in the manuscript's redistribution-safe reproducibility archive. Patient-cohort raw files are downloaded from their source custodians rather than redistributed here. The fixed name-to-DrugBank mapping used for the submitted cohort analyses is included in that archive.

## Repository map

```text
DDI-BurdenMap/
├── assets/              # project-facing scientific visuals
├── code/                # topology, DDInter and patient-cohort scripts
├── data/README.md       # input provenance and licensing notes
├── out/                 # aggregate machine-readable outputs
├── CITATION.cff
├── .zenodo.json
├── requirements.txt
└── DATA_LICENSE_NOTICE.md
```

## Reproducibility design

Random procedures use fixed seeds recorded in the scripts. The candidate-network ranking uses pair membership only; pharmacological labels and patient data do not enter hub ranking. Patient cohorts are used only after the watchlist has been fixed.

## Citation

Until the journal article receives its final bibliographic citation, use [`CITATION.cff`](CITATION.cff). A Zenodo DOI can be added after a versioned release is archived; no DOI is invented here.

## Author

**Adeeb Noor**  
Department of Information Technology, Faculty of Computing and Information Technology  
King Abdulaziz University, Jeddah, Saudi Arabia

---

<div align="center">

**DDI-BurdenMap - from pair lists to an auditable map of interaction burden.**

</div>
