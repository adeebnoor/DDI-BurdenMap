<div align="center">

# DDI-BurdenMap

### Mapping where drug-drug interaction burden concentrates

**A severity- and mechanism-aware network analysis across candidate, reference, and independently curated DDI resources.**

[![Python 3.12](https://img.shields.io/badge/Python-3.12-3776AB?logo=python&logoColor=white)](requirements.txt)
[![Reproducibility](https://img.shields.io/badge/reproducibility-fixed%20seeds-2E8B57)](code/)
[![External validation](https://img.shields.io/badge/external%20validation-DDInter-0E7490)](#independent-ddinter-validation)
[![Zenodo](https://img.shields.io/badge/DOI-Zenodo%20release%20pending-6F42C1)](.zenodo.json)
[![DDInter boundary](https://img.shields.io/badge/DDInter-CC%20BY--NC--SA%204.0-important)](DATA_LICENSE_NOTICE.md)

*Companion reproducibility repository for:*  
**Hub drugs concentrate drug-drug interaction burden across candidate and curated networks: a severity- and mechanism-aware network analysis**

</div>

---

## The question

Most DDI research asks whether a **pair** of drugs interacts. DDI-BurdenMap asks a complementary systems-level question:

> **Is interaction burden spread broadly across drugs, or does it concentrate around a relatively small set of hub drugs?**

The project maps that structure, tests whether it persists across datasets and pharmacological mechanisms, and separately asks whether topology enriches for severity.

<p align="center">
  <img src="assets/concentration_curve.svg" width="820" alt="Concentration summary for DDI burden across candidate, reference and DDInter networks">
</p>

## Results at a glance

| Network | Drugs | Pairs | Gini | Top 10% endpoint share | Giant component |
|---|---:|---:|---:|---:|---:|
| Filtered candidate network | 1,114 | 16,316 | **0.637** | **42.3%** | 99.5% |
| Construction-reference network | 1,599 | 21,897 | **0.654** | **45.1%** | 98.4% |
| Independent DDInter export | 1,939 | 160,235 | **0.503** | **31.6%** | 99.7% |

**Interpretation:** DDI burden is consistently hub-concentrated, but **degree is a workload/connectivity signal - not a patient-level risk score and not a stand-alone severity score**.

<p align="center">
  <img src="assets/ddinter_validation.svg" width="900" alt="Independent DDInter severity-aware validation summary">
</p>

## What is here

- Analysis scripts for candidate/reference topology, mechanism stratification, power-law testing, DDInter severity analysis, and figure generation.
- Machine-readable aggregate results under `out/`.
- A transparent licensing boundary: DDInter raw and pair-level derivative files are never redistributed here.
- `CITATION.cff` and `.zenodo.json` prepared for a versioned archival release.
- Fixed random seeds and pinned dependencies for reproducible execution.

## Study workflow

<p align="center">
  <img src="assets/workflow.svg" width="900" alt="Analytical workflow from fixed candidate DDI pairs to topology and independent validation">
</p>

The candidate set is a **fixed analytical input**. The construction-reference network is an **internal consistency / mechanism-stratified check** because it shares upstream provenance. DDInter provides the **independent external structural replication**.

## Reproduce the analyses

```bash
git clone https://github.com/adeebnoor/DDI-BurdenMap.git
cd DDI-BurdenMap

python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

The two redistributable fixed inputs used for the candidate/reference analyses are distributed in the manuscript's redistribution-safe reproducibility archive. Place them at:

```text
data/d3_candidate_ddi_pairs.csv
data/GoldD3R.txt
```

Their expected SHA-256 digests are documented in [`data/README.md`](data/README.md). Then run:

```bash
python code/candidate_network_analysis.py data/d3_candidate_ddi_pairs.csv
python code/confirmed_network_analysis.py data/GoldD3R.txt
python code/mechanism_stratified_analysis.py \
  data/GoldD3R.txt out/mechanism_topology.csv out/mechanism_topology.json
python code/reference_network_powerlaw.py data/GoldD3R.txt
```

## Independent DDInter validation

DDInter source CSVs and the processed pair-level derivative are **not committed** because DDInter states a **CC BY-NC-SA 4.0** license. Reconstruct the analytical network locally from its public category downloads:

```bash
python code/prepare_ddinter_from_public.py \
  --download \
  --raw-dir data/raw_ddinter \
  --output data/ddinter2_unique.csv

python code/ddinter_severity_analysis.py \
  data/ddinter2_unique.csv \
  out/ddinter_severity_results.json
```

`data/raw_ddinter/` and `data/ddinter2_unique.csv` are blocked by `.gitignore` to reduce accidental redistribution. See [`DATA_LICENSE_NOTICE.md`](DATA_LICENSE_NOTICE.md).

## Repository map

```text
DDI-BurdenMap/
├── assets/              # project-facing scientific visuals
├── code/                # analysis + figure-generation scripts
├── data/README.md       # input provenance, paths, SHA-256 digests
├── out/                 # aggregate machine-readable outputs
├── CITATION.cff
├── .zenodo.json
├── requirements.txt
└── DATA_LICENSE_NOTICE.md
```

## Reproducibility design

Random procedures use fixed seeds recorded in the scripts. The candidate-network ranking uses pair membership only; pharmacological labels do not enter hub ranking. Mechanism-stratified analyses are reported separately, and DDInter severity is tested against explicit permutation nulls.

## Data and license boundary

This repository does **not** relicense third-party databases. In particular, DDInter source files and DDInter-derived pair-level data are not redistributed. The redistribution-safe archive supplied with the manuscript contains only material that is appropriate to share and is prepared for permanent archival.

## Citation

Until the journal article receives its final bibliographic citation, use [`CITATION.cff`](CITATION.cff). The repository metadata are prepared for a `v1.0.0` archival release; a Zenodo DOI can be added after that release is archived.

## Author

**Adeeb Noor**  
Department of Information Technology, Faculty of Computing and Information Technology  
King Abdulaziz University, Jeddah, Saudi Arabia

---

<div align="center">

**DDI-BurdenMap - from pair lists to an auditable map of interaction burden.**

</div>
