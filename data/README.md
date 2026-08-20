# Analytical inputs

DDI-BurdenMap uses two fixed, redistributable inputs for the candidate and construction-reference analyses. The release/archive should contain these files at the paths below:

| Expected path | Purpose | SHA-256 |
|---|---|---|
| `data/d3_candidate_ddi_pairs.csv` | Fixed filtered candidate-pair input | `a697fb1cc6a0731c16c6ec6c11a7e1327405a64f9c7d602ef83f87e24187ceec` |
| `data/GoldD3R.txt` | Construction-reference / mechanism-labelled input | `c6c189ac407ff0c9b9aed6de53dc7b33f2b240522fefdbfebcf198e3a0fccca9` |

Verify them before analysis:

```bash
sha256sum data/d3_candidate_ddi_pairs.csv data/GoldD3R.txt
```

The analysis scripts fail visibly if the expected files are absent; no synthetic replacement is used.

## DDInter is different

DDInter is an independent third-party source and is **not redistributed in this repository**. In particular, do not commit either:

- `data/raw_ddinter/`
- `data/ddinter2_unique.csv`

Instead reconstruct the DDInter analytical edge list locally from DDInter's public downloads:

```bash
python code/prepare_ddinter_from_public.py \
  --download \
  --raw-dir data/raw_ddinter \
  --output data/ddinter2_unique.csv
```

Then run the severity-aware independent validation:

```bash
python code/ddinter_severity_analysis.py \
  data/ddinter2_unique.csv \
  out/ddinter_severity_results.json \
  figures/Fig6_ddinter_severity_corrected.png
```

See [`../DATA_LICENSE_NOTICE.md`](../DATA_LICENSE_NOTICE.md) for the DDInter licensing boundary.
