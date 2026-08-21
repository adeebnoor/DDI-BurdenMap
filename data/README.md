# Analytical inputs and data boundary

The two fixed redistributable network inputs are supplied in the manuscript's redistribution-safe archive:

| Expected path | Purpose | SHA-256 |
|---|---|---|
| `data/d3_candidate_ddi_pairs.csv` | Fixed filtered candidate-pair input | `a697fb1cc6a0731c16c6ec6c11a7e1327405a64f9c7d602ef83f87e24187ceec` |
| `data/GoldD3R.txt` | Construction-reference / mechanism-labelled input | `c6c189ac407ff0c9b9aed6de53dc7b33f2b240522fefdbfebcf198e3a0fccca9` |

NHANES, MIMIC-IV, DDInter, and ONC source files are **not redistributed** in this repository. Obtain them from the original providers and comply with their terms. Only aggregate cohort results are stored under `out/`.

See `../DATA_LICENSE_NOTICE.md` for the full boundary.
