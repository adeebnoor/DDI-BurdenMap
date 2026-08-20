# Redistributable analytical inputs

The two fixed inputs used for the candidate and construction-reference analyses are stored in `data/source/` as deterministic gzip + Base64 text fragments. Reconstruct them with:

```bash
python code/bootstrap_inputs.py
```

The script verifies the reconstructed bytes against the archived SHA-256 digests before analysis.

| File | SHA-256 |
|---|---|
| `GoldD3R.txt` | `c6c189ac407ff0c9b9aed6de53dc7b33f2b240522fefdbfebcf198e3a0fccca9` |
| `d3_candidate_ddi_pairs.csv` | `a697fb1cc6a0731c16c6ec6c11a7e1327405a64f9c7d602ef83f87e24187ceec` |

DDInter is **not** stored here. See [`../DATA_LICENSE_NOTICE.md`](../DATA_LICENSE_NOTICE.md) and reconstruct it directly from DDInter with `code/prepare_ddinter_from_public.py`.
