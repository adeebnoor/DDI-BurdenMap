#!/usr/bin/env python3
"""Reconstruct the two redistributable fixed analytical inputs from compact repository fragments."""
from __future__ import annotations
import base64, gzip, hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "data" / "source"
DATA = ROOT / "data"
FILES = {
    "GoldD3R.txt": {
        "pattern": "GoldD3R_txt.gz.b64.chunk*",
        "sha256": "c6c189ac407ff0c9b9aed6de53dc7b33f2b240522fefdbfebcf198e3a0fccca9",
    },
    "d3_candidate_ddi_pairs.csv": {
        "pattern": "d3_candidate_ddi_pairs_csv.gz.b64.chunk*",
        "sha256": "a697fb1cc6a0731c16c6ec6c11a7e1327405a64f9c7d602ef83f87e24187ceec",
    },
}

for name, spec in FILES.items():
    parts = sorted(SRC.glob(spec["pattern"]))
    if not parts:
        raise FileNotFoundError(f"No source fragments found for {name}")
    encoded = "".join(p.read_text(encoding="ascii").strip() for p in parts)
    raw = gzip.decompress(base64.b64decode(encoded))
    digest = hashlib.sha256(raw).hexdigest()
    if digest != spec["sha256"]:
        raise RuntimeError(f"SHA-256 mismatch for {name}: {digest}")
    out = DATA / name
    out.write_bytes(raw)
    print(f"reconstructed {out.relative_to(ROOT)} ({len(raw):,} bytes; sha256={digest})")
