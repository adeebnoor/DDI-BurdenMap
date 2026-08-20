#!/usr/bin/env python3
"""Reconstruct the DDInter analytical edge list from DDInter's public category downloads.

The DDInter source CSVs are NOT redistributed with this archive because DDInter states
that its data are licensed CC BY-NC-SA 4.0. This script downloads the eight files from
the public DDInter download endpoint at runtime and creates the de-duplicated local
analysis file expected by ddinter_severity_analysis.py.

Source page: https://ddinter.scbdd.com/download/
Terms:       https://ddinter.scbdd.com/terms/

Access date used in the manuscript: 2026-08-16.
"""
from __future__ import annotations

import argparse
from pathlib import Path
from urllib.request import Request, urlopen
import shutil
import pandas as pd

BASE = "https://ddinter.scbdd.com/static/media/download"
CODES = ("A", "B", "D", "H", "L", "P", "R", "V")
URLS = {c: f"{BASE}/ddinter_downloads_code_{c}.csv" for c in CODES}
REQUIRED = ["DDInterID_A", "Drug_A", "DDInterID_B", "Drug_B", "Level"]
SEVERITY_RANK = {"Unknown": 0, "Minor": 1, "Moderate": 2, "Major": 3}


def fetch(url: str, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    req = Request(url, headers={"User-Agent": "DDI-BurdenMap/1.0"})
    with urlopen(req, timeout=120) as r, path.open("wb") as f:
        shutil.copyfileobj(r, f)


def load_category(code: str, raw_dir: Path, download: bool) -> pd.DataFrame:
    path = raw_dir / f"ddinter_downloads_code_{code}.csv"
    if download:
        print(f"Downloading {URLS[code]}")
        fetch(URLS[code], path)
    if not path.exists():
        raise FileNotFoundError(
            f"Missing {path}. Re-run with --download or place the public DDInter file there."
        )
    df = pd.read_csv(path, dtype=str)
    missing = [c for c in REQUIRED if c not in df.columns]
    if missing:
        raise ValueError(f"{path.name}: missing required columns {missing}; found {list(df.columns)}")
    df = df[REQUIRED].copy()
    df["atc_cat"] = code
    return df


def canonicalize(df: pd.DataFrame) -> pd.DataFrame:
    df = df.dropna(subset=["DDInterID_A", "DDInterID_B"]).copy()
    df["Level"] = df["Level"].fillna("Unknown").replace("", "Unknown")
    unknown = sorted(set(df["Level"]) - set(SEVERITY_RANK))
    if unknown:
        print("Warning: unrecognized severity labels mapped to Unknown:", unknown)
        df.loc[df["Level"].isin(unknown), "Level"] = "Unknown"

    swap = df["DDInterID_A"] > df["DDInterID_B"]
    for ca, cb in (("DDInterID_A", "DDInterID_B"), ("Drug_A", "Drug_B")):
        a = df.loc[swap, ca].copy()
        df.loc[swap, ca] = df.loc[swap, cb].values
        df.loc[swap, cb] = a.values

    df = df[df["DDInterID_A"] != df["DDInterID_B"]].copy()
    df["_sev"] = df["Level"].map(SEVERITY_RANK).fillna(0).astype(int)
    # Stable category order means ties retain the first public export in A/B/D/H/L/P/R/V order.
    df["_cat_order"] = pd.Categorical(df["atc_cat"], categories=list(CODES), ordered=True).codes
    df = df.sort_values(
        ["DDInterID_A", "DDInterID_B", "_sev", "_cat_order"],
        ascending=[True, True, False, True],
        kind="stable",
    )
    out = df.drop_duplicates(["DDInterID_A", "DDInterID_B"], keep="first").copy()
    return out[REQUIRED + ["atc_cat"]].sort_values(
        ["DDInterID_A", "DDInterID_B"], kind="stable"
    ).reset_index(drop=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-dir", default="data/raw_ddinter", help="directory for the 8 source CSVs")
    ap.add_argument("--output", default="data/ddinter2_unique.csv", help="local reconstructed edge list")
    ap.add_argument("--download", action="store_true", help="download/overwrite the eight public files")
    args = ap.parse_args()

    raw_dir = Path(args.raw_dir)
    frames = [load_category(c, raw_dir, args.download) for c in CODES]
    raw = pd.concat(frames, ignore_index=True)
    out = canonicalize(raw)
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, index=False)

    print(f"raw rows: {len(raw):,}")
    print(f"unique unordered pairs: {len(out):,}")
    print(f"unique drugs: {len(set(out.DDInterID_A) | set(out.DDInterID_B)):,}")
    print("severity counts:", out["Level"].value_counts(dropna=False).to_dict())
    print(f"wrote local derivative: {out_path}")
    print("Do not redistribute this derivative outside DDInter's license terms.")


if __name__ == "__main__":
    main()
