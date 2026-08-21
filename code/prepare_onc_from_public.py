#!/usr/bin/env python3
"""Fetch the ONC expert-consensus DDI reference lists from their public source.

The two ONC lists analysed in this study are distributed as drug-level expanded,
DrugBank-mapped CSVs inside the public potential-DDI reference collection
maintained by the University of Pittsburgh DBMI group. That repository carries no
declared licence, so its files are NOT redistributed with this article or its
reproducibility archive. This script retrieves them directly from the source at a
pinned commit, so the analysis remains exactly reproducible without our
redistributing third-party material.

Underlying consensus lists (the intellectual content) are published as:
  Phansalkar S, et al. High-priority drug-drug interactions for use in electronic
    health records. J Am Med Inform Assoc. 2012;19(5):735-743.
  Phansalkar S, et al. Drug-drug interactions that should be non-interruptive in
    order to reduce alert fatigue in electronic health records.
    J Am Med Inform Assoc. 2013;20(3):489-493.

Usage:
  python prepare_onc_from_public.py [--out-dir PDDI-Datasets] [--commit SHA]

Then:
  python onc_clinical_relevance_analysis.py \
      data/d3_candidate_ddi_pairs.csv  PDDI-Datasets  out/onc_results.json
"""
from __future__ import annotations

import argparse
import hashlib
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = "https://github.com/dbmi-pitt/public-PDDI-analysis.git"
PINNED_COMMIT = "8199ee66b60bcb337f777889a210dd0d72a96e8f"

WANTED = [
    ("PDDI-Datasets/ONC-High-Priority/ONC_High_Priority_Mapped.csv",
     "ONC-High-Priority/ONC_High_Priority_Mapped.csv"),
    ("PDDI-Datasets/ONC-Non-Interuptive/ONC_Non_Interuptive_List_Mapped_Formatted.csv",
     "ONC-Non-Interuptive/ONC_Non_Interuptive_List_Mapped_Formatted.csv"),
    ("PDDI-Datasets/ONC-Non-Interuptive/ONC_Non_Interuptive_Mapped.csv",
     "ONC-Non-Interuptive/ONC_Non_Interuptive_Mapped.csv"),
]

# The pinned git commit fixes all file identities. Additional SHA-256 checks are
# recorded for the High-Priority and class-annotated formatted Non-Interruptive
# files used by the knowledge-base analysis.
EXPECTED_SHA256 = {
    "ONC-High-Priority/ONC_High_Priority_Mapped.csv":
        "6478cc0ff57244ed66f84a6f534596f8d3f6267af5fb2bad869c60aa53ef37e7",
    "ONC-Non-Interuptive/ONC_Non_Interuptive_List_Mapped_Formatted.csv":
        "8d5d5702a3c8ad46b1b2d5c1e5536a3d690113cfeb1522559fdfe96818b6ee1d",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default="PDDI-Datasets")
    ap.add_argument("--commit", default=PINNED_COMMIT)
    ap.add_argument("--skip-verify", action="store_true",
                    help="do not fail if recorded checksums differ")
    a = ap.parse_args()

    out = Path(a.out_dir)
    with tempfile.TemporaryDirectory() as tmp:
        print(f"Cloning {REPO} at {a.commit[:12]} ...")
        subprocess.run(["git", "clone", "--quiet", REPO, tmp], check=True)
        subprocess.run(["git", "-C", tmp, "checkout", "--quiet", a.commit], check=True)
        for src, dst in WANTED:
            s = Path(tmp) / src
            if not s.exists():
                sys.exit(f"expected file missing upstream at this commit: {src}")
            d = out / dst
            d.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(s, d)
            got = sha256(d)
            want = EXPECTED_SHA256.get(dst)
            status = "pinned-commit identity"
            if want:
                status = "sha256 ok"
                if got != want:
                    status = "CHECKSUM MISMATCH"
                    if not a.skip_verify:
                        sys.exit(f"{dst}: expected {want}, got {got}")
            print(f"  {dst}  sha256={got[:16]}...  {status}")

    print(f"\nWrote the ONC reference lists to: {out.resolve()}")
    print("These third-party files are deliberately not redistributed with the article.")


if __name__ == "__main__":
    main()
