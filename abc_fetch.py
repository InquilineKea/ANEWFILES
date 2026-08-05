#!/usr/bin/env python3
"""
Download Allen Brain Cell Atlas metadata directly from its public S3 bucket.

The ABC atlas bucket is anonymously readable over plain HTTPS, so this needs no
AWS credentials and no boto3 -- just urllib. That matters in locked-down
environments, where the ABC bucket is often reachable when the usual portal
hosts (10x, HuBMAP, figshare, cellxgene's release directory) are not.

    python3 abc_fetch.py --list
    python3 abc_fetch.py --dataset merfish        # 1.6 GB, has spatial coords
    python3 abc_fetch.py --dataset xenium-marmoset
    python3 abc_fetch.py --dataset aging-mouse

Then:
    python3 abc_frequencies.py merfish_cells.csv
"""

from __future__ import annotations

import argparse
import os
import sys
import urllib.request

BUCKET = "https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com"

# Curated set. `spatial` records whether the dataset carries per-cell tissue
# coordinates; `donors` is what actually constrains age-related questions.
DATASETS: dict[str, dict] = {
    "merfish": {
        "key": "metadata/MERFISH-C57BL6J-638850-CCF/20231215/views/"
               "cell_metadata_with_parcellation_annotation.csv",
        "out": "merfish_cells.csv",
        "mb": 1607,
        "spatial": True,
        "donors": 1,
        "desc": "Whole mouse brain MERFISH, 3.74M cells, CCF parcellation "
                "(gives you HPF/CA1/CA3/DG) + x/y/z coordinates. Single animal.",
    },
    "xenium-marmoset": {
        "key": "metadata/HMBA-Xenium-CJ23.56.004-BG/20250630/cell_metadata.csv",
        "out": "xenium_marmoset_cells.csv",
        "mb": 297,
        "spatial": True,
        "donors": 1,
        "desc": "Xenium, marmoset basal ganglia, one donor aged 5.7 yr. "
                "Coordinates are x/y/experiment space. No hippocampus, no age range.",
    },
    "xenium-marmoset-anno": {
        "key": "metadata/HMBA-Xenium-CJ23.56.004-BG/20250630/"
               "cell_anatomical_annotations.csv",
        "out": "xenium_marmoset_anatomy.csv",
        "mb": 116,
        "spatial": False,
        "donors": 1,
        "desc": "Anatomical labels to join onto xenium-marmoset by cell_label.",
    },
    "aging-mouse": {
        "key": "metadata/Zeng-Aging-Mouse-10Xv3/20250131/cell_metadata.csv",
        "out": "aging_mouse_cells.csv",
        "mb": 419,
        "spatial": False,
        "donors": 108,
        "desc": "Aging mouse brain 10Xv3, 108 donors in two age groups "
                "(~P56 adult vs ~P550 aged). Real age contrast, NO coordinates.",
    },
    "aging-mouse-donor": {
        "key": "metadata/Zeng-Aging-Mouse-10Xv3/20250131/donor.csv",
        "out": "aging_mouse_donors.csv",
        "mb": 1,
        "spatial": False,
        "donors": 108,
        "desc": "Donor table for aging-mouse: age, sex, genotype. Tiny; start here.",
    },
}


def human(n: int) -> str:
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024:
            return "%.1f %s" % (n, unit)
        n /= 1024.0
    return "%.1f TB" % n


def download(key: str, out: str, expect_mb: float) -> None:
    url = "%s/%s" % (BUCKET, key)
    if os.path.exists(out):
        print("already present: %s (%s) -- delete to re-fetch"
              % (out, human(os.path.getsize(out))))
        return
    print("GET %s" % url)
    print("    -> %s (~%d MB)" % (out, expect_mb))
    tmp = out + ".part"
    last = [0]

    def hook(blocks: int, bs: int, total: int) -> None:
        got = blocks * bs
        pct = 100.0 * got / total if total > 0 else 0
        if got - last[0] > 50 * 1024 * 1024:
            last[0] = got
            print("    %5.1f%%  %s" % (pct, human(got)), flush=True)

    urllib.request.urlretrieve(url, tmp, reporthook=hook)
    os.replace(tmp, out)
    print("    done: %s" % human(os.path.getsize(out)))


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--dataset", choices=sorted(DATASETS))
    p.add_argument("--list", action="store_true")
    p.add_argument("--outdir", default=".")
    args = p.parse_args()

    if args.list or not args.dataset:
        print("Available ABC atlas metadata downloads:\n")
        for name, d in DATASETS.items():
            print("  %-22s %6d MB  spatial=%-5s donors=%-4d"
                  % (name, d["mb"], d["spatial"], d["donors"]))
            print("      %s" % d["desc"])
        print("\nNote the trade-off: the datasets with coordinates have one donor;")
        print("the one with 108 donors has no coordinates.")
        return 0

    d = DATASETS[args.dataset]
    os.makedirs(args.outdir, exist_ok=True)
    download(d["key"], os.path.join(args.outdir, d["out"]), d["mb"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
