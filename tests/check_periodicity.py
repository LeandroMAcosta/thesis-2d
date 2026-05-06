#!/usr/bin/env python3
"""2D periodicity test: verify that running the alfa=0 / sigma_l=0
sim for a deterministic period produces snapshots with bit-identical
energy. Checks the first and last snapshot.
"""

import argparse
import re
import sys
from glob import glob
from pathlib import Path


def Et(path):
    with open(path) as f:
        line = f.readline()
    m = re.search(r"Et\s*=\s*([\d.+\-eE]+)", line)
    return float(m.group(1)) if m else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("data_dir", type=Path)
    args = ap.parse_args()

    files = sorted(glob(str(args.data_dir / "X*.dat")))
    if len(files) < 2:
        print("periodicity: need at least 2 snapshots", file=sys.stderr)
        sys.exit(2)

    Et0 = Et(files[0])
    Et1 = Et(files[-1])
    if Et0 is None or Et1 is None:
        print("periodicity: cannot parse Et", file=sys.stderr); sys.exit(2)

    print(f"=== 2D Periodicity check: t=0 → t={files[-1].split('X')[-1].split('.dat')[0]} ===")
    if Et0 == Et1:
        print(f"  \033[92m✓\033[0m Energy bit-identical — Et = {Et0:.9e} J at both snapshots")
        sys.exit(0)
    rel = abs(Et1 - Et0) / abs(Et0)
    if rel < 1e-12:
        print(f"  \033[92m✓\033[0m Energy preserved within reduction-noise — rel={rel:.2e}")
        sys.exit(0)
    print(f"  \033[91m✗\033[0m Energy diverged — rel={rel:.2e}")
    sys.exit(1)


if __name__ == "__main__":
    main()
