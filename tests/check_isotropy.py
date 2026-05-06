#!/usr/bin/env python3
"""2D-specific test: isotropy of momentum distribution.

After relaxation, the simulation should produce a distribution where
both momentum components have the same variance (no privileged
direction). Check: <p_x²> / <p_y²> ∈ [1−ε, 1+ε] for ε small.

Reads the header of the LAST X*.dat snapshot in a dir (header has
the precomputed isotropy from particle arrays, more accurate than
binning into histograms).
"""

import argparse
import re
import sys
from glob import glob
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("data_dir", type=Path)
    ap.add_argument("--tol", type=float, default=0.02,
                    help="absolute tolerance: |1 − iso| ≤ tol (default 0.02)")
    args = ap.parse_args()

    files = sorted(glob(str(args.data_dir / "X*.dat")))
    if not files:
        print(f"isotropy: no X*.dat in {args.data_dir}", file=sys.stderr)
        sys.exit(2)

    # Parse last snapshot header
    with open(files[-1]) as f:
        line = f.readline()
    m = re.search(r"isotropy\s*=\s*([\d.+\-eE]+)", line)
    if not m:
        print(f"isotropy: no isotropy field in {files[-1]}", file=sys.stderr)
        sys.exit(2)
    iso = float(m.group(1))

    deviation = abs(1.0 - iso)
    if deviation <= args.tol:
        print(f"\033[92m✓\033[0m isotropy <p_x²>/<p_y²> = {iso:.5f}, "
              f"|1−iso| = {deviation:.4f} ≤ {args.tol}")
        sys.exit(0)
    else:
        print(f"\033[91m✗\033[0m isotropy <p_x²>/<p_y²> = {iso:.5f}, "
              f"|1−iso| = {deviation:.4f} > {args.tol}")
        sys.exit(1)


if __name__ == "__main__":
    main()
