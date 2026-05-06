#!/usr/bin/env python3
"""2D-specific test: independence of p_x and p_y.

The two momentum components are statistically independent at
equilibrium (each is its own thermal Gaussian). The Pearson
correlation coefficient

   corr = <p_x · p_y> / √(<p_x²> <p_y²>)

should be ≈ 0.

For finite N, the standard error of corr under H_0 is ~1/√N. With
N = 2^18 = 262144, the 1-sigma noise is ~2e-3. We use 5e-2 as a
generous bound.
"""

import argparse
import re
import sys
from glob import glob
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("data_dir", type=Path)
    ap.add_argument("--tol", type=float, default=0.05,
                    help="bound on |corr| (default 0.05)")
    args = ap.parse_args()

    files = sorted(glob(str(args.data_dir / "X*.dat")))
    if not files:
        print(f"independence: no X*.dat in {args.data_dir}", file=sys.stderr)
        sys.exit(2)

    with open(files[-1]) as f:
        line = f.readline()
    m = re.search(r"corr\s*=\s*([\d.+\-eE]+)", line)
    if not m:
        print(f"independence: no corr field in {files[-1]}", file=sys.stderr)
        sys.exit(2)
    corr = float(m.group(1))

    if abs(corr) <= args.tol:
        print(f"\033[92m✓\033[0m correlation <p_x·p_y>/√(<p_x²><p_y²>) "
              f"= {corr:+.4e}, |corr| ≤ {args.tol}")
        sys.exit(0)
    else:
        print(f"\033[91m✗\033[0m correlation {corr:+.4e} exceeds tolerance "
              f"{args.tol}")
        sys.exit(1)


if __name__ == "__main__":
    main()
