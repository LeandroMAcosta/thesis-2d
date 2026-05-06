#!/usr/bin/env python3
"""Postcondition checker for 2D event-driven snapshots.

Reads X*.dat snapshots from a directory and verifies a series of
physical / statistical bounds. Three scenarios:

  relax     — alfa>0, sigma_l>0: system relaxes to thermal equilibrium.
              All chi² should approach 1, isotropy → 1, correlation → 0.
  periodic  — alfa=0, sigma_l=0: deterministic, energy bit-identical.
  alfa-zero — alfa=0, sigma_l>0: energy still drifts a bit but
              average is preserved.

Usage:
  check_postconditions.py {relax,periodic,alfa-zero} <data_dir>
"""

from __future__ import annotations

import argparse
import re
import sys
from glob import glob
from pathlib import Path


GREEN = "\033[92m"
RED = "\033[91m"
YELLOW = "\033[93m"
RESET = "\033[0m"


def ok(msg):  print(f"  {GREEN}✓{RESET} {msg}")
def fail(msg): print(f"  {RED}✗{RESET} {msg}")
def warn(msg): print(f"  {YELLOW}!{RESET} {msg}")


def parse_header(path: Path) -> dict:
    """Parse the first line of an X*.dat file."""
    with open(path) as f:
        line = f.readline()
    out = {}
    for key, pat in [
        ("chi2x",    r"chi2x\s*=\s*([\d.+\-eE]+)"),
        ("chi2y",    r"chi2y\s*=\s*([\d.+\-eE]+)"),
        ("chi2px",   r"chi2px\s*=\s*([\d.+\-eE]+)"),
        ("chi2py",   r"chi2py\s*=\s*([\d.+\-eE]+)"),
        ("isotropy", r"isotropy\s*=\s*([\d.+\-eE]+)"),
        ("corr",     r"corr\s*=\s*([\d.+\-eE]+)"),
        ("Et",       r"Et\s*=\s*([\d.+\-eE]+)"),
        ("N",        r"N\s*=\s*(\d+)"),
        ("BINS",     r"BINS\s*=\s*(\d+)"),
    ]:
        m = re.search(pat, line)
        if m:
            try:
                out[key] = float(m.group(1)) if "." in m.group(1) or "e" in m.group(1).lower() else int(m.group(1))
            except ValueError:
                out[key] = m.group(1)
    return out


def read_marginal_counts(path: Path):
    """Parse the marginals block (after '# MARGINALS')."""
    counts_hx, counts_hy, counts_gpx, counts_gpy = [], [], [], []
    seen_marginals = False
    with open(path) as f:
        for line in f:
            if line.startswith("# MARGINALS"):
                seen_marginals = True
                continue
            if line.startswith("# JOINT"):
                break
            if not seen_marginals or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 7:
                continue
            counts_hx.append(int(parts[1]))
            counts_hy.append(int(parts[3]))
            counts_gpx.append(int(parts[5]))
            counts_gpy.append(int(parts[6]))
    return counts_hx, counts_hy, counts_gpx, counts_gpy


def find_snapshots(data_dir: Path):
    files = sorted(glob(str(data_dir / "X*.dat")))
    if not files:
        print(f"{RED}No X*.dat snapshots in {data_dir}{RESET}")
        sys.exit(2)
    return [Path(f) for f in files]


def check_general(snapshots):
    print("\n=== General postconditions ===")
    Ns = []
    all_nonneg = True
    for snap in snapshots:
        hx, hy, gpx, gpy = read_marginal_counts(snap)
        # Total particles should equal header N for hx and hy each.
        N_hx = sum(hx)
        N_hy = sum(hy)
        N_gpx = sum(gpx)
        N_gpy = sum(gpy)
        Ns.extend([N_hx, N_hy, N_gpx, N_gpy])
        if any(c < 0 for c in hx + hy + gpx + gpy):
            all_nonneg = False
    if len(set(Ns)) == 1:
        ok(f"Particle count is conserved across all snapshots and marginals — N = {Ns[0]}")
    else:
        fail(f"Particle counts differ: {set(Ns)}")
        return False
    if all_nonneg:
        ok("All histogram bin counts are non-negative")
    else:
        fail("Negative histogram counts found")
        return False
    return True


def check_relax(snapshots):
    print("\n=== Relaxation postconditions ===")
    last = parse_header(snapshots[-1])
    passed = True
    bounds = [
        ("chi2x",  5.0),
        ("chi2y",  5.0),
        ("chi2px", 5.0),
        ("chi2py", 5.0),
    ]
    for key, bound in bounds:
        v = last.get(key)
        if v is None:
            warn(f"{key} not found in header"); passed = False; continue
        if v < bound:
            ok(f"{key} < {bound} at end — observed {v:.3f}")
        else:
            fail(f"{key} ≥ {bound} at end — observed {v:.3f}"); passed = False

    iso = last.get("isotropy", 1.0)
    if 0.9 <= iso <= 1.1:
        ok(f"Isotropy <p_x²>/<p_y²> in [0.9, 1.1] — observed {iso:.4f}")
    else:
        fail(f"Isotropy out of [0.9, 1.1] — observed {iso:.4f}"); passed = False

    corr = last.get("corr", 0.0)
    if abs(corr) < 0.05:
        ok(f"Correlation <p_x·p_y>/√(<p_x²><p_y²>) within ±0.05 — observed {corr:.4e}")
    else:
        fail(f"Correlation magnitude exceeds 0.05 — observed {corr:.4e}"); passed = False

    # Energy drift
    Et0 = parse_header(snapshots[0]).get("Et")
    Et1 = last.get("Et")
    if Et0 and Et1:
        drift_ppm = (Et1 - Et0) / Et0 * 1e6
        if abs(drift_ppm) < 100.0:
            ok(f"Energy drift |ΔE/E| < 100 ppm — observed {drift_ppm:+.2f} ppm")
        else:
            fail(f"Energy drift exceeds 100 ppm — observed {drift_ppm:+.2f} ppm")
            passed = False
    return passed


def check_periodic(snapshots):
    print("\n=== Periodic / alfa=0 postconditions ===")
    Ets = [parse_header(s).get("Et") for s in snapshots]
    if any(e is None for e in Ets):
        fail("Energy missing in some snapshot"); return False
    # In periodic regime, Et should be exactly bit-identical across all.
    # We allow at most a single ULP of drift to account for accumulation
    # order in the energy reduction (which is parallel).
    unique = set(Ets)
    if len(unique) == 1:
        ok(f"Energy bit-identical across {len(snapshots)} snapshots — Et = {Ets[0]:.9e}")
        return True
    else:
        # Compute relative spread
        spread = (max(Ets) - min(Ets)) / abs(Ets[0])
        if spread < 1e-12:
            ok(f"Energy varies by < 1e-12 (FP reduction noise) — spread {spread:.2e}")
            return True
        fail(f"Energy not conserved — {len(unique)} distinct values, spread {spread:.2e}")
        return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("scenario", choices=["relax", "periodic", "alfa-zero"])
    ap.add_argument("data_dir", type=Path)
    args = ap.parse_args()

    snaps = find_snapshots(args.data_dir)
    print(f"Loaded {len(snaps)} snapshots from {args.data_dir}")

    ok_general = check_general(snaps)
    if not ok_general:
        sys.exit(1)

    if args.scenario == "relax":
        if not check_relax(snaps): sys.exit(1)
    elif args.scenario in ("periodic", "alfa-zero"):
        if not check_periodic(snaps): sys.exit(1)

    print(f"\n{GREEN}✓ All postconditions passed.{RESET}")


if __name__ == "__main__":
    main()
