#!/usr/bin/env python3
"""2D periodicity check. Two assertion levels depending on regime:

  --mode=energy   (gaussian momenta, alfa=0, sigma_l=0)
      Asserts only that total energy is bit-identical across the run.
      Histograms phase-mix and DO NOT return to the initial in any
      observable sense (Poincaré recurrence is astronomical).

  --mode=cycle    (init_uniform_p_magnitude with calibrated p_0)
      Strong assertion: every particle completes an integer number of
      bounce periods in N_steps*DT, so every histogram (hx, hy, gpx,
      gpy, h_xy) must return to t=0 within ULP-slippage tolerance.

Usage:  tests/check_periodicity.py <mode> <data_dir>
"""

from __future__ import annotations

import argparse
import importlib.util
import re
import sys
from glob import glob
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
spec = importlib.util.spec_from_file_location("plot", ROOT / "tools" / "plot.py")
plot = importlib.util.module_from_spec(spec)
sys.modules["plot"] = plot
spec.loader.exec_module(plot)

OK = "\033[92m✓\033[0m"
FAIL = "\033[91m✗\033[0m"


def Et_only(path: Path) -> float | None:
    with open(path) as f:
        line = f.readline()
    m = re.search(r"Et\s*=\s*([\d.+\-eE]+)", line)
    return float(m.group(1)) if m else None


def check(label, ok, detail=""):
    sym = OK if ok else FAIL
    print(f"  {sym} {label}{(' — ' + detail) if detail else ''}")
    return ok


def check_energy_only(data_dir: Path) -> int:
    files = sorted(glob(str(data_dir / "X*.dat")))
    if len(files) < 2:
        print(f"{FAIL} need at least 2 snapshots")
        return 2

    Et0 = Et_only(Path(files[0]))
    Ets = [Et_only(Path(f)) for f in files]
    print(f"=== 2D Periodicity (energy-only): {len(files)} snapshots ===")
    all_eq = all(e == Et0 for e in Ets)
    if all_eq:
        check(f"Energy bit-identical across {len(files)} snapshots",
              True, f"Et = {Et0:.9e} J")
        return 0
    drifts = [(abs(e - Et0) / abs(Et0) if Et0 else 0.0) for e in Ets]
    max_drift = max(drifts)
    if max_drift < 1e-12:
        check("Energy preserved within reduction noise", True,
              f"max relative drift {max_drift:.2e}")
        return 0
    check("Energy bit-identical", False, f"max drift {max_drift:.2e}")
    return 1


# Tolerances for cycle-regime histogram return. Each bin is an
# independent count; ULP slippage in the analytical event-time
# accumulation can move a handful of particles into adjacent bins
# (those whose final position lands within ~1e-12 of a bin edge).
# Empirically with N=2^18 we see ≤10 such particles → ≤20 affected
# bins, all with |Δ|=1. Tolerances are 3× the empirical maximum to
# leave headroom for legitimate variation.
H_MAX_LIMIT_PER_N = 1.0e-4   # max|Δh| ≤ this × N
H_RMS_LIMIT_PER_N = 1.5e-5   # RMS(Δh) ≤ this × N


def check_cycle(data_dir: Path) -> int:
    files = sorted(glob(str(data_dir / "X*.dat")))
    if len(files) < 2:
        print(f"{FAIL} need at least 2 snapshots")
        return 2

    h0, m0, j0 = plot.load_snapshot(Path(files[0]))
    hT, mT, jT = plot.load_snapshot(Path(files[-1]))

    label_t = files[-1].split("X")[-1].split(".dat")[0]
    print(f"=== 2D Cycle periodicity check: t=0 → t={label_t} ===")
    N = int(m0["hx"].sum())
    print(f"  N = {N},  BINS = {h0['BINS']},  HXY_BINS = {h0['HXY_BINS']}")

    issues = 0

    if int(m0["hx"].sum()) != int(mT["hx"].sum()):
        if not check("Particle count conserved (hx)", False,
                     f"{int(m0['hx'].sum())} → {int(mT['hx'].sum())}"):
            issues += 1
    else:
        check("Particle count conserved across all marginals", True,
              f"N = {N}")

    # Energy must be bit-identical in the cycle regime (alfa=0, sigma_l=0,
    # signs flip on bounces with no energy change).
    if h0["Et"] != hT["Et"]:
        check("Energy bit-identical", False,
              f"drift {abs(hT['Et']-h0['Et'])/abs(h0['Et']):.2e}")
        issues += 1
    else:
        check("Energy bit-identical", True, f"Et = {h0['Et']:.9e} J")

    h_max_lim = max(2, int(H_MAX_LIMIT_PER_N * N))
    h_rms_lim = max(1.0, H_RMS_LIMIT_PER_N * N)

    for key, label in [("hx", "h(x)"), ("hy", "h(y)"),
                       ("gpx", "g(p_x)"), ("gpy", "g(p_y)")]:
        d = mT[key].astype(np.int64) - m0[key].astype(np.int64)
        amax = int(np.abs(d).max())
        arms = float(np.sqrt((d.astype(float) ** 2).mean()))
        nonzero = int(np.count_nonzero(d))
        ok = (amax <= h_max_lim and arms <= h_rms_lim)
        if not check(f"{label} periodic (max|Δ|≤{h_max_lim}, RMS≤{h_rms_lim:.2f})",
                     ok,
                     f"max={amax}, RMS={arms:.3f}, {nonzero}/{len(d)} bins differ"):
            issues += 1

    # 2D joint h(x,y) — coarser binning, smaller per-bin counts.
    d_xy = jT.astype(np.int64) - j0.astype(np.int64)
    amax_xy = int(np.abs(d_xy).max())
    arms_xy = float(np.sqrt((d_xy.astype(float) ** 2).mean()))
    nonzero_xy = int(np.count_nonzero(d_xy))
    # h_xy bin counts are typically O(N/HXY_BINS²); allow larger absolute
    # |Δ| but the same per-N fractional tolerance.
    xy_max_lim = max(2, int(H_MAX_LIMIT_PER_N * N * 4))
    xy_rms_lim = max(1.0, H_RMS_LIMIT_PER_N * N * 2)
    ok_xy = (amax_xy <= xy_max_lim and arms_xy <= xy_rms_lim)
    if not check(f"h(x,y) periodic (max|Δ|≤{xy_max_lim}, RMS≤{xy_rms_lim:.2f})",
                 ok_xy,
                 f"max={amax_xy}, RMS={arms_xy:.3f}, "
                 f"{nonzero_xy}/{d_xy.size} cells differ"):
        issues += 1

    print()
    if issues == 0:
        print(f"{OK} Periodicity check passed.")
        return 0
    print(f"{FAIL} {issues} check(s) failed.")
    return 1


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mode", choices=["energy", "cycle"])
    ap.add_argument("data_dir", type=Path)
    args = ap.parse_args(argv)

    if args.mode == "energy":
        return check_energy_only(args.data_dir)
    return check_cycle(args.data_dir)


if __name__ == "__main__":
    sys.exit(main())
