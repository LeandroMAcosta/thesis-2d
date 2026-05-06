#!/usr/bin/env python3
"""Plotting toolkit for the 2D event-driven simulation.

Subcommands:
  marginals  — 4-panel: hx, hy, gpx, gpy with theoretical overlays
  joint      — 2-panel heatmap: h(x,y) + g(p_x, p_y) (the latter from
               particles binned at runtime)
  radial     — distribution of |p|, vs Rayleigh reference
  angular    — distribution of arctan2(py, px), vs uniform reference
  energy     — Et vs evolution snapshot index
  scatter    — scatter plot of (x, y) and (p_x, p_y) at a snapshot,
               sub-sampled if N is large
  dashboard  — single PNG with all the above for one snapshot

Reads X*.dat snapshots produced by main-event. For per-particle data
(needed for radial / angular / scatter) reads graba.dmp checkpoints.
"""

from __future__ import annotations

import argparse
import re
import sys
from glob import glob
from pathlib import Path
from typing import Tuple

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# ---------------------------------------------------------------------------
# Constants from constants.h (kept in sync manually).
# ---------------------------------------------------------------------------

P_SIGMA = 5.24684e-24
P_RANGE = 3.0e-23
P_DELTA = 6.0e-26
HXY_BINS_DEFAULT = 64

# ---------------------------------------------------------------------------
# Snapshot parser
# ---------------------------------------------------------------------------


def parse_header(line: str) -> dict:
    out = {}
    for key, pat in [
        ("chi2x", r"chi2x\s*=\s*([\d.+\-eE]+)"),
        ("chi2y", r"chi2y\s*=\s*([\d.+\-eE]+)"),
        ("chi2px", r"chi2px\s*=\s*([\d.+\-eE]+)"),
        ("chi2py", r"chi2py\s*=\s*([\d.+\-eE]+)"),
        ("isotropy", r"isotropy\s*=\s*([\d.+\-eE]+)"),
        ("corr", r"corr\s*=\s*([\d.+\-eE]+)"),
        ("Et", r"Et\s*=\s*([\d.+\-eE]+)"),
        ("N", r"N\s*=\s*(\d+)"),
        ("BINS", r"BINS\s*=\s*(\d+)"),
        ("HXY_BINS", r"HXY_BINS\s*=\s*(\d+)"),
    ]:
        m = re.search(pat, line)
        if m:
            try:
                v = m.group(1)
                out[key] = int(v) if "." not in v and "e" not in v.lower() else float(v)
            except ValueError:
                out[key] = m.group(1)
    return out


def load_snapshot(path: Path):
    """Parse a single .dat snapshot.

    Returns (header_dict, marginals_dict, joint_array)
    where:
      marginals_dict has keys 'x', 'hx', 'y', 'hy', 'p', 'gpx', 'gpy'
        as numpy arrays.
      joint_array is a 2D numpy array (HXY_BINS × HXY_BINS) of int counts.
    """
    text = path.read_text()
    lines = text.splitlines()
    header = parse_header(lines[0])

    BINS = header["BINS"]
    HXY = header["HXY_BINS"]
    n_marg = 2 * BINS + 4

    # Find the start of the marginals block
    i = 0
    while i < len(lines) and not lines[i].startswith("# MARGINALS"):
        i += 1
    i += 2  # skip "# MARGINALS" and "# x_val ..." header

    rows = []
    while i < len(lines) and not lines[i].startswith("# JOINT"):
        if lines[i].startswith("#") or not lines[i].strip():
            i += 1
            continue
        rows.append(lines[i].split())
        i += 1
    arr = np.array(rows, dtype=float)
    marginals = {
        "x": arr[:, 0],
        "hx": arr[:, 1].astype(int),
        "y": arr[:, 2],
        "hy": arr[:, 3].astype(int),
        "p": arr[:, 4],
        "gpx": arr[:, 5].astype(int),
        "gpy": arr[:, 6].astype(int),
    }

    # Joint block
    i += 1  # skip "# JOINT" line
    joint = np.zeros((HXY, HXY), dtype=int)
    for j in range(HXY):
        if i >= len(lines):
            break
        joint[j] = [int(v) for v in lines[i].split()][:HXY]
        i += 1

    return header, marginals, joint


def load_dump(path: Path, N: int) -> Tuple[np.ndarray, ...]:
    """Read a graba.dmp binary checkpoint. Returns (x, y, px, py, evolution)."""
    with open(path, "rb") as f:
        evolution = np.frombuffer(f.read(4), dtype=np.uint32)[0]
        x = np.frombuffer(f.read(8 * N), dtype=np.float64).copy()
        y = np.frombuffer(f.read(8 * N), dtype=np.float64).copy()
        px = np.frombuffer(f.read(8 * N), dtype=np.float64).copy()
        py = np.frombuffer(f.read(8 * N), dtype=np.float64).copy()
    return x, y, px, py, int(evolution)


# ---------------------------------------------------------------------------
# Theoretical overlays
# ---------------------------------------------------------------------------


def maxwell_boltzmann_1d(p_grid, N, sigma=P_SIGMA, p_delta=P_DELTA):
    """Number of counts per bin for a 1D MB distribution scaled to N particles."""
    return (p_delta * N / (sigma * np.sqrt(2 * np.pi))) * np.exp(
        -(p_grid ** 2) / (2 * sigma ** 2)
    )


def rayleigh(p_grid, N, sigma=P_SIGMA):
    """Number of counts per bin (linear bin spacing) for the magnitude
    of a 2D MB momentum: |p| Rayleigh-distributed."""
    # f(p) = p / sigma^2 * exp(-p^2 / 2sigma^2)
    return N * (p_grid / sigma ** 2) * np.exp(-(p_grid ** 2) / (2 * sigma ** 2))


# ---------------------------------------------------------------------------
# Plot subcommands
# ---------------------------------------------------------------------------


def cmd_marginals(args):
    """4-panel: hx, hy, gpx, gpy with theoretical overlays."""
    snap = Path(args.snapshot)
    header, marg, _ = load_snapshot(snap)
    N = header["N"]
    BINS = header["BINS"]

    # Drop the 4 boundary bins of x marginals
    interior = slice(2, 2 * BINS + 2)

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))

    # Top-left: h(x)
    ax = axes[0, 0]
    ax.bar(marg["x"][interior], marg["hx"][interior],
           width=1.0 / BINS, color="steelblue", alpha=0.8, label="h(x)")
    ax.axhline(N / (2 * BINS), color="crimson", linestyle="--",
               label=f"uniform reference = N/(2·BINS)")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("counts")
    ax.set_title(rf"h(x):  $\chi^2 = ${header['chi2x']:.4f}")
    ax.legend(loc="lower center", fontsize=9)
    ax.grid(alpha=0.3)

    # Top-right: h(y)
    ax = axes[0, 1]
    ax.bar(marg["y"][interior], marg["hy"][interior],
           width=1.0 / BINS, color="seagreen", alpha=0.8, label="h(y)")
    ax.axhline(N / (2 * BINS), color="crimson", linestyle="--", label="uniform ref.")
    ax.set_xlabel("y [m]")
    ax.set_ylabel("counts")
    ax.set_title(rf"h(y):  $\chi^2 = ${header['chi2y']:.4f}")
    ax.legend(loc="lower center", fontsize=9)
    ax.grid(alpha=0.3)

    # Bottom-left: g(p_x)
    ax = axes[1, 0]
    ax.bar(marg["p"], marg["gpx"], width=P_DELTA, color="orange", alpha=0.8,
           label=r"$g(p_x)$")
    th = maxwell_boltzmann_1d(marg["p"], N)
    ax.plot(marg["p"], th, color="black", linewidth=2, label="Maxwell-Boltzmann")
    ax.set_xlabel(r"$p_x$ [kg·m/s]")
    ax.set_ylabel("counts")
    ax.set_title(rf"$g(p_x):  \chi^2 = ${header['chi2px']:.4f}")
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(alpha=0.3)

    # Bottom-right: g(p_y)
    ax = axes[1, 1]
    ax.bar(marg["p"], marg["gpy"], width=P_DELTA, color="purple", alpha=0.8,
           label=r"$g(p_y)$")
    ax.plot(marg["p"], th, color="black", linewidth=2, label="MB ref.")
    ax.set_xlabel(r"$p_y$ [kg·m/s]")
    ax.set_ylabel("counts")
    ax.set_title(rf"$g(p_y):  \chi^2 = ${header['chi2py']:.4f}")
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(alpha=0.3)

    fig.suptitle(
        f"{snap.name}   N={N}   isotropy={header['isotropy']:.4f}   "
        f"corr={header['corr']:+.2e}   Et={header['Et']:.3e}",
        fontsize=11,
    )
    plt.tight_layout()
    plt.savefig(args.output, dpi=120, bbox_inches="tight")
    print(f"wrote {args.output}")


def cmd_joint(args):
    """Heatmap of h(x,y). If a graba.dmp is also given, plots g(p_x,p_y) too."""
    snap = Path(args.snapshot)
    header, marg, joint = load_snapshot(snap)

    n_panels = 1 + (1 if args.dump else 0)
    fig, axes = plt.subplots(1, n_panels, figsize=(7 * n_panels, 6),
                             squeeze=False)
    axes = axes.flatten()

    # Position joint
    ax = axes[0]
    extent = [-0.5, 0.5, -0.5, 0.5]
    im = ax.imshow(joint.T, origin="lower", extent=extent,
                   cmap="viridis", aspect="auto")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_title(f"h(x, y) — {snap.name}")
    plt.colorbar(im, ax=ax, label="counts per cell")

    # Momentum joint from particles
    if args.dump:
        N = header["N"]
        x, y, px, py, evol = load_dump(Path(args.dump), N)
        ax = axes[1]
        H, xedges, yedges = np.histogram2d(
            px, py, bins=64, range=[[-P_RANGE, P_RANGE], [-P_RANGE, P_RANGE]]
        )
        im = ax.imshow(H.T, origin="lower",
                       extent=[-P_RANGE, P_RANGE, -P_RANGE, P_RANGE],
                       cmap="magma", aspect="auto")
        ax.set_xlabel(r"$p_x$ [kg·m/s]")
        ax.set_ylabel(r"$p_y$ [kg·m/s]")
        ax.set_title(rf"$g(p_x, p_y)$ — {evol} steps")
        plt.colorbar(im, ax=ax, label="counts per cell")

    plt.tight_layout()
    plt.savefig(args.output, dpi=120, bbox_inches="tight")
    print(f"wrote {args.output}")


def cmd_radial(args):
    """Radial distribution |p| from a checkpoint, vs Rayleigh."""
    if not args.dump:
        print("radial requires --dump <graba.dmp>", file=sys.stderr); sys.exit(2)
    snap_header = None
    if args.snapshot:
        snap_header, _, _ = load_snapshot(Path(args.snapshot))
        N = snap_header["N"]
    else:
        N = args.n if args.n else None
        if N is None:
            print("radial: need --n N or --snapshot to know N", file=sys.stderr)
            sys.exit(2)
    x, y, px, py, evol = load_dump(Path(args.dump), N)
    p_mag = np.sqrt(px ** 2 + py ** 2)

    fig, ax = plt.subplots(figsize=(9, 6))
    counts, edges, _ = ax.hist(p_mag, bins=60, range=(0, 3 * P_SIGMA),
                                color="steelblue", alpha=0.8, label=r"$|\vec p|$")
    centers = 0.5 * (edges[:-1] + edges[1:])
    bin_w = edges[1] - edges[0]
    th = rayleigh(centers, N) * bin_w
    ax.plot(centers, th, color="crimson", linewidth=2,
            label=r"Rayleigh: $f(p) = (p/\sigma^2) \exp(-p^2/2\sigma^2)$")
    ax.set_xlabel(r"$|\vec p|$ [kg·m/s]")
    ax.set_ylabel("counts")
    ax.set_title(rf"Distribución radial de momentos — evol={evol}, N={N}")
    ax.legend(fontsize=10)
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(args.output, dpi=120, bbox_inches="tight")
    print(f"wrote {args.output}")


def cmd_angular(args):
    """Angular distribution θ = arctan2(py, px), should be uniform."""
    if not args.dump:
        print("angular requires --dump <graba.dmp>", file=sys.stderr); sys.exit(2)
    if args.snapshot:
        h, _, _ = load_snapshot(Path(args.snapshot))
        N = h["N"]
    else:
        N = args.n if args.n else None
    x, y, px, py, evol = load_dump(Path(args.dump), N)
    theta = np.arctan2(py, px)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Cartesian histogram
    ax = axes[0]
    ax.hist(theta, bins=60, range=(-np.pi, np.pi), color="purple",
            alpha=0.8, label=r"$\theta$")
    ax.axhline(N / 60, color="crimson", linestyle="--", label="uniform ref.")
    ax.set_xlabel(r"$\theta$ [rad]")
    ax.set_ylabel("counts")
    ax.set_title("Distribución angular")
    ax.legend(fontsize=10)
    ax.grid(alpha=0.3)

    # Polar / circular histogram
    ax = plt.subplot(1, 2, 2, projection="polar")
    counts, edges = np.histogram(theta, bins=60, range=(-np.pi, np.pi))
    ax.bar((edges[:-1] + edges[1:]) / 2, counts,
           width=2 * np.pi / 60, color="purple", alpha=0.8)
    ax.set_title("Polar view", pad=20)

    plt.suptitle(rf"$\theta = \arctan(p_y/p_x)$ — evol={evol}, N={N}", fontsize=11)
    plt.tight_layout()
    plt.savefig(args.output, dpi=120, bbox_inches="tight")
    print(f"wrote {args.output}")


def cmd_scatter(args):
    """Scatter plots of (x, y) and (p_x, p_y) for a checkpoint."""
    if not args.dump:
        print("scatter requires --dump <graba.dmp>", file=sys.stderr); sys.exit(2)
    if args.snapshot:
        h, _, _ = load_snapshot(Path(args.snapshot))
        N = h["N"]
    else:
        N = args.n if args.n else None
    x, y, px, py, evol = load_dump(Path(args.dump), N)

    # Sub-sample for visibility if N large
    n_show = min(args.n_sub, N)
    if n_show < N:
        rng = np.random.default_rng(0)
        idx = rng.choice(N, size=n_show, replace=False)
        x, y, px, py = x[idx], y[idx], px[idx], py[idx]

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    axes[0].scatter(x, y, s=1, alpha=0.3, color="steelblue")
    axes[0].set_xlabel("x [m]"); axes[0].set_ylabel("y [m]")
    axes[0].set_xlim(-0.5, 0.5); axes[0].set_ylim(-0.5, 0.5)
    axes[0].set_title(f"Posiciones ({n_show} de {N} partículas)")
    axes[0].set_aspect("equal")
    axes[0].grid(alpha=0.3)

    axes[1].scatter(px, py, s=1, alpha=0.3, color="orange")
    axes[1].set_xlabel(r"$p_x$ [kg·m/s]"); axes[1].set_ylabel(r"$p_y$ [kg·m/s]")
    axes[1].set_xlim(-P_RANGE, P_RANGE); axes[1].set_ylim(-P_RANGE, P_RANGE)
    axes[1].set_title(rf"Momentos $(p_x, p_y)$")
    axes[1].set_aspect("equal")
    axes[1].grid(alpha=0.3)

    plt.suptitle(f"Scatter — evol={evol}", fontsize=11)
    plt.tight_layout()
    plt.savefig(args.output, dpi=120, bbox_inches="tight")
    print(f"wrote {args.output}")


def cmd_energy(args):
    """Et vs evolution snapshot index across all snapshots in a dir."""
    files = sorted(glob(str(Path(args.dir) / "X*.dat")))
    if not files:
        print(f"no X*.dat in {args.dir}", file=sys.stderr); sys.exit(2)
    Et = []
    evols = []
    for f in files:
        h, _, _ = load_snapshot(Path(f))
        Et.append(h["Et"])
        m = re.search(r"X(\d+)\.dat", Path(f).name)
        evols.append(int(m.group(1)) if m else 0)
    Et = np.array(Et)
    evols = np.array(evols)
    drift = (Et - Et[0]) / Et[0] * 1e6  # ppm

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(evols, drift, marker="o", color="crimson")
    ax.axhline(0, color="black", linewidth=0.5)
    ax.set_xlabel("evolution (steps)")
    ax.set_ylabel(r"$(E_t - E_0)/E_0$ [ppm]")
    ax.set_title(f"Energy drift — {len(files)} snapshots, max |drift| = {abs(drift).max():.2f} ppm")
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(args.output, dpi=120, bbox_inches="tight")
    print(f"wrote {args.output}")


def cmd_dashboard(args):
    """One PNG with marginals + joint + scatter + energy drift."""
    snap = Path(args.snapshot)
    header, marg, joint = load_snapshot(snap)
    N = header["N"]
    BINS = header["BINS"]
    interior = slice(2, 2 * BINS + 2)

    fig = plt.figure(figsize=(20, 12))
    gs = fig.add_gridspec(3, 4, hspace=0.4, wspace=0.35)

    # Row 1: marginals
    ax = fig.add_subplot(gs[0, 0])
    ax.bar(marg["x"][interior], marg["hx"][interior], width=1.0/BINS,
           color="steelblue")
    ax.axhline(N/(2*BINS), color="crimson", linestyle="--")
    ax.set_title(rf"h(x): $\chi^2={header['chi2x']:.3f}$")
    ax.grid(alpha=0.3)

    ax = fig.add_subplot(gs[0, 1])
    ax.bar(marg["y"][interior], marg["hy"][interior], width=1.0/BINS,
           color="seagreen")
    ax.axhline(N/(2*BINS), color="crimson", linestyle="--")
    ax.set_title(rf"h(y): $\chi^2={header['chi2y']:.3f}$")
    ax.grid(alpha=0.3)

    ax = fig.add_subplot(gs[0, 2])
    ax.bar(marg["p"], marg["gpx"], width=P_DELTA, color="orange")
    ax.plot(marg["p"], maxwell_boltzmann_1d(marg["p"], N),
            color="black", linewidth=2)
    ax.set_title(rf"g(p$_x$): $\chi^2={header['chi2px']:.3f}$")
    ax.grid(alpha=0.3)

    ax = fig.add_subplot(gs[0, 3])
    ax.bar(marg["p"], marg["gpy"], width=P_DELTA, color="purple")
    ax.plot(marg["p"], maxwell_boltzmann_1d(marg["p"], N),
            color="black", linewidth=2)
    ax.set_title(rf"g(p$_y$): $\chi^2={header['chi2py']:.3f}$")
    ax.grid(alpha=0.3)

    # Row 2: joints
    ax = fig.add_subplot(gs[1, :2])
    im = ax.imshow(joint.T, origin="lower", extent=[-0.5, 0.5, -0.5, 0.5],
                   cmap="viridis", aspect="auto")
    ax.set_title("h(x, y)")
    ax.set_xlabel("x"); ax.set_ylabel("y")
    plt.colorbar(im, ax=ax)

    if args.dump:
        x, y, px, py, evol = load_dump(Path(args.dump), N)
        ax = fig.add_subplot(gs[1, 2:])
        H, _, _ = np.histogram2d(px, py, bins=64,
                                 range=[[-P_RANGE, P_RANGE], [-P_RANGE, P_RANGE]])
        im = ax.imshow(H.T, origin="lower",
                       extent=[-P_RANGE, P_RANGE, -P_RANGE, P_RANGE],
                       cmap="magma", aspect="auto")
        ax.set_title(r"g(p$_x$, p$_y$) — from particles")
        ax.set_xlabel(r"$p_x$"); ax.set_ylabel(r"$p_y$")
        plt.colorbar(im, ax=ax)

        # Row 3: radial + angular + energy drift
        ax = fig.add_subplot(gs[2, 0])
        p_mag = np.sqrt(px**2 + py**2)
        c, ed, _ = ax.hist(p_mag, bins=60, range=(0, 3*P_SIGMA),
                            color="steelblue", alpha=0.8)
        ce = 0.5*(ed[:-1] + ed[1:])
        bw = ed[1] - ed[0]
        ax.plot(ce, rayleigh(ce, N) * bw, color="crimson", linewidth=2)
        ax.set_title(r"|$\vec p$| (Rayleigh)")
        ax.grid(alpha=0.3)

        ax = fig.add_subplot(gs[2, 1], projection="polar")
        theta = np.arctan2(py, px)
        cnt, ed = np.histogram(theta, bins=60, range=(-np.pi, np.pi))
        ax.bar((ed[:-1]+ed[1:])/2, cnt, width=2*np.pi/60,
               color="purple", alpha=0.8)
        ax.set_title(r"$\theta$ uniform", pad=15)

    # Energy drift across all snapshots in same dir
    snap_dir = snap.parent
    files = sorted(glob(str(snap_dir / "X*.dat")))
    if len(files) > 1:
        Ets = []; evs = []
        for f in files:
            h, _, _ = load_snapshot(Path(f))
            Ets.append(h["Et"])
            m = re.search(r"X(\d+)\.dat", Path(f).name)
            evs.append(int(m.group(1)) if m else 0)
        Ets = np.array(Ets); evs = np.array(evs)
        drift = (Ets - Ets[0]) / Ets[0] * 1e6
        ax = fig.add_subplot(gs[2, 2:])
        ax.plot(evs, drift, marker="o", color="crimson")
        ax.axhline(0, color="black", linewidth=0.5)
        ax.set_title(f"Energy drift (max |drift| = {abs(drift).max():.2f} ppm)")
        ax.set_xlabel("evolution"); ax.set_ylabel(r"$\Delta E/E$ [ppm]")
        ax.grid(alpha=0.3)

    fig.suptitle(
        f"{snap.name}   N={N}   isotropy={header['isotropy']:.4f}   "
        f"corr={header['corr']:+.2e}   Et={header['Et']:.4e}",
        fontsize=12,
    )
    plt.savefig(args.output, dpi=110, bbox_inches="tight")
    print(f"wrote {args.output}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main():
    ap = argparse.ArgumentParser(description="2D event-driven plotting toolkit")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("marginals")
    p.add_argument("snapshot")
    p.add_argument("-o", "--output", default="marginals.png")
    p.set_defaults(func=cmd_marginals)

    p = sub.add_parser("joint")
    p.add_argument("snapshot")
    p.add_argument("--dump", help="optional graba.dmp for the momentum joint")
    p.add_argument("-o", "--output", default="joint.png")
    p.set_defaults(func=cmd_joint)

    p = sub.add_parser("radial")
    p.add_argument("--snapshot", help="optional, for N inference")
    p.add_argument("--dump", required=True)
    p.add_argument("--n", type=int, help="N if no snapshot given")
    p.add_argument("-o", "--output", default="radial.png")
    p.set_defaults(func=cmd_radial)

    p = sub.add_parser("angular")
    p.add_argument("--snapshot")
    p.add_argument("--dump", required=True)
    p.add_argument("--n", type=int)
    p.add_argument("-o", "--output", default="angular.png")
    p.set_defaults(func=cmd_angular)

    p = sub.add_parser("scatter")
    p.add_argument("--snapshot")
    p.add_argument("--dump", required=True)
    p.add_argument("--n", type=int)
    p.add_argument("--n-sub", type=int, default=20000,
                   help="max particles to plot (default 20k)")
    p.add_argument("-o", "--output", default="scatter.png")
    p.set_defaults(func=cmd_scatter)

    p = sub.add_parser("energy")
    p.add_argument("dir")
    p.add_argument("-o", "--output", default="energy.png")
    p.set_defaults(func=cmd_energy)

    p = sub.add_parser("dashboard")
    p.add_argument("snapshot")
    p.add_argument("--dump", help="optional graba.dmp for momentum / radial / angular")
    p.add_argument("-o", "--output", default="dashboard.png")
    p.set_defaults(func=cmd_dashboard)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
