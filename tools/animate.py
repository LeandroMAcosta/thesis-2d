#!/usr/bin/env python3
"""Generate animated GIFs of the 2D simulation's evolution.

Reads multiple X*.dat snapshots (and optionally graba.dmp.<evolution>
checkpoints) from a directory and produces:

  joint     — h(x,y) heatmap evolving frame by frame
  marginals — 4-panel (hx, hy, gpx, gpy) evolving with theoretical overlays
  scatter   — particle positions evolving (sub-sampled)
  pscatter  — particle (p_x, p_y) evolving (sub-sampled)
  full      — combined dashboard-style animation

Usage:
  animate.py {joint,marginals,scatter,pscatter,full} <snapshots_dir> -o out.gif
"""

from __future__ import annotations

import argparse
import re
import sys
from glob import glob
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import imageio.v2 as imageio

# Reuse parsers from plot.py
sys.path.insert(0, str(Path(__file__).parent))
from plot import (
    load_snapshot,
    load_dump,
    maxwell_boltzmann_1d,
    rayleigh,
    P_SIGMA,
    P_RANGE,
    P_DELTA,
)


def list_snapshots(d: Path):
    files = sorted(glob(str(d / "X*.dat")))
    if not files:
        print(f"animate: no snapshots in {d}", file=sys.stderr); sys.exit(2)
    out = []
    for f in files:
        m = re.search(r"X(\d+)\.dat", Path(f).name)
        evol = int(m.group(1)) if m else 0
        out.append((evol, Path(f)))
    return out


def list_dumps(d: Path):
    files = sorted(glob(str(d / "graba.dmp.*")))
    out = {}
    for f in files:
        m = re.search(r"graba\.dmp\.(\d+)", Path(f).name)
        if m:
            out[int(m.group(1))] = Path(f)
    return out


def render_frame_to_array(fig):
    fig.canvas.draw()
    buf = np.asarray(fig.canvas.buffer_rgba())
    return buf[:, :, :3]  # drop alpha


def cmd_joint(args):
    d = Path(args.dir)
    snaps = list_snapshots(d)
    frames = []
    # Find global vmax for consistent colour scale
    vmax = 0
    for evol, path in snaps:
        _, _, joint = load_snapshot(path)
        vmax = max(vmax, int(joint.max()))
    print(f"animate joint: {len(snaps)} frames, vmax={vmax}")

    for i, (evol, path) in enumerate(snaps):
        header, _, joint = load_snapshot(path)
        fig, ax = plt.subplots(figsize=(8, 7))
        im = ax.imshow(joint.T, origin="lower", extent=[-0.5, 0.5, -0.5, 0.5],
                       cmap="viridis", vmin=0, vmax=vmax, aspect="auto")
        ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]")
        ax.set_title(
            f"h(x, y) — evolution = {evol:>7d}    N = {header['N']}\n"
            rf"$\chi^2_x={header['chi2x']:.3f}$  $\chi^2_y={header['chi2y']:.3f}$  "
            rf"isotropy={header['isotropy']:.4f}  corr={header['corr']:+.2e}"
        )
        plt.colorbar(im, ax=ax, label="counts per cell")
        plt.tight_layout()
        frames.append(render_frame_to_array(fig))
        plt.close(fig)
        print(f"  frame {i+1}/{len(snaps)}")
    imageio.mimsave(args.output, frames, duration=args.duration, loop=0)
    print(f"wrote {args.output} ({len(frames)} frames)")


def cmd_marginals(args):
    d = Path(args.dir)
    snaps = list_snapshots(d)
    # Find max counts for stable y-axes
    h_max = 0; g_max = 0
    for evol, path in snaps:
        header, marg, _ = load_snapshot(path)
        BINS = header["BINS"]
        interior = slice(2, 2 * BINS + 2)
        h_max = max(h_max, marg["hx"][interior].max(), marg["hy"][interior].max())
        g_max = max(g_max, marg["gpx"].max(), marg["gpy"].max())

    frames = []
    for i, (evol, path) in enumerate(snaps):
        header, marg, _ = load_snapshot(path)
        N = header["N"]; BINS = header["BINS"]
        interior = slice(2, 2 * BINS + 2)
        fig, axes = plt.subplots(2, 2, figsize=(13, 9))

        ax = axes[0, 0]
        ax.bar(marg["x"][interior], marg["hx"][interior], width=1.0/BINS,
               color="steelblue", alpha=0.8)
        ax.axhline(N/(2*BINS), color="crimson", linestyle="--")
        ax.set_ylim(0, h_max * 1.15)
        ax.set_xlabel("x"); ax.set_ylabel("counts")
        ax.set_title(rf"$h(x): \chi^2 = {header['chi2x']:.4f}$")
        ax.grid(alpha=0.3)

        ax = axes[0, 1]
        ax.bar(marg["y"][interior], marg["hy"][interior], width=1.0/BINS,
               color="seagreen", alpha=0.8)
        ax.axhline(N/(2*BINS), color="crimson", linestyle="--")
        ax.set_ylim(0, h_max * 1.15)
        ax.set_xlabel("y"); ax.set_ylabel("counts")
        ax.set_title(rf"$h(y): \chi^2 = {header['chi2y']:.4f}$")
        ax.grid(alpha=0.3)

        ax = axes[1, 0]
        ax.bar(marg["p"], marg["gpx"], width=P_DELTA, color="orange", alpha=0.8)
        ax.plot(marg["p"], maxwell_boltzmann_1d(marg["p"], N),
                color="black", linewidth=2)
        ax.set_ylim(0, g_max * 1.15)
        ax.set_xlabel(r"$p_x$"); ax.set_ylabel("counts")
        ax.set_title(rf"$g(p_x): \chi^2 = {header['chi2px']:.4f}$")
        ax.grid(alpha=0.3)

        ax = axes[1, 1]
        ax.bar(marg["p"], marg["gpy"], width=P_DELTA, color="purple", alpha=0.8)
        ax.plot(marg["p"], maxwell_boltzmann_1d(marg["p"], N),
                color="black", linewidth=2)
        ax.set_ylim(0, g_max * 1.15)
        ax.set_xlabel(r"$p_y$"); ax.set_ylabel("counts")
        ax.set_title(rf"$g(p_y): \chi^2 = {header['chi2py']:.4f}$")
        ax.grid(alpha=0.3)

        fig.suptitle(
            f"Evolution = {evol:>7d}    isotropy = {header['isotropy']:.4f}    "
            f"corr = {header['corr']:+.3e}    Et = {header['Et']:.3e}",
            fontsize=11,
        )
        plt.tight_layout()
        frames.append(render_frame_to_array(fig))
        plt.close(fig)
        print(f"  marginals frame {i+1}/{len(snaps)}")
    imageio.mimsave(args.output, frames, duration=args.duration, loop=0)
    print(f"wrote {args.output} ({len(frames)} frames)")


def cmd_scatter(args):
    d = Path(args.dir)
    snaps = list_snapshots(d)
    dumps = list_dumps(d)
    if not dumps:
        print("scatter animation needs graba.dmp.<evolution> files (keep_snapshots=true)",
              file=sys.stderr); sys.exit(2)

    rng = np.random.default_rng(0)
    n_show = args.n_sub
    sample_idx = None  # set on first frame for consistency

    frames = []
    for i, (evol, path) in enumerate(snaps):
        if evol not in dumps:
            continue
        header, _, _ = load_snapshot(path)
        N = header["N"]
        if sample_idx is None:
            n_show = min(n_show, N)
            sample_idx = rng.choice(N, size=n_show, replace=False)
        x, y, _, _, _ = load_dump(dumps[evol], N)
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.scatter(x[sample_idx], y[sample_idx], s=1.5, alpha=0.4, color="steelblue")
        ax.set_xlim(-0.5, 0.5); ax.set_ylim(-0.5, 0.5)
        ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]")
        ax.set_aspect("equal")
        ax.set_title(f"evol = {evol:>7d}    {n_show} de {N} partículas")
        ax.grid(alpha=0.3)
        plt.tight_layout()
        frames.append(render_frame_to_array(fig))
        plt.close(fig)
        print(f"  scatter frame {i+1}/{len(snaps)}")
    imageio.mimsave(args.output, frames, duration=args.duration, loop=0)
    print(f"wrote {args.output}")


def cmd_pscatter(args):
    d = Path(args.dir)
    snaps = list_snapshots(d)
    dumps = list_dumps(d)
    if not dumps:
        print("pscatter needs graba.dmp.<evolution> (keep_snapshots=true)",
              file=sys.stderr); sys.exit(2)

    rng = np.random.default_rng(0)
    n_show = args.n_sub
    sample_idx = None
    frames = []
    for i, (evol, path) in enumerate(snaps):
        if evol not in dumps:
            continue
        header, _, _ = load_snapshot(path)
        N = header["N"]
        if sample_idx is None:
            n_show = min(n_show, N)
            sample_idx = rng.choice(N, size=n_show, replace=False)
        _, _, px, py, _ = load_dump(dumps[evol], N)
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.scatter(px[sample_idx], py[sample_idx], s=1.5, alpha=0.4, color="orange")
        ax.set_xlim(-P_RANGE, P_RANGE); ax.set_ylim(-P_RANGE, P_RANGE)
        ax.set_xlabel(r"$p_x$"); ax.set_ylabel(r"$p_y$")
        ax.set_aspect("equal")
        ax.set_title(f"evol = {evol:>7d}    isotropy = {header['isotropy']:.4f}")
        ax.grid(alpha=0.3)
        plt.tight_layout()
        frames.append(render_frame_to_array(fig))
        plt.close(fig)
        print(f"  pscatter frame {i+1}/{len(snaps)}")
    imageio.mimsave(args.output, frames, duration=args.duration, loop=0)
    print(f"wrote {args.output}")


def main():
    ap = argparse.ArgumentParser(description="2D event-driven animation toolkit")
    ap.add_argument("--duration", type=float, default=0.5,
                    help="seconds per frame (default 0.5)")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("joint")
    p.add_argument("dir")
    p.add_argument("-o", "--output", default="anim_joint.gif")
    p.set_defaults(func=cmd_joint)

    p = sub.add_parser("marginals")
    p.add_argument("dir")
    p.add_argument("-o", "--output", default="anim_marginals.gif")
    p.set_defaults(func=cmd_marginals)

    p = sub.add_parser("scatter")
    p.add_argument("dir")
    p.add_argument("--n-sub", type=int, default=5000)
    p.add_argument("-o", "--output", default="anim_scatter.gif")
    p.set_defaults(func=cmd_scatter)

    p = sub.add_parser("pscatter")
    p.add_argument("dir")
    p.add_argument("--n-sub", type=int, default=5000)
    p.add_argument("-o", "--output", default="anim_pscatter.gif")
    p.set_defaults(func=cmd_pscatter)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
