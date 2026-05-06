#!/usr/bin/env python3
"""Live TUI monitor for the 2D event-driven simulation.

Spawns ./main-event as a subprocess, parses stdout for energy /
batch progress, runs lightweight reads of the freshest snapshot when
keep_snapshots is on, and renders a multi-panel terminal view with
plotext.

Usage:
  ./tools/tui [config.toml]

Default: tests/config.verify.relax.toml from the repo root.

When the simulation finishes, the final frame stays on screen until
you press q / Q / Esc / Ctrl-C.
"""

from __future__ import annotations

import argparse
import os
import re
import select
import shutil
import signal
import subprocess
import sys
import termios
import time
import tty
from glob import glob
from pathlib import Path

import plotext as plt

ROOT = Path(__file__).resolve().parent.parent

# ---------------------------------------------------------------------------
# Snapshot reader (slim subset of plot.py)
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
            v = m.group(1)
            try:
                out[key] = int(v) if "." not in v and "e" not in v.lower() else float(v)
            except ValueError:
                out[key] = v
    return out


def latest_snapshot(work_dir: Path):
    files = sorted(glob(str(work_dir / "X*.dat")))
    return Path(files[-1]) if files else None


def read_marginals(snap: Path):
    if snap is None or not snap.exists():
        return None, None, None
    text = snap.read_text()
    lines = text.splitlines()
    if not lines:
        return None, None, None
    header = parse_header(lines[0])
    BINS = header.get("BINS", 500)
    HXY = header.get("HXY_BINS", 64)
    # Find marginals
    i = 0
    while i < len(lines) and not lines[i].startswith("# MARGINALS"):
        i += 1
    i += 2
    rows = []
    while i < len(lines) and not lines[i].startswith("# JOINT"):
        if lines[i].startswith("#") or not lines[i].strip():
            i += 1; continue
        parts = lines[i].split()
        if len(parts) >= 7:
            rows.append((float(parts[0]), int(parts[1]), int(parts[3]),
                         float(parts[4]), int(parts[5]), int(parts[6])))
        i += 1
    if not rows:
        return header, None, None
    # Joint
    i += 1
    joint = []
    for _ in range(HXY):
        if i >= len(lines): break
        joint.append([int(v) for v in lines[i].split()][:HXY])
        i += 1
    return header, rows, joint


# ---------------------------------------------------------------------------
# Render
# ---------------------------------------------------------------------------


def color_for_chi2(v):
    if v is None: return "white"
    if 0.5 < v < 2.0: return "green"
    if 0.2 < v < 5.0: return "yellow"
    return "red"


def color_for_iso(v):
    if v is None: return "white"
    d = abs(1 - v)
    if d < 0.05: return "green"
    if d < 0.15: return "yellow"
    return "red"


def color_for_corr(v):
    if v is None: return "white"
    a = abs(v)
    if a < 0.02: return "green"
    if a < 0.1: return "yellow"
    return "red"


def render(progress_lines, snap_header, snap_marginals, snap_joint, status):
    plt.clf()
    plt.theme("clear")

    rows, cols = shutil.get_terminal_size((100, 30))

    # Make a simple 2-column layout: left for marginals, right for joint
    # plotext's matrix() suits the joint heatmap.
    half_cols = max(40, cols // 2 - 2)

    # Top: status line
    print("\033[2J\033[H", end="")  # clear screen
    print(f"\033[1;36m=== thesis-2d live monitor ===\033[0m   {status}\n")

    if snap_header:
        def cv(c, s):
            cmap = {"red": "31", "green": "32", "yellow": "33", "white": "37"}
            return f"\033[1;{cmap.get(c,'37')}m{s}\033[0m"

        chi2x = snap_header.get("chi2x", 0.0)
        chi2y = snap_header.get("chi2y", 0.0)
        chi2px = snap_header.get("chi2px", 0.0)
        chi2py = snap_header.get("chi2py", 0.0)
        iso = snap_header.get("isotropy", 0.0)
        corr = snap_header.get("corr", 0.0)
        Et = snap_header.get("Et", 0.0)
        N = snap_header.get("N", "?")

        s_x = cv(color_for_chi2(chi2x), f"{chi2x:.3f}")
        s_y = cv(color_for_chi2(chi2y), f"{chi2y:.3f}")
        s_px = cv(color_for_chi2(chi2px), f"{chi2px:.3f}")
        s_py = cv(color_for_chi2(chi2py), f"{chi2py:.3f}")
        s_iso = cv(color_for_iso(iso), f"{iso:.4f}")
        s_corr = cv(color_for_corr(corr), f"{corr:+.3e}")

        print(f"  N={N:>9}     chi2x={s_x}  chi2y={s_y}  chi2px={s_px}  chi2py={s_py}")
        print(f"  Et={Et:.4e}   isotropy={s_iso}  corr={s_corr}")
    print()

    # Plot marginals
    if snap_marginals:
        xs = [r[0] for r in snap_marginals]
        hxs = [r[1] for r in snap_marginals]
        hys = [r[2] for r in snap_marginals]
        ps = [r[3] for r in snap_marginals]
        gpxs = [r[4] for r in snap_marginals]
        gpys = [r[5] for r in snap_marginals]

        plt.clf()
        plt.subplots(2, 2)
        plt.subplot(1, 1).plot_size(half_cols, max(8, rows//4))
        plt.bar(xs, hxs, label="h(x)", color="blue")
        plt.title("h(x)")

        plt.subplot(1, 2).plot_size(half_cols, max(8, rows//4))
        plt.bar(xs, hys, label="h(y)", color="green")
        plt.title("h(y)")

        plt.subplot(2, 1).plot_size(half_cols, max(8, rows//4))
        plt.bar(ps, gpxs, label="g(px)", color="orange")
        plt.title("g(p_x)")

        plt.subplot(2, 2).plot_size(half_cols, max(8, rows//4))
        plt.bar(ps, gpys, label="g(py)", color="magenta")
        plt.title("g(p_y)")

        plt.show()

    print()
    # Last few progress lines from main-event stdout
    print("\033[2;37m--- progress (last 8 lines from main-event) ---\033[0m")
    for line in progress_lines[-8:]:
        print(f"  \033[37m{line}\033[0m")


# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("config", nargs="?",
                    default=str(ROOT / "tests" / "config.verify.relax.toml"))
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--refresh", type=float, default=0.5,
                    help="seconds between renders (default 0.5)")
    args = ap.parse_args()

    bin_path = ROOT / "main-event"
    if not bin_path.exists():
        print("tui: main-event not built. Run `make`.", file=sys.stderr); sys.exit(2)

    # Set up working dir with config + dump enabled
    import tempfile
    work_dir = Path(tempfile.mkdtemp(prefix="tui_2d_"))
    shutil.copy(args.config, work_dir / "config.toml")

    # Force dump=true so we can read snapshots
    cfg_path = work_dir / "config.toml"
    text = cfg_path.read_text()
    if "dump" in text:
        text = re.sub(r"^dump\s*=\s*\w+", "dump = true", text, flags=re.MULTILINE)
    text += "\nkeep_snapshots = true\n"
    cfg_path.write_text(text)

    env = os.environ.copy()
    env["SIM_SEED"] = str(args.seed)
    env["OMP_NUM_THREADS"] = str(args.threads)

    proc = subprocess.Popen(
        [str(bin_path), "config.toml"],
        cwd=work_dir,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        env=env,
        text=True,
        bufsize=1,
    )

    progress_lines = []
    print("\033[?25l", end="")  # hide cursor
    last_render = 0.0
    try:
        while True:
            # Read available stdout (non-blocking)
            ready = select.select([proc.stdout], [], [], 0.05)[0]
            if ready:
                line = proc.stdout.readline()
                if line:
                    progress_lines.append(line.rstrip())
                else:
                    # EOF
                    if proc.poll() is not None:
                        break
            now = time.monotonic()
            if now - last_render > args.refresh:
                snap = latest_snapshot(work_dir)
                header, marginals, joint = read_marginals(snap) if snap else (None, None, None)
                status = "RUNNING" if proc.poll() is None else "DONE"
                render(progress_lines, header, marginals, joint, status)
                last_render = now

            if proc.poll() is not None and not ready:
                break

        # Final render
        snap = latest_snapshot(work_dir)
        header, marginals, joint = read_marginals(snap) if snap else (None, None, None)
        render(progress_lines, header, marginals, joint, "DONE — press q to exit")
        # Wait for keypress
        old = termios.tcgetattr(sys.stdin)
        try:
            tty.setraw(sys.stdin.fileno())
            while True:
                ch = sys.stdin.read(1)
                if ch in ("q", "Q", "\x1b", "\x03"):
                    break
        finally:
            termios.tcsetattr(sys.stdin, termios.TCSADRAIN, old)

    finally:
        print("\033[?25h", end="")  # show cursor
        if proc.poll() is None:
            proc.send_signal(signal.SIGTERM)
            proc.wait(timeout=5)
        # Don't auto-rm work_dir so user can inspect


if __name__ == "__main__":
    main()
