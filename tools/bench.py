#!/usr/bin/env python3
"""Lightweight benchmark recorder for the 2D event-driven simulation.

Stores runs in a local SQLite DB (bench.db at the repo root). Each run
captures: timestamp, config hash, N, total steps, threads, seed,
wall_seconds, steps_per_sec, final chi² (4 of them), isotropy,
correlation, energy drift, host, git commit, notes.

Subcommands:
  run     CONFIG [--threads N] [--seed N] [--notes "..."]
  list    [--limit N]
  show    ID
  compare           — group by config_hash, compare across runs
  export  csv|md
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import os
import re
import socket
import sqlite3
import subprocess
import sys
import tempfile
import time
from glob import glob
from pathlib import Path
from shutil import copy

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_DB = ROOT / "bench.db"

SCHEMA = """
CREATE TABLE IF NOT EXISTS runs (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    timestamp TEXT NOT NULL,
    binary TEXT NOT NULL,
    config_path TEXT NOT NULL,
    config_hash TEXT NOT NULL,
    n_part INTEGER NOT NULL,
    total_steps INTEGER NOT NULL,
    n_threads INTEGER,
    seed INTEGER,
    wall_seconds REAL NOT NULL,
    steps_per_sec REAL,
    final_chi2x REAL, final_chi2y REAL, final_chi2px REAL, final_chi2py REAL,
    final_isotropy REAL, final_corr REAL,
    final_Et REAL,
    energy_drift_ppm REAL,
    git_commit TEXT,
    host TEXT,
    cpu_model TEXT,
    notes TEXT
);
CREATE INDEX IF NOT EXISTS idx_runs_config_hash ON runs(config_hash);
"""


def detect_git_commit() -> str:
    try:
        r = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                           cwd=ROOT, capture_output=True, text=True, check=True)
        return r.stdout.strip()
    except Exception:
        return "unknown"


def detect_cpu() -> str:
    try:
        with open("/proc/cpuinfo") as f:
            for line in f:
                if line.startswith("model name"):
                    return line.split(":", 1)[1].strip()
    except Exception:
        pass
    return "unknown"


def hash_config(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()[:16]


def parse_config(path: Path):
    """Extract n_part and sum of steps from a TOML config."""
    text = path.read_text()
    n_part = int(re.search(r"n_part\s*=\s*(\d+)", text).group(1))
    steps_raw = re.search(r"steps\s*=\s*\[(.*?)\]", text, re.DOTALL).group(1)
    steps = [int(s.strip()) for s in re.findall(r"\d+", steps_raw)]
    total = sum(steps)
    return n_part, total


def parse_final_metrics(work_dir: Path):
    """From the latest X*.dat in work_dir, extract chi², isotropy, corr, Et."""
    files = sorted(glob(str(work_dir / "X*.dat")))
    if not files:
        return None
    with open(files[-1]) as f:
        line = f.readline()
    out = {}
    for k, p in [("chi2x", r"chi2x\s*=\s*([\d.+\-eE]+)"),
                 ("chi2y", r"chi2y\s*=\s*([\d.+\-eE]+)"),
                 ("chi2px", r"chi2px\s*=\s*([\d.+\-eE]+)"),
                 ("chi2py", r"chi2py\s*=\s*([\d.+\-eE]+)"),
                 ("isotropy", r"isotropy\s*=\s*([\d.+\-eE]+)"),
                 ("corr", r"corr\s*=\s*([\d.+\-eE]+)"),
                 ("Et", r"Et\s*=\s*([\d.+\-eE]+)")]:
        m = re.search(p, line)
        if m:
            out[k] = float(m.group(1))
    # also Et at first snapshot for drift
    with open(files[0]) as f:
        line0 = f.readline()
    m0 = re.search(r"Et\s*=\s*([\d.+\-eE]+)", line0)
    if m0:
        Et0 = float(m0.group(1))
        if out.get("Et"):
            out["energy_drift_ppm"] = (out["Et"] - Et0) / Et0 * 1e6
    return out


def open_db(path):
    db = sqlite3.connect(path)
    db.executescript(SCHEMA)
    db.row_factory = sqlite3.Row
    return db


# ---------------------------------------------------------------------------


def cmd_run(args):
    binary = ROOT / "main-event"
    if not binary.exists():
        print(f"bench: {binary} not built. Run `make` first.", file=sys.stderr)
        sys.exit(2)

    cfg = Path(args.config).resolve()
    n_part, total_steps = parse_config(cfg)

    work_dir = Path(tempfile.mkdtemp(prefix="bench_2d_"))
    copy(cfg, work_dir / "config.toml")

    env = os.environ.copy()
    env["SIM_SEED"] = str(args.seed)
    env["OMP_NUM_THREADS"] = str(args.threads)

    print(f"bench: running {binary.name}  config={cfg.name}  N={n_part}  steps={total_steps}  "
          f"threads={args.threads}  seed={args.seed}")
    t0 = time.monotonic()
    proc = subprocess.run([str(binary), "config.toml"],
                          cwd=work_dir, env=env,
                          capture_output=True, text=True, timeout=600)
    wall = time.monotonic() - t0
    if proc.returncode != 0:
        print(proc.stdout, file=sys.stderr)
        print(proc.stderr, file=sys.stderr)
        print(f"bench: binary exited {proc.returncode}", file=sys.stderr)
        sys.exit(1)

    metrics = parse_final_metrics(work_dir) or {}
    record = {
        "timestamp": dt.datetime.utcnow().isoformat() + "Z",
        "binary": str(binary),
        "config_path": str(cfg),
        "config_hash": hash_config(cfg),
        "n_part": n_part,
        "total_steps": total_steps,
        "n_threads": args.threads,
        "seed": args.seed,
        "wall_seconds": wall,
        "steps_per_sec": total_steps / wall if wall > 0 else None,
        "final_chi2x": metrics.get("chi2x"),
        "final_chi2y": metrics.get("chi2y"),
        "final_chi2px": metrics.get("chi2px"),
        "final_chi2py": metrics.get("chi2py"),
        "final_isotropy": metrics.get("isotropy"),
        "final_corr": metrics.get("corr"),
        "final_Et": metrics.get("Et"),
        "energy_drift_ppm": metrics.get("energy_drift_ppm"),
        "git_commit": detect_git_commit(),
        "host": socket.gethostname(),
        "cpu_model": detect_cpu(),
        "notes": args.notes,
    }
    db = open_db(args.db)
    cols = ", ".join(record.keys())
    placeholders = ", ".join("?" * len(record))
    cur = db.execute(f"INSERT INTO runs ({cols}) VALUES ({placeholders})",
                     list(record.values()))
    db.commit()
    rid = cur.lastrowid
    print(f"  ✓ recorded as run #{rid}")
    print(f"    wall          : {wall:.2f} s")
    print(f"    steps/sec     : {record['steps_per_sec']/1e6:.2f} M")
    if metrics.get("energy_drift_ppm") is not None:
        print(f"    energy drift  : {metrics['energy_drift_ppm']:+.2f} ppm")
    print(f"    chi²x/y/px/py : {metrics.get('chi2x',0):.3f} / "
          f"{metrics.get('chi2y',0):.3f} / "
          f"{metrics.get('chi2px',0):.3f} / "
          f"{metrics.get('chi2py',0):.3f}")
    print(f"    isotropy      : {metrics.get('isotropy',0):.4f}")
    print(f"    correlation   : {metrics.get('corr',0):+.4e}")
    # Cleanup work_dir
    import shutil as sh; sh.rmtree(work_dir)


def cmd_list(args):
    db = open_db(args.db)
    rows = db.execute(
        f"SELECT id, timestamp, n_part, total_steps, n_threads, "
        f"wall_seconds, steps_per_sec, notes FROM runs "
        f"ORDER BY id DESC LIMIT ?", (args.limit,)).fetchall()
    print(f" id  ts                 N         steps      thr  wall(s)   Mst/s  notes")
    print("-" * 100)
    for r in rows:
        ts = r["timestamp"][:19]
        print(f" {r['id']:>3}  {ts}  {r['n_part']:>9,}  {r['total_steps']:>9,}  "
              f"{r['n_threads'] or '?':>3}  {r['wall_seconds']:>7.2f}  "
              f"{(r['steps_per_sec'] or 0)/1e6:>5.2f}  {(r['notes'] or '')[:40]}")


def cmd_show(args):
    db = open_db(args.db)
    r = db.execute("SELECT * FROM runs WHERE id = ?", (args.id,)).fetchone()
    if not r:
        print(f"run #{args.id} not found", file=sys.stderr); sys.exit(1)
    for k in r.keys():
        print(f"  {k:<22} {r[k]}")


def cmd_compare(args):
    db = open_db(args.db)
    rows = db.execute(
        "SELECT config_hash, MIN(wall_seconds) AS best, COUNT(*) AS n, "
        "n_part, total_steps FROM runs GROUP BY config_hash ORDER BY n_part, total_steps").fetchall()
    if not rows:
        print("no runs"); return
    print(f"config_hash       N         steps     best wall(s)  Mst/s  runs")
    for r in rows:
        mst = r["total_steps"] / r["best"] / 1e6 if r["best"] > 0 else 0
        print(f"{r['config_hash']:<16}  {r['n_part']:>9,}  {r['total_steps']:>9,}  "
              f"{r['best']:>10.3f}  {mst:>5.2f}  {r['n']}")


def cmd_export(args):
    db = open_db(args.db)
    rows = db.execute("SELECT * FROM runs ORDER BY id").fetchall()
    cols = list(rows[0].keys()) if rows else []
    if args.fmt == "csv":
        import csv
        w = csv.writer(sys.stdout)
        w.writerow(cols)
        for r in rows:
            w.writerow([r[k] for k in cols])
    elif args.fmt == "md":
        print("| " + " | ".join(cols) + " |")
        print("|" + "|".join("---" for _ in cols) + "|")
        for r in rows:
            print("| " + " | ".join(str(r[k]) for k in cols) + " |")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", default=str(DEFAULT_DB))
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("run")
    p.add_argument("config")
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--notes", default="")
    p.set_defaults(func=cmd_run)

    p = sub.add_parser("list")
    p.add_argument("--limit", type=int, default=20)
    p.set_defaults(func=cmd_list)

    p = sub.add_parser("show")
    p.add_argument("id", type=int)
    p.set_defaults(func=cmd_show)

    p = sub.add_parser("compare")
    p.set_defaults(func=cmd_compare)

    p = sub.add_parser("export")
    p.add_argument("fmt", choices=["csv", "md"])
    p.set_defaults(func=cmd_export)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
