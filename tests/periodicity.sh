#!/usr/bin/env bash
# 2D strict periodicity test (alfa=0, sigma_l=0).

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BIN="$ROOT/main-event"
PY="${PY:-python3}"

if [[ ! -x "$BIN" ]]; then echo "periodicity: $BIN not built" >&2; exit 2; fi

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
cp "$ROOT/tests/config.verify.periodic.toml" "$TMP/config.toml"

echo "[periodicity] running 1M-step periodicity sim..."
(cd "$TMP" && SIM_SEED=42 OMP_NUM_THREADS=8 "$BIN" config.toml > sim.log 2>&1)

"$PY" "$ROOT/tests/check_periodicity.py" "$TMP"
