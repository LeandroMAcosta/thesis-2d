#!/usr/bin/env bash
# Run the relax config and check p_x, p_y independence specifically.

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BIN="$ROOT/main-event"
PY="${PY:-python3}"

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
cp "$ROOT/tests/config.verify.relax.toml" "$TMP/config.toml"

echo "[independence] running relax sim..."
(cd "$TMP" && SIM_SEED=42 OMP_NUM_THREADS=8 "$BIN" config.toml > sim.log 2>&1)

"$PY" "$ROOT/tests/check_independence.py" "$TMP"
