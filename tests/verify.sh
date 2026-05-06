#!/usr/bin/env bash
# Run all postcondition checks against the relax + periodic configs.

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BIN="$ROOT/main-event"
PY="${PY:-python3}"

if [[ ! -x "$BIN" ]]; then echo "verify: $BIN not built" >&2; exit 2; fi

run_scenario() {
    local name="$1" config="$2" scenario="$3"
    echo "[verify] running $name ($config) ..."
    local TMP
    TMP=$(mktemp -d)
    cp "$ROOT/tests/$config" "$TMP/config.toml"
    (cd "$TMP" && SIM_SEED=42 OMP_NUM_THREADS=8 "$BIN" config.toml > sim.log 2>&1)
    "$PY" "$ROOT/tests/check_postconditions.py" "$scenario" "$TMP"
    rm -rf "$TMP"
}

run_scenario relax    config.verify.relax.toml    relax
run_scenario periodic config.verify.periodic.toml periodic

echo
echo "verify: all scenarios passed."
