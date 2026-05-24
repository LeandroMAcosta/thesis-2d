#!/usr/bin/env bash
# 2D periodicity tests. Two regimes:
#   1. periodic-gaussian (alfa=0, sigma_l=0, gaussian momenta)
#      → only checks energy bit-conservation
#   2. cycle (init_uniform_p_magnitude calibrated for 175 periods)
#      → checks every histogram returns to t=0 within ULP tolerance

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BIN="$ROOT/main-event"
# Prefer venv python (has numpy); fall back to plain python3.
VENV_PY="$ROOT/tools/.venv/bin/python"
if [[ -x "$VENV_PY" ]]; then
    PY="$VENV_PY"
else
    PY="${PY:-python3}"
fi

if [[ ! -x "$BIN" ]]; then echo "periodicity: $BIN not built" >&2; exit 2; fi

run_case() {
    local mode="$1" config="$2"
    local TMP=$(mktemp -d)
    trap "rm -rf '$TMP'" RETURN
    cp "$config" "$TMP/config.toml"
    echo "[periodicity:$mode] running sim (config=$(basename "$config"))..."
    (cd "$TMP" && SIM_SEED=42 OMP_NUM_THREADS=8 "$BIN" config.toml > sim.log 2>&1)
    "$PY" "$ROOT/tests/check_periodicity.py" "$mode" "$TMP"
    rm -rf "$TMP"
    trap - RETURN
}

failed=0
run_case energy "$ROOT/tests/config.verify.periodic.toml" || failed=$((failed+1))
run_case cycle  "$ROOT/tests/config.cycle.toml"           || failed=$((failed+1))

if (( failed > 0 )); then
    echo "periodicity: $failed regime(s) failed."
    exit 1
fi
echo "periodicity: all regimes passed."
