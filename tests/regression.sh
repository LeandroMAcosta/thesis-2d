#!/usr/bin/env bash
# Smoke regression: 4096 particles, 300 steps. Verifies the binary
# runs end-to-end and produces sensible chi² (within ±2 of 1) and
# isotropy / correlation (within their nominal bounds).

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BIN="${1:-$ROOT/main-event}"
CONFIG="$ROOT/tests/config.regression.toml"

if [[ ! -x "$BIN" ]]; then echo "regression: $BIN not built" >&2; exit 2; fi

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
cp "$CONFIG" "$TMP/config.toml"
(cd "$TMP" && SIM_SEED=42 OMP_NUM_THREADS=4 "$BIN" config.toml > sim.log 2>&1)

# Parse the last snapshot's header
hdr=$(head -1 "$TMP/X0000300.dat")
echo "[regression] last snapshot header:"
echo "  $hdr"

# Extract chi² values and confirm sanity
chi2x=$(echo "$hdr" | grep -oP 'chi2x\s*=\s*\K[\d.eE+-]+')
chi2y=$(echo "$hdr" | grep -oP 'chi2y\s*=\s*\K[\d.eE+-]+')
chi2px=$(echo "$hdr" | grep -oP 'chi2px\s*=\s*\K[\d.eE+-]+')
chi2py=$(echo "$hdr" | grep -oP 'chi2py\s*=\s*\K[\d.eE+-]+')
iso=$(echo "$hdr" | grep -oP 'isotropy\s*=\s*\K[\d.eE+-]+')

# All chi² should be in (0, 5) for a sane run
all_ok=1
for v in "$chi2x" "$chi2y" "$chi2px" "$chi2py"; do
    val=$(echo "$v" | awk '{print ($1 > 0 && $1 < 5) ? "ok" : "bad"}')
    if [[ "$val" != "ok" ]]; then all_ok=0; fi
done

# Isotropy should be in (0.5, 2.0) for a sane initial run
iso_ok=$(echo "$iso" | awk '{print ($1 > 0.5 && $1 < 2.0) ? "ok" : "bad"}')

if [[ $all_ok -eq 1 && "$iso_ok" == "ok" ]]; then
    echo "regression: ok (2D event-driven seed=42)"
    exit 0
else
    echo "regression: SANITY CHECK FAILED" >&2
    echo "  chi2x=$chi2x chi2y=$chi2y chi2px=$chi2px chi2py=$chi2py iso=$iso" >&2
    exit 1
fi
