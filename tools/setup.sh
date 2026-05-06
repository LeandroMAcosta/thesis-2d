#!/usr/bin/env bash
# One-time setup of the Python venv used by tools/.

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT/tools"

if ! command -v python3 >/dev/null; then
    echo "setup.sh: python3 not found in PATH" >&2
    exit 2
fi

if [[ ! -d .venv ]]; then
    echo "Creating venv at $ROOT/tools/.venv ..."
    python3 -m venv .venv
fi

# shellcheck disable=SC1091
source .venv/bin/activate
python3 -m pip install --quiet --upgrade pip
python3 -m pip install --quiet -r requirements.txt
echo "venv ready at $ROOT/tools/.venv. Activate with: source $ROOT/tools/.venv/bin/activate"
