#!/usr/bin/env bash
# CI script: runs all tests in sequence.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

echo "===== build ====="
make >/dev/null

echo "===== regression ====="
bash tests/regression.sh

echo
echo "===== verify (relax + periodic) ====="
bash tests/verify.sh

echo
echo "===== periodicity ====="
bash tests/periodicity.sh

echo
echo "===== isotropy ====="
bash tests/isotropy.sh

echo
echo "===== independence ====="
bash tests/independence.sh

echo
echo "===== ALL TESTS PASSED ====="
