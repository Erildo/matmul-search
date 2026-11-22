#!/usr/bin/env bash
# Run the whole test suite with a hard outer limit. Never hangs:
#  - each test is capped at 10s (tests/conftest.py)
#  - the whole run is killed after ${MAX_SECONDS:-120}s
set -euo pipefail
cd "$(dirname "$0")"
exec timeout --signal=KILL "${MAX_SECONDS:-120}" python3 -m pytest "$@"
