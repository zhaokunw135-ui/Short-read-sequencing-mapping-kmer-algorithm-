#!/bin/bash
# Unit tests

PROJECT_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$PROJECT_ROOT" || exit 1

echo "Running Unit Tests"
echo "---------------------------------------------------"
PYTHON_BIN="python3"
if [ -x "$PROJECT_ROOT/.venv/bin/python" ]; then
    PYTHON_BIN="$PROJECT_ROOT/.venv/bin/python"
fi
PYTHONPATH="$PROJECT_ROOT/scripts${PYTHONPATH:+:$PYTHONPATH}" \
    "$PYTHON_BIN" -m unittest discover -s tests -v
status=$?

echo "---------------------------------------------------"
echo "Unit tests complete!"
exit "$status"


