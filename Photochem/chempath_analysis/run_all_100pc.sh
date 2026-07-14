#!/usr/bin/env bash
# Chempath pathway analysis for saved Photochem 100pc atmospheres.
# Uses base conda Photochem — no need to re-run simulations.
set -euo pipefail

PYTHON="${PYTHON:-/Users/gregcooke/anaconda3/bin/python3}"
DIR="$(cd "$(dirname "$0")" && pwd)"
LAYERS="${LAYERS:-20 30 40}"
SPECIES="${SPECIES:-O3}"
TIMEOUT="${TIMEOUT:-3600}"

for SZA in 45 48.2 60; do
  echo "=== Converting SZA ${SZA} ==="
  "$PYTHON" "$DIR/steadystate_to_chempath.py" --pal 100pc --sza "$SZA" --layers $LAYERS
  echo "=== Pathways SZA ${SZA} ==="
  "$PYTHON" "$DIR/run_pathways.py" --pal 100pc --sza "$SZA" --layers $LAYERS --species $SPECIES --timeout "$TIMEOUT"
done

echo "Done. See $DIR/pathways/"
