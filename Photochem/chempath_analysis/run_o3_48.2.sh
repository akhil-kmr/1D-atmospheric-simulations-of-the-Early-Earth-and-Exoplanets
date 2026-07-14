#!/usr/bin/env bash
# O3 reaction budget + Chempath pathway analysis for 48.2 deg SZA.
# Runs 100pc and 1pc PAL cases using saved steady-state atmospheres.
set -euo pipefail

PYTHON="${PYTHON:-/Users/gregcooke/anaconda3/bin/python3}"
DIR="$(cd "$(dirname "$0")" && pwd)"
SZA="${SZA:-48.2}"
TIMEOUT="${TIMEOUT:-3600}"
PALS="${PALS:-100pc 1pc}"

echo "=== O3 reaction budget (SZA ${SZA}) ==="
"$PYTHON" "$DIR/o3_reaction_budget_48.2.py" --pal $PALS

for PAL in $PALS; do
  LAYERS_FILE="$DIR/o3_${SZA}_${PAL}/pathway_layers.csv"
  if [[ ! -f "$LAYERS_FILE" ]]; then
    echo "Missing $LAYERS_FILE" >&2
    exit 1
  fi
  # Skip header; join layer indices for --layers
  LAYERS="$(tail -n +2 "$LAYERS_FILE" | tr '\n' ' ')"

  echo "=== Convert ${PAL} SZA ${SZA} (layers: ${LAYERS}) ==="
  "$PYTHON" "$DIR/steadystate_to_chempath.py" --pal "$PAL" --sza "$SZA" --layers $LAYERS

  echo "=== Pathways ${PAL} SZA ${SZA} ==="
  "$PYTHON" "$DIR/run_pathways.py" --pal "$PAL" --sza "$SZA" --layers $LAYERS --species O3 --timeout "$TIMEOUT"
done

echo "Done. See $DIR/o3_${SZA}_*/ and $DIR/pathways/"
