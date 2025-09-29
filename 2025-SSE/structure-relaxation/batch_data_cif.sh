#!/usr/bin/env bash
# Convert LAMMPS .data files to .cif using lammps_data_to_cif.py

set -uo pipefail  # (no -e; we continue on failures per file)

# ====== SETTINGS ======
BASE_DIR="/home/khw/solid-state-electrolyte/paper/relaxation/uspex/cal/after"
CONVERTER_PY="$BASE_DIR/data_cif.py"  
AFTER_DIR="$BASE_DIR/cif"
PATTERN="${PATTERN:-*.data}"

PYTHON_BIN="${PYTHON_BIN:-python3}"
# ======================

mkdir -p "$AFTER_DIR"

if [[ ! -f "$CONVERTER_PY" ]]; then
  echo "[ERR] Converter not found: $CONVERTER_PY" >&2
  exit 1
fi

found=0; ok=0; skip=0; fail=0

find "$BASE_DIR" -type f -name "$PATTERN" -print0 |
while IFS= read -r -d '' src; do
  ((found++))
  up="$src"
  uspex="cal"  # fallback
  while [[ "$up" != "$BASE_DIR" ]]; do
    b="$(basename "$up")"
    if [[ "$b" == USPEX* ]]; then
      uspex="$b"
      break
    fi
    up="$(dirname "$up")"
  done

  base="$(basename "$src")"
  stem="${base%.data}"
  out="$AFTER_DIR/${uspex}__${stem}.cif"

  if [[ -s "$out" ]]; then
    echo "[SKIP] exists: $(basename "$out")"
    ((skip++))
    continue
  fi

  echo "[RUN ] $src -> $out"
  if "$PYTHON_BIN" "$CONVERTER_PY" "$src" -o "$out"; then
    ((ok++))
  else
    echo "[FAIL] $src" >&2
    ((fail++))
    rm -f "$out" >/dev/null 2>&1 || true
  fi
done

echo "----------------------------------------"
echo "Pattern   : $PATTERN"
echo "Found     : $found"
echo "Converted : $ok"
echo "Skipped   : $skip (already existed)"
echo "Failed    : $fail"
echo "Output dir: $AFTER_DIR"

