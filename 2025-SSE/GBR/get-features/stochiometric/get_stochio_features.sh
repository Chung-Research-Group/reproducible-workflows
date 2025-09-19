#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   bash ./get_stochio_features.sh  [input.txt] [output.csv]
# Example:
#   bash ./get_stochio_features.sh  chemparse_result ./stochiometric_features.csv

if [[ $# -lt 1 || $# -gt 2 ]]; then
  echo "Usage: $0 <input.txt> [output.csv]" 1>&2
  exit 1
fi

infile="$1"
outfile="${2:-stochiometric_features.csv}"

python3 - <<'PY' "$infile" "$outfile"
import sys, re, ast, csv

infile=sys.argv[1]
outfile=sys.argv[2]

with open(infile, 'r', encoding='utf-8') as f:
    text = f.read()

pattern = r"\s*([^:]+):\s*(\{[^}]*\})"
matches = list(re.finditer(pattern, text))

if not matches:
    raise SystemExit("No '<formula>: { ... }' pairs found in input.")

records = []
elements = set()

for m in matches:
    formula = m.group(1).strip()
    comp_str = m.group(2)
    comp = ast.literal_eval(comp_str)
    comp = {str(k): float(v) for k, v in comp.items()}
    records.append((formula, comp))
    elements.update(comp.keys())

elements = sorted(elements)

with open(outfile, 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(["Formula", *elements])
    for formula, comp in records:
        total = sum(comp.values())
        row = [formula]
        for el in elements:
            val = comp.get(el, 0.0)
            frac = (val/total) if total else 0.0
            row.append(f"{frac:.6f}")
        w.writerow(row)
PY

echo "Wrote ${outfile}"

