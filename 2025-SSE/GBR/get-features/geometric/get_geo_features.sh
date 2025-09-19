#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   bash ./get_geo_features.sh /path/to/cif_dir [output.csv] [zeopp_network_binary]
#
# Example:
#   bash ./get_geo_features.sh ./geometric_features.csv /home/program/zeo++-0.3/network

CIF_DIR="${1:-.}"
OUT_CSV="${2:-geometric_features.csv}"
ZEOPP_BIN="${3:-/home/khw/program/zeo++-0.3/network}"

PROBE_R=0
PROBE_R_HA=0
N_SAMPLES=50000

printf "name,Density (g/cm3),POAV (cm3/g),a,b,c\n" > "$OUT_CSV"

shopt -s nullglob
IFS=$'\n'

for cif in "$CIF_DIR"/*.cif "$CIF_DIR"/*.CIF; do
  [[ -e "$cif" ]] || continue

  name="$(basename "$cif")"
  name="${name%.*}"   # USPEX_0.cif -> USPEX_0
  volpo="${cif%.*}.volpo"

  echo "[*] Processing: $name"

  "$ZEOPP_BIN" -ha -volpo "$PROBE_R" "$PROBE_R_HA" "$N_SAMPLES" "$cif" >/dev/null 2>&1 || true

  density="NA"
  poav="NA"
  if [[ -f "$volpo" ]]; then
    density_val="$(awk '/^@/{for(i=1;i<=NF;i++) if($i=="Density:"){print $(i+1); exit}}' "$volpo")"
    poav_val="$(awk '/^@/{for(i=1;i<=NF;i++) if($i=="POAV_cm^3/g:"){print $(i+1); exit}}' "$volpo")"
    [[ -n "${density_val:-}" ]] && density="$density_val"
    [[ -n "${poav_val:-}"    ]] && poav="$poav_val"
  fi

  a="$(awk 'tolower($1)=="_cell_length_a"{print $2; exit}' "$cif")"
  b="$(awk 'tolower($1)=="_cell_length_b"{print $2; exit}' "$cif")"
  c="$(awk 'tolower($1)=="_cell_length_c"{print $2; exit}' "$cif")"

  clean() {
    local v="${1:-}"
    v="${v%%(*}"
    v="${v//\"/}"; v="${v//\'/}"
    v="${v//,/}"
    printf "%s" "$v" | tr -d '[:space:]'
  }
  a="$(clean "$a")"; b="$(clean "$b")"; c="$(clean "$c")"
  [[ -z "$a" ]] && a="NA"; [[ -z "$b" ]] && b="NA"; [[ -z "$c" ]] && c="NA"

  printf "%s,%s,%s,%s,%s,%s\n" "$name" "$density" "$poav" "$a" "$b" "$c" >> "$OUT_CSV"
done

echo "Done. Wrote: $OUT_CSV"

