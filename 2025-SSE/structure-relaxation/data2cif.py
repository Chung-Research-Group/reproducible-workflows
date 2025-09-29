#!/usr/bin/env python3
"""
Convert a LAMMPS data file (write_data output) to a CIF using pymatgen.

Usage:
    python data2cif.py input.data
    python data2cif.py input.data --out-dir /path/to/after
"""

from __future__ import annotations
import argparse
import re
from pathlib import Path
from typing import Dict, List, Tuple

from pymatgen.core import Lattice, Structure
from pymatgen.io.cif import CifWriter

FLOAT = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"

def _find_bounds(lines: List[str]) -> Tuple[float,float,float,float,float,float,float,float,float,bool]:
    xlo = xhi = ylo = yhi = zlo = zhi = None
    xy = xz = yz = 0.0
    triclinic = False

    rx = re.compile(rf"^\s*({FLOAT})\s+({FLOAT})\s+xlo\s+xhi\b", re.I)
    ry = re.compile(rf"^\s*({FLOAT})\s+({FLOAT})\s+ylo\s+yhi\b", re.I)
    rz = re.compile(rf"^\s*({FLOAT})\s+({FLOAT})\s+zlo\s+zhi\b", re.I)
    rtilt = re.compile(rf"^\s*({FLOAT})\s+({FLOAT})\s+({FLOAT})\s+xy\s+xz\s+yz\b", re.I)

    for ln in lines:
        if xlo is None:
            m = rx.match(ln)
            if m:
                xlo, xhi = float(m.group(1)), float(m.group(2))
                continue
        if ylo is None:
            m = ry.match(ln)
            if m:
                ylo, yhi = float(m.group(1)), float(m.group(2))
                continue
        if zlo is None:
            m = rz.match(ln)
            if m:
                zlo, zhi = float(m.group(1)), float(m.group(2))
                continue
        if not triclinic:
            m = rtilt.match(ln)
            if m:
                xy, xz, yz = (float(m.group(i)) for i in (1, 2, 3))
                triclinic = True

    if None in (xlo, xhi, ylo, yhi, zlo, zhi):
        raise ValueError("Failed to parse x/y/z bounds (xlo xhi, ylo yhi, zlo zhi).")

    return xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz, triclinic

def _parse_atom_type_labels(lines: List[str]) -> Dict[int, str]:
    mapping: Dict[int, str] = {}
    try:
        i = next(idx for idx, ln in enumerate(lines) if ln.strip().lower() == "atom type labels")
    except StopIteration:
        return mapping
    i += 1
    while i < len(lines) and lines[i].strip() == "":
        i += 1
    while i < len(lines):
        ln = lines[i].strip()
        if not ln:
            break
        if re.match(r"^(Masses|Atoms\b|Velocities\b|Bonds\b|Angles\b)", ln, re.I):
            break
        toks = ln.split()
        if len(toks) >= 2 and toks[0].isdigit():
            mapping[int(toks[0])] = toks[1]
        else:
            break
        i += 1
    return mapping

def _parse_atoms_atomic(lines: List[str]) -> Tuple[List[int], List[int], List[List[float]]]:
    start = None
    for idx, ln in enumerate(lines):
        if ln.strip().lower().startswith("atoms # atomic"):
            start = idx; break
        if ln.strip() == "Atoms":
            start = idx; break
    if start is None:
        raise ValueError("Could not find 'Atoms # atomic' (or 'Atoms') section.")

    i = start + 1
    if i < len(lines) and lines[i].strip() == "":
        i += 1

    ids: List[int] = []
    types: List[int] = []
    coords: List[List[float]] = []

    while i < len(lines):
        ln = lines[i].strip()
        if not ln:
            break
        if re.match(r"^(Velocities|Bonds|Angles|Dihedrals|Impropers|Masses|Pair Coeffs|Atom Type Labels)", ln):
            break
        ln = ln.split("#", 1)[0].strip()
        if not ln:
            i += 1; continue

        toks = ln.split()
        if len(toks) < 5:
            if len(toks) >= 6 and toks[0].isdigit() and toks[2].isdigit():
                aid = int(toks[0]); atype = int(toks[2])
                x, y, z = map(float, toks[3:6])
            else:
                raise ValueError(f"Bad atom line: {ln}")
        else:
            aid = int(toks[0]); atype = int(toks[1])
            x, y, z = map(float, toks[2:5])

        ids.append(aid)
        types.append(atype)
        coords.append([x, y, z])
        i += 1

    return ids, types, coords

def _build_lattice(xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz) -> Lattice:
    lx = xhi - xlo
    ly = yhi - ylo
    lz = zhi - zlo
    a = [lx, 0.0, 0.0]
    b = [xy, ly, 0.0]
    c = [xz, yz, lz]
    return Lattice([a, b, c])

def data_to_cif(in_path: Path, out_path: Path | None = None, out_dir: Path | None = None) -> Path:
    lines = in_path.read_text().splitlines()

    xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz, _ = _find_bounds(lines)
    lattice = _build_lattice(xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz)

    t2el = _parse_atom_type_labels(lines)
    if not t2el:
        raise ValueError("Missing 'Atom Type Labels' section with type→element mapping.")

    ids, types, coords = _parse_atoms_atomic(lines)

    shifted = [[x - xlo, y - ylo, z - zlo] for (x, y, z) in coords]
    species = [t2el[int(t)] for t in types]
    struct = Structure(lattice, species, shifted, coords_are_cartesian=True)
    
    if out_dir is not None:
        out_path = Path(out_dir) / (in_path.stem + ".cif")
    elif out_path is None:
        out_path = in_path.with_suffix(".cif")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    CifWriter(struct, symprec=None).write_file(str(out_path))
    return out_path

def main():
    ap = argparse.ArgumentParser(description="Convert LAMMPS 'write_data' file to CIF (pymatgen).")
    ap.add_argument("data", type=Path, help="Input LAMMPS data file")
    ap.add_argument("-o", "--out", type=Path, help="Explicit output CIF path (optional)")
    ap.add_argument("--out-dir", type=Path, help="Directory to save CIF (basename preserved)")
    args = ap.parse_args()

    out = data_to_cif(args.data, args.out, args.out_dir)
    print(f"Wrote: {out}")

if __name__ == "__main__":
    main()

