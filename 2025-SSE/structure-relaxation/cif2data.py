from ovito.io import import_file, export_file
from ovito.modifiers import UnwrapTrajectoriesModifier, PythonScriptModifier
from ovito.data import ParticleType
import numpy as np

with open("search_list.txt", "r") as f:
    structure_list = [line.strip() for line in f if line.strip()]

for struc in structure_list:
    try:
        pipeline = import_file(f"/home/khw/solid-state-electrolyte/paper/relaxation/search/cif/{struc}.cif", input_format="cif")
        export_file(
            pipeline,
            f"/home/khw/solid-state-electrolyte/paper/relaxation/search/data/{struc}.data",
            "lammps/data",restricted_triclinic=True,
            export_type_names=True
        )
        print(f"success: {struc}.cif → {struc}.data")
    except Exception as e:
        print(f"fail: {struc}.cif ({e})")
