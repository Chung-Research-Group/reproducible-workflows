from ase.io import read,write
import os
import sys

struc = sys.argv[1]
T = sys.argv[2]
folder_name = f'{struc}_{T}'

save_path = os.path.join('cif_deform_structures',folder_name)
work_path = os.path.join('deform_structures',folder_name)

if not os.path.exists('cif_deform_structures'): os.mkdir('cif_deform_structures')
if not os.path.exists(save_path): os.mkdir(save_path)

strucs_list = [f for f in os.listdir(work_path) if '.data' in f]
for s in strucs_list:
    s_raw = s.split('.data')[0]
    atom=read(os.path.join(work_path,f'{s_raw}.data'),format='lammps-data')
    write(os.path.join(save_path,f'{s_raw}.cif'),atom,format='cif')
    print(f'Saved: {s_raw}.cif')

