import os
import numpy as np
from ase.io import read, write
from tqdm import tqdm
import sys

struc_name = sys.argv[1]
output_file_dir = 'npt_output'
traj_dir = 'traj'
opt_output_dir = 'opt_npt'
struc_type = struc_name.split('_')[0]
element_mapping = {
    'PIM-1':{'O':'O', 'H':'C', 'He':'C', 'Li':'C', 'Be':'C', 'B':'C', 'C':'O', 'N':'C', 'F':'C', 'Ne':'N', 'Na':'H','Mg':'C'},
    '6FDA-SBF': {'O':'O', 'H':'C', 'He':'C', 'Li':'C', 'Be':'C', 'B':'N', 'C':'C', 'N':'N', 'F':'H', 'Ne':'F'},
    '6FDA-BSBF': {'O':'O', 'H':'C', 'He':'C', 'Li':'C', 'Be':'C', 'B':'N', 'C':'C', 'N':'N', 'F':'H', 'Ne':'Br', 'Na':'F'},
    '6FBBA-SBF': {'O':'O', 'H':'C', 'He':'C', 'Li':'C', 'Be':'C', 'B':'N', 'C':'C', 'N':'N', 'F':'H', 'Ne':'F', 'Na':'H'},
    '6FBBA-BSBF': {'O':'O', 'H':'C', 'He':'C', 'Li':'C', 'Be':'C', 'B':'N', 'C':'C', 'N':'N', 'F':'H', 'Ne':'Br', 'Na':'F','Mg':'H'}
}



if not os.path.exists(opt_output_dir): os.mkdir(opt_output_dir)

# if not os.path.exists(os.path.join(opt_output_dir, f'{struc_name}_opt.cif')):
if 1:
    output_file = os.path.join(output_file_dir,f'{struc_name}_output.txt')
    trajectory_file = os.path.join(traj_dir,f'{struc_name}.lammpstrj')

    with open(output_file, "r") as file:
        lines = file.readlines()

    vol_change = np.array([float(line.split()[6]) for line in lines[1:]])
    hist, bin_edges = np.histogram(vol_change, bins=2000, density=True)
    max_bin_index = np.argmax(hist)
    mode_value = (bin_edges[max_bin_index] + bin_edges[max_bin_index + 1]) / 2
    frame_number = (np.abs(vol_change - mode_value)).argmin()
    step_number = [float(line.split()[0]) for line in lines[1:]][frame_number]
    atoms = read(trajectory_file, index=frame_number)
    new_symbols = [element_mapping[struc_type].get(symbol, symbol) for symbol in atoms.get_chemical_symbols()]
    atoms.set_chemical_symbols(new_symbols)
    print(f'{struc_name},{frame_number},{step_number}')
    write(os.path.join(opt_output_dir,f'{struc_name}_opt.cif'), atoms)




