#!/usr/bin/env python3
import json
import os, sys, argparse
import numpy as np
from io_helper import replace_string_in_file, append_line_to_file

# Add the parent directory to sys.path
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../ParseQeOutput/'))
sys.path.insert(0, parent_dir)
from parse_qe_output import parse_qe_output

# import IO helpers
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../ExtXyz/'))
sys.path.insert(0, parent_dir)
from ase_lammps_io import read_xyz, write_xyz

kAtmToGpa = 0.000101325

def main():
    parser = argparse.ArgumentParser(description="ASE deformation via deformation gradient F (LAMMPS-compatible).")
    parser.add_argument("--case", type=str, default=None, help="Tag of the calculation")
    args = parser.parse_args()

    case_ = args.case

    print("case:", case_)
    script_run = case_

    for dir_idx in range(1, 7):
        for sign in ('-', '+'):
            idx = f"{dir_idx}{sign}"
            case_idx = "%s%s" %(case_, idx)
            script_out = f"{case_idx}.out"

            print()
            print('####################################')
            print('Parsing file {:s}'.format(script_out))
            frame_def, status = parse_qe_output(script_out, xyzfmt = 'ase', debug = False)

            for key in frame_def['cell_properties']:
                print(key, frame_def['cell_properties'][key])

            stress_tensor_atm = frame_def['cell_properties']['stress']

            stress_tensor_GPa_dict = {}
            for i in range(9):
                idx = str(i + 1)
                stress_tensor_GPa_dict[idx] = stress_tensor_atm[i]*kAtmToGpa
            
            fname = f"{case_}{dir_idx}{sign}.json"
            with open(fname, 'w') as f: 
                json.dump(stress_tensor_GPa_dict, f)

if __name__ == "__main__":
    main()
