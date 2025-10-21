#!/usr/bin/env python3
import os, sys, argparse
import numpy as np
from io_helper import replace_string_in_file, append_line_to_file

# import IO helpers
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../ExtXyz/'))
sys.path.insert(0, parent_dir)
from ase_lammps_io import read_xyz, write_xyz

def main():
    parser = argparse.ArgumentParser(description="ASE deformation via deformation gradient F (LAMMPS-compatible).")
    parser.add_argument("--case", type=str, default=None, help="Tag of the calculation")
    parser.add_argument("--run", action="store_true", help="Conduct QE calculation")

    args = parser.parse_args()

    case_ = args.case

    print("case:", case_)
    script_template = "qe.template"
    temp_dir = "/scratch/asgouros/"
    restart_mode = "restart"
    #restart_mode = "from_scratch"

    for dir_idx in range(1, 7):
        for sign in ('-', '+'):
            idx = f"{dir_idx}{sign}"
            case_idx = "%s%s" %(case_, idx)

            if restart_mode == "restart":
                prefix = case_
            else:
                prefix = case_idx

            xyz_def = f"{case_idx}.xyz"
            atoms_def = read_xyz(xyz_def, frame=0)
            #write_xyz(atoms_def, xyz_def+"_debug")


            # Set the prefix and outdir
            outdir = "%s/%s/" %(temp_dir, prefix)
            os.system('mkdir -p %s' %(outdir))

            script_name = f"{case_idx}"
            script_run = f"{script_name}.in"
            script_out = f"{script_name}.out"

            if not os.path.isfile(script_template):
               print('Cannot find script template: %s' %(script_template))
               sys.exit()
            os.system('cp %s %s' %(script_template, script_run))

            replace_string_in_file(script_run, "F_prefix", prefix)
            replace_string_in_file(script_run, "F_outdir", outdir)
            replace_string_in_file(script_run, "F_restart_mode", restart_mode)

            # Generate cell parameters
            append_line_to_file(script_run, "CELL_PARAMETERS {angstrom}")
            cell = atoms_def.get_cell().array  # rows: a, b, c
            for vec in cell:
                append_line_to_file(script_run, f"{vec[0]:.16f} {vec[1]:.16f} {vec[2]:.16f}")

            # Generate atomic positions
            append_line_to_file(script_run, "ATOMIC_POSITIONS {angstrom}")
            symbols = atoms_def.get_chemical_symbols()
            positions = atoms_def.get_positions()  # Cartesian Å
            for s, (x, y, z) in zip(symbols, positions):
                append_line_to_file(script_run, f"{s} {x:.16f} {y:.16f} {z:.16f}")

            PARA_PREFIX="mpirun -np 6"
            PARA_POSTFIX=" -npool 1 "
            cmd = f"{PARA_PREFIX} pw.x {PARA_POSTFIX} < {script_run} > {script_out} # run command"
            print(cmd)
            if args.run:
                os.system(cmd)

if __name__ == "__main__":
    main()
