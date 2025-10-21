#!/usr/bin/env python3
import os, sys, argparse
import numpy as np

# import IO helpers
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../ExtXyz/'))
sys.path.insert(0, parent_dir)
from ase_lammps_io import read_xyz, write_xyz

def deformation_gradient(dir_idx: int, s: float) -> np.ndarray:
    """
    Build F for Voigt dir: 1=xx, 2=yy, 3=zz, 4=yz, 5=xz, 6=xy.
    s = ±up (engineering strain / shear).
    """
    F = np.eye(3)
    if   dir_idx == 1:  # xx
        F[0, 0] += s
    elif dir_idx == 2:  # yy
        F[1, 1] += s
    elif dir_idx == 3:  # zz
        F[2, 2] += s
    elif dir_idx == 4:  # yz: y' += s * z
        F[1, 2] += s
    elif dir_idx == 5:  # xz: x' += s * z
        F[0, 2] += s
    elif dir_idx == 6:  # xy: x' += s * y
        F[0, 1] += s
    else:
        raise ValueError("dir_idx must be 1..6")
    return F

def apply_deformation_F(atoms, F: np.ndarray):
    """
    Affine map using F with ASE's row-vector cell convention:
    C_new = C_old @ F.T, positions remapped consistently.
    """
    deformed = atoms.copy()
    C = deformed.get_cell().array  # rows = a,b,c
    C_new = C @ F.T                # <-- key fix (right-multiply by F^T)
    deformed.set_cell(C_new, scale_atoms=True)
    return deformed

def main():
    parser = argparse.ArgumentParser(description="ASE deformation via deformation gradient F (LAMMPS-compatible).")
    parser.add_argument("--epsilon", type=float, default=1e-6, help="Deformation magnitude (epsilon).")
    parser.add_argument("--xyz_ref", type=str, default=None, help="Path of nre reference xyz file")
    parser.add_argument("--case", type=str, default=None, help="Tag of the calculation")
    args = parser.parse_args()

    epsilon = args.epsilon
    xyz_ref = args.xyz_ref
    case_ = args.case

    print("epsilon: ", epsilon)
    print("xyz_ref:", xyz_ref)
    print("case:", case_)

    # Reference config

    atoms0 = read_xyz(xyz_ref, frame=0)
    #write_xyz(atoms0, xyz_ref+"_debug")

    C0 = atoms0.get_cell().array

    for dir_idx in range(1, 7):
        for sign in ('-', '+'):
            if sign == '-':
                s = -epsilon
            else: # sign = '+'
                s = epsilon

            F = deformation_gradient(dir_idx, s)
            atoms_def = apply_deformation_F(atoms0, F)

            idx = f"{dir_idx}{sign}"
            case_idx = "%s%s" %(case_, idx)

            out_data = f"{case_idx}.xyz"
            write_xyz(atoms_def, out_data)
            print(f"Wrote {out_data}")

if __name__ == "__main__":
    main()

