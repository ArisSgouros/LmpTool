#!/usr/bin/env python3
import os, sys
import numpy as np

# Make ../../../../ExtXyz/src importable (adjust if needed)
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../../ExtXyz/src/'))
sys.path.insert(0, parent_dir)

from ase_lammps_io import read_lammps, write_lammps

# ----- user knob -----
up = 1e-6   # deformation magnitude (engineering strain); negative applied below

def deformation_gradient(dir_idx: int, up: float) -> np.ndarray:
    """
    Build deformation gradient F for the *negative* deformation in 'elastic.in' style.
    Voigt dir: 1=xx, 2=yy, 3=zz, 4=yz, 5=xz, 6=xy
    """
    F = np.eye(3)
    s = up
    if dir_idx == 1:       # xx
        F[0, 0] += s
    elif dir_idx == 2:     # yy
        F[1, 1] += s
    elif dir_idx == 3:     # zz
        F[2, 2] += s
    elif dir_idx == 4:     # yz shear: c_y += s * lz  (equivalent to LAMMPS 'yz delta -up*lz0')
        F[1, 2] += s
    elif dir_idx == 5:     # xz shear: c_x += s * lz
        F[0, 2] += s
    elif dir_idx == 6:     # xy shear: b_x += s * ly
        F[0, 1] += s
    else:
        raise ValueError("dir_idx must be 1..6")
    return F

def apply_deformation(atoms, F: np.ndarray):
    """Return a copy with cell and positions affinely deformed: r' = F r, C' = F C."""
    deformed = atoms.copy()
    C = deformed.get_cell().array  # 3x3
    C_new = F @ C
    deformed.set_cell(C_new, scale_atoms=True)  # affine remap (like 'remap units box')
    return deformed

def main():
    # 0) read reference
    atoms0 = read_lammps("in.data", style="atomic")

    # 1) export undeformed for sanity (optional)
    write_lammps(atoms0, "o.pos_ref.data", style="atomic", units="metal")

    # 2) loop over Voigt directions (negative deformation only)
    for dir_idx in range(1, 7):
        for mode in [-1, 1]:
            tag = f"{dir_idx}"
            if mode == -1: tag += '-'
            if mode == 1 : tag += '+'
            F = deformation_gradient(dir_idx, up*mode)
            atoms_def = apply_deformation(atoms0, F)

            # Export LAMMPS data files
            out_data = f"o.pos_{tag}.data"
            write_lammps(atoms_def, out_data, style="atomic", units="metal")
            print(f"Wrote {out_data}")

            # Run LAMMPS
            os.system(f"lmp_mpi -in in.stress -v tag {tag} > o.log_{tag}")
            print(f"Eval stress {out_data}")

if __name__ == "__main__":
    main()

