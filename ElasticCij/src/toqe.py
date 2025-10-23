#!/usr/bin/env python3
import os, sys, argparse
from typing import Dict, Optional, Union, List
import numpy as np
from pathlib import Path
from io_helper import replace_string_in_file, append_line_to_file
from ase.io import read, write
from ase.data import atomic_numbers
from ase import Atoms

def read_xyz(input_file: str, frame: Optional[int] = 0) -> Union[Atoms, List[Atoms]]:
    """
    Read an extended XYZ file.

    Parameters
    ----------
    input_file : str
        Path to .xyz (extended XYZ).
    frame : int | None
        - int (default 0): return that frame as a single Atoms.
        - None: return all frames as a list[Atoms].

    Returns
    -------
    Atoms or list[Atoms]
    """
    in_path = Path(input_file)
    if not in_path.exists():
        raise FileNotFoundError(f"Input file not found: {in_path}")

    if frame is None:
        return read(str(in_path), format="extxyz", index=":")
    else:
        return read(str(in_path), format="extxyz", index=frame)

def write_xyz(atoms_or_list: Union[Atoms, List[Atoms]], output_file: str) -> Path:
    """
    Write an extended XYZ file.

    Parameters
    ----------
    atoms_or_list : Atoms | list[Atoms]
        Single structure or a trajectory.
    output_file : str
        Path to output .xyz

    Returns
    -------
    Path
    """
    out_path = Path(output_file)
    write(str(out_path), atoms_or_list, format="extxyz")
    return out_path

def main():
    parser = argparse.ArgumentParser(description="ASE script for quantum epsresso files")
    parser.add_argument("--header", default="header.in", help="Quantum espresso header file")
    parser.add_argument("--xyz", default="", help="name of xyz file")
    parser.add_argument("--cmd", default="", help="Quantum espresso cmd; e.g., mpirun -np 6 pw.x -npool 1")
    args = parser.parse_args()

    header = args.header
    xyz_file = args.xyz
    cmd_first = args.cmd
    name = xyz_file.replace(".xyz", "")
    script_in = f"{name}.in"
    script_out = f"{name}.out"

    atoms_def = read_xyz(xyz_file, frame=0)
    #write_xyz(atoms_def, xyz_def+"_debug")

    if not os.path.isfile(header):
       print(f'Cannot find header file: {header}')
       sys.exit()
    os.system(f'cp {header} {name}.in')

    # Generate cell parameters
    append_line_to_file(script_in, "CELL_PARAMETERS {angstrom}")
    cell = atoms_def.get_cell().array  # rows: a, b, c
    for vec in cell:
        append_line_to_file(script_in, f"{vec[0]:.16f} {vec[1]:.16f} {vec[2]:.16f}")

    # Generate atomic positions
    append_line_to_file(script_in, "ATOMIC_POSITIONS {angstrom}")
    symbols = atoms_def.get_chemical_symbols()
    positions = atoms_def.get_positions()  # Cartesian Å
    for s, (x, y, z) in zip(symbols, positions):
        append_line_to_file(script_in, f"{s} {x:.16f} {y:.16f} {z:.16f}")

    if cmd_first:
        cmd_full = f"{cmd_first} -in {script_in} > {script_out}"
        print(cmd_full)
    #    os.system(cmd)

if __name__ == "__main__":
    main()
