#!/usr/bin/env python3

from __future__ import annotations
import os, sys, argparse
import numpy as np
from pathlib import Path
from typing import Dict, Optional, Union, List
from ase.io import read, write
from ase.data import atomic_numbers
from ase import Atoms

# ------------------------- helpers -------------------------

def _parse_type_map(s: str) -> Dict[int, str]:
    """Parse '1=Si,2=O' -> {1:'Si', 2:'O'}; whitespace tolerated."""
    mapping: Dict[int, str] = {}
    if not s:
        return mapping
    for item in s.split(","):
        item = item.strip()
        if not item:
            continue
        k, v = item.split("=")
        mapping[int(k.strip())] = v.strip()
    return mapping


def read_lammps(input_file: str, style: str = "atomic", types: str = "") -> Atoms:
    """
    Read a LAMMPS data file.

    Parameters
    ----------
    input_file : str
        Path to LAMMPS data file.
    style : str
        LAMMPS atom style used in the file ('atomic', 'charge', 'molecular', 'full').
    types : str
        Optional mapping like "1=Si,2=O" to map integer types -> element symbols.

    Returns
    -------
    Atoms
    """
    in_path = Path(input_file)
    if not in_path.exists():
        raise FileNotFoundError(f"Input file not found: {in_path}")

    symbol_of_type = _parse_type_map(types)
    Z_of_type = {itype: atomic_numbers[sym] for itype, sym in symbol_of_type.items()} if symbol_of_type else None

    atoms = read(
        str(in_path),
        format="lammps-data",
        style=style,
        Z_of_type=Z_of_type
    )
    return atoms

def write_lammps(
    atoms: Atoms,
    output_file: str,
    style: str = "atomic",
    units: str = "metal",
    force_skew: bool = False,
    write_velocities: bool = True,
) -> Path:
    """
    Write a LAMMPS data file.

    Parameters
    ----------
    atoms : Atoms
        Structure to write.
    output_file : str
        Path to output .data
    style : str
        LAMMPS atom style ('atomic', 'charge', 'molecular', 'full').
    units : str
        LAMMPS units ('metal', 'real', 'si', 'lj', 'cgs', 'electron').
    force_skew : bool
        Allow triclinic (skewed) cell output if present.
    write_velocities : bool
        Include velocities if present.

    Returns
    -------
    Path
    """
    out_path = Path(output_file)
    has_vel = atoms.get_velocities() is not None

    kw = {
        "format": "lammps-data",
        "atom_style": style,
        "units": units,
        "force_skew": force_skew,
    }

    if has_vel and write_velocities:
        try:
            write(str(out_path), atoms, **kw, velocities=True)
        except TypeError:
            # Older ASE versions may not accept 'velocities' kwarg
            write(str(out_path), atoms, **kw)
    else:
        write(str(out_path), atoms, **kw)

    return out_path

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
    parser = argparse.ArgumentParser(description="ASE io examples")
    parser.add_argument("--file_in", default="", help="Input data file.")
    parser.add_argument("--file_out", default="", help="Output data file.")
    parser.add_argument("--type_in", default="xyz", help="Input data type.")
    parser.add_argument("--type_out", default="xyz", help="Output data type.")


    args = parser.parse_args()
    file_in: str = args.file_in
    type_in: str = args.type_in
    file_out: str = args.file_out
    type_out: str = args.type_out

    # Reference config
    if type_in == 'lammps':
        atoms0 = read_lammps(file_in, style="atomic")
    elif type_in == 'xyz':
        atoms0 = read_xyz(file_in, frame=0)
    else:
        print(f"Unknown type_in {type_in}")
        exit()

    if type_out == 'lammps':
        write_lammps(atoms0, file_out, style="atomic", units="metal", force_skew=True)
    elif type_out == 'xyz':
        write_xyz(atoms0, file_out)
    else:
        print(f"Unknown type_out {type_out}")
        exit()

if __name__ == "__main__":
    main()
