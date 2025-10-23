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
def is_number(s: str) -> bool:
    try:
        float(s)
        return True
    except ValueError:
        return False

def parse_csv_list(s: str) -> List[str]:
    # robust split: "a, b,,c" -> ["a","b","c"]
    return [x.strip() for x in s.split(",") if x.strip()]

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
    C_new = C @ F.T                # right-multiply by F^T
    deformed.set_cell(C_new, scale_atoms=True)
    return deformed

def main():
    parser = argparse.ArgumentParser(description="ASE deformation via deformation gradient F (LAMMPS-compatible).")
    parser.add_argument("--data_in", default="", help="Input data file.")
    parser.add_argument("--prefix_out", default="deform", help="Output file prefix.")
    parser.add_argument("--component_fmt", default="voigt", help="cart or voigt")
    parser.add_argument("--in_filetype", default="lammps", help="lammps, xyz")
    parser.add_argument("--out_filetype", default="lammps", help="lammps, xyz")

    parser.add_argument(
        "--components",
        type=str,
        default="xx,yy,zz,yz,xz,xy",
        help="Comma-separated components (subset of xx,yy,zz,yz,xz,xy).",
    )
    parser.add_argument(
        "--epsilon",
        type=str,
        required=True,
        help="Comma-separated list of strain magnitudes as strings (e.g. '0.01,-0.01').",
    )

    args = parser.parse_args()
    data_in: str = args.data_in
    prefix_out: str = args.prefix_out
    component_fmt: str = args.component_fmt
    in_filetype: str = args.in_filetype
    out_filetype: str = args.out_filetype

    # Parse and validate lists
    epsilon_list_str = parse_csv_list(args.epsilon)
    if not epsilon_list_str:
        print("Error: --epsilon is empty after parsing.", file=sys.stderr)
        sys.exit(2)

    # Convert epsilon strings to floats; keep original strings for filenames
    try:
        epsilon_floats = [float(e) for e in epsilon_list_str]
    except ValueError as e:
        print(f"Error: could not parse --epsilon values as float: {e}", file=sys.stderr)
        sys.exit(2)
    if any(e == 0.0 for e in epsilon_floats):
        print("Error: epsilon values must be non-zero.", file=sys.stderr)
        sys.exit(2)

    components_list_str = parse_csv_list(args.components)
    if not components_list_str:
        print("Error: --components is empty after parsing.", file=sys.stderr)
        sys.exit(2)

    if component_fmt == "cart":
        # Convert components strings to Voigt if
        to_voigt = {"xx":1, "yy":2, "zz":3, "yz":4, "xz":5, "xy":6}
        try:
            components_ints = [to_voigt[c] for c in components_list_str]
        except ValueError as e:
            print(f"Error: could not convert cart to voigt", file=sys.stderr)
            sys.exit(2)
    elif component_fmt == "voigt":
        try:
            components_ints = [int(c) for c in components_list_str]
        except ValueError as c:
            print(f"Error: could not parse --component values as float: {c}", file=sys.stderr)
            sys.exit(2)
        if not all(1 <= x <= 6 for x in components_ints):
            print("Error: components values must be either [xx,yy,zz,yz,xz,xy] or [1,2,3,4,5,6]", file=sys.stderr)
            sys.exit(2)

    # Echo
    print(f"data_in: {data_in}")
    print(f"prefix_out: {prefix_out}")
    print(f"epsilon_list: {epsilon_list_str}")
    print(f"component_list (tag): {components_list_str}")
    print(f"component_list (voigt): {components_ints}")

    # Reference config
    if in_filetype == 'lammps':
        atoms0 = read_lammps(data_in, style="atomic")
        write_lammps(atoms0, f"{prefix_out}_ref.data", style="atomic", units="metal", force_skew=True)
    elif in_filetype == 'xyz':
        atoms0 = read_xyz(data_in, frame=0)
        write_xyz(atoms0, f"{prefix_out}_ref.xyz")

    for eps_str, eps_val in zip(epsilon_list_str, epsilon_floats):
        for comp_str, comp_int in zip(components_list_str, components_ints):
            F = deformation_gradient(comp_int, eps_val)
            atoms_def = apply_deformation_F(atoms0, F)
            if in_filetype == 'lammps':
                write_lammps(atoms_def, f"{prefix_out}_{comp_str}_{eps_str}.data", style="atomic", units="metal", force_skew=True)
            elif in_filetype == 'xyz':
                write_xyz(atoms_def, f"{prefix_out}_{comp_str}_{eps_str}.xyz")

if __name__ == "__main__":
    main()

