#!/usr/bin/env python3
"""
ase_lammps_io.py
================

Four simple I/O functions built on ASE:

    atoms_or_list = read_xyz(input_filename, frame=0)
    atoms          = read_lammps(input_filename, style="atomic", types="")
    write_xyz(atoms_or_list, output_filename)
    write_lammps(atoms, output_filename, style="atomic", units="metal",
                 force_skew=False, write_velocities=True)

Notes
-----
- `types` is an optional string mapping for LAMMPS integer types, e.g. "1=Si,2=O".
- `read_xyz(..., frame=None)` returns all frames (list[Atoms]); an integer frame returns a single Atoms.
- `write_lammps` will include velocities if present unless `write_velocities=False`. For older ASE, the
  `velocities` kwarg fallback is handled.
"""

from __future__ import annotations
from pathlib import Path
from typing import Dict, Optional, Union, List

from ase.io import read, write
from ase.data import atomic_numbers
from ase import Atoms


__all__ = [
    "read_xyz",
    "read_lammps",
    "write_xyz",
    "write_lammps",
]


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


# -------------------------- READ ---------------------------

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


# ------------------------- WRITE ---------------------------

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


# --------------------------- CLI (optional) ---------------------------

def _cli():
    import argparse
    p = argparse.ArgumentParser(description="Simple ASE-based I/O (XYZ ↔ LAMMPS data).")
    sub = p.add_subparsers(dest="cmd", required=True)

    # xyz -> data
    p_a = sub.add_parser("to-data", help="XYZ → LAMMPS data")
    p_a.add_argument("xyz", help="Input extended XYZ")
    p_a.add_argument("data", help="Output LAMMPS data")
    p_a.add_argument("--style", default="atomic",
                    choices=["atomic", "charge", "molecular", "full"])
    p_a.add_argument("--units", default="metal",
                    choices=["metal", "real", "si", "lj", "cgs", "electron"])
    p_a.add_argument("--force-skew", action="store_true")
    p_a.add_argument("--no-vel", action="store_true")
    p_a.add_argument("--frame", type=int, default=0,
                     help="Frame index to read from XYZ (use -1 for last, or run with --all then pick).")
    p_a.add_argument("--all", action="store_true",
                     help="Write all frames as separate DATA files (suffix .f{idx}.data).")

    # data -> xyz
    p_b = sub.add_parser("to-xyz", help="LAMMPS data → XYZ")
    p_b.add_argument("data", help="Input LAMMPS data")
    p_b.add_argument("xyz", help="Output XYZ")
    p_b.add_argument("--style", default="atomic",
                    choices=["atomic", "charge", "molecular", "full"])
    p_b.add_argument("--types", default="", help='Optional map "1=Si,2=O".')

    args = p.parse_args()

    if args.cmd == "to-data":
        if args.all:
            frames = read_xyz(args.xyz, frame=None)  # list[Atoms]
            for i, at in enumerate(frames):
                out = args.data if args.data.endswith(".data") else f"{args.data}.data"
                out_i = out.replace(".data", f".f{i}.data")
                write_lammps(
                    at, out_i, style=args.style, units=args.units,
                    force_skew=args.force_skew, write_velocities=not args.no_vel
                )
                print("Wrote", out_i)
        else:
            at = read_xyz(args.xyz, frame=args.frame)  # Atoms
            out = args.data if args.data.endswith(".data") else f"{args.data}.data"
            write_lammps(
                at, out, style=args.style, units=args.units,
                force_skew=args.force_skew, write_velocities=not args.no_vel
            )
            print("Wrote", out)

    elif args.cmd == "to-xyz":
        at = read_lammps(args.data, style=args.style, types=args.types)
        out = args.xyz if args.xyz.endswith(".xyz") else f"{args.xyz}.xyz"
        write_xyz(at, out)
        print("Wrote", out)


if __name__ == "__main__":
    _cli()

