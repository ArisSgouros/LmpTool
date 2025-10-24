#!/usr/bin/env python3
from __future__ import annotations
import argparse, json
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
from ase import Atoms
from ase.io import write

def _is_num_seq(x) -> bool:
    return isinstance(x, (list, tuple)) and all(isinstance(t, (int, float)) for t in x)

def _sanitize_info(info: dict) -> dict:
    out = {}
    for k, v in info.items():
        if _is_num_seq(v):
            arr = np.asarray(v, dtype=float)
            out[k] = arr
        else:
            out[k] = v
    return out

def build_atoms_from_json(
    data: Dict[str, Any],
    *,
    include_cell: bool = True,
    include_forces: bool = True,
    include_info: bool = True,
    info_keys: Optional[List[str]] = None,
    pbc_mode: str = "auto",   # "auto" | "on" | "off"
) -> Atoms:
    atoms_data = data.get("atoms_data")
    if not isinstance(atoms_data, list) or len(atoms_data) == 0:
        raise ValueError("JSON must contain a non-empty 'atoms_data' list.")

    symbols = [a["symbols"] for a in atoms_data]
    positions = np.array([a["positions"] for a in atoms_data], dtype=float)

    kwargs: Dict[str, Any] = {"symbols": symbols, "positions": positions}

    # Cell + PBC
    if include_cell and ("cell_data" in data) and data["cell_data"] is not None:
        cell_list = data["cell_data"]
        if not (isinstance(cell_list, list) and len(cell_list) == 9):
            raise ValueError("cell_data must be a list of 9 numbers.")
        cell = np.array(cell_list, dtype=float).reshape(3, 3)
        kwargs["cell"] = cell
        kwargs["pbc"] = (True, True, True) if pbc_mode != "off" else (False, False, False)
    else:
        kwargs["pbc"] = (True, True, True) if pbc_mode == "on" else (False, False, False)

    atoms = Atoms(**kwargs)

    # Optional per-atom forces
    if include_forces:
        have_all_forces = all(
            ("forces" in a) and (a["forces"] is not None) and (len(a["forces"]) == 3)
            for a in atoms_data
        )
        if have_all_forces:
            forces = np.array([a["forces"] for a in atoms_data], dtype=float)
            atoms.new_array("forces", forces)  # this becomes 'force' in EXTXYZ

    # Global info (sanitize lists -> np.ndarray; 9->6 for stress/virial)
    if include_info and isinstance(data.get("info"), dict):
        info_src = data["info"]
        if info_keys:
            info_src = {k: info_src[k] for k in info_keys if k in info_src}
        atoms.info.update(_sanitize_info(info_src))

    return atoms

def main() -> None:
    p = argparse.ArgumentParser(description="Convert JSON (atoms/cell/info) to EXTXYZ using ASE.")
    p.add_argument("input", help="Path to input JSON file.")
    p.add_argument("-o", "--output", help="Path to output .xyz (default: stdout).")
    p.add_argument("--no-cell", action="store_true", help="Exclude cell (no Lattice in header).")
    p.add_argument("--no-forces", action="store_true", help="Exclude per-atom forces array.")
    p.add_argument("--no-info", action="store_true", help="Exclude global info tags.")
    p.add_argument("--info-keys", nargs="+", help="Only include these info keys (e.g., energy stress virial).")
    p.add_argument("--pbc", choices=["auto", "on", "off"], default="auto",
                   help="Periodic boundary conditions. 'auto' uses PBC if cell is present.")
    args = p.parse_args()

    data = json.loads(Path(args.input).read_text(encoding="utf-8"))
    atoms = build_atoms_from_json(
        data,
        include_cell=not args.no_cell,
        include_forces=not args.no_forces,
        include_info=not args.no_info,
        info_keys=args.info_keys,
        pbc_mode=args.pbc,
    )

    if args.output:
        write(args.output, atoms, format="extxyz")
    else:
        from sys import stdout
        write(stdout, atoms, format="extxyz")

if __name__ == "__main__":
    main()

