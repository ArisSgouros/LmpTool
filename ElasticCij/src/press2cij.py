#!/usr/bin/env python3
from __future__ import annotations
import json
import argparse
from pathlib import Path
from typing import Dict, List, Tuple
import numpy as np
import sys

def parse_csv_list(s: str) -> List[str]:
    # robust split: "a, b,,c" -> ["a","b","c"]
    return [x.strip() for x in s.split(",") if x.strip()]

def read_json(path: Path) -> dict:
    try:
        with path.open("r") as f:
            return json.load(f)
    except FileNotFoundError:
        raise FileNotFoundError(f"Missing file: {path}")
    except json.JSONDecodeError as e:
        raise ValueError(f"Invalid JSON in {path}: {e}")

def main():
    parser = argparse.ArgumentParser(
        description="Build stiffness matrix from press JSON files."
    )
    parser.add_argument("--prefix_in", default="", help="Input file prefix (no suffix).")
    parser.add_argument("--prefix_out", default="press2cij", help="Output file prefix.")
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

    prefix_in: str = args.prefix_in
    prefix_out: str = args.prefix_out

    # Parse and validate lists
    epsilon_list_str = parse_csv_list(args.epsilon)
    if not epsilon_list_str:
        print("Error: --epsilon is empty after parsing.", file=sys.stderr)
        sys.exit(2)

    components = parse_csv_list(args.components)
    if not components:
        print("Error: --components is empty after parsing.", file=sys.stderr)
        sys.exit(2)

    # Convert epsilon strings to floats (only for math); keep original strings for filenames
    try:
        epsilon_floats = [float(e) for e in epsilon_list_str]
    except ValueError as e:
        print(f"Error: could not parse --epsilon values as float: {e}", file=sys.stderr)
        sys.exit(2)
    if any(e == 0.0 for e in epsilon_floats):
        print("Error: epsilon values must be non-zero.", file=sys.stderr)
        sys.exit(2)

    # Echo
    print(f"prefix_in: {prefix_in}")
    print(f"prefix_out: {prefix_out}")
    print(f"epsilon_list: {epsilon_list_str}")
    print(f"component_list: {components}")

    # Load reference pressures
    ref_path = Path(f"{prefix_in}_ref.json")
    press_ref = read_json(ref_path)

    # Load deformed pressures: press_def[epsilon_str][direction] = dict of component->value
    press_def: Dict[str, Dict[str, Dict[str, float]]] = {}
    for eps_str in epsilon_list_str:
        press_def[eps_str] = {}
        for direction in components:
            path = Path(f"{prefix_in}_{direction}_{eps_str}.json")
            press_def[eps_str][direction] = read_json(path)

    # Compute stiffness per epsilon using finite difference:
    # C_{ij} = - (P_j(def) - P_j(ref)) / ε_i, where i=direction, j=press component
    stiffness_eps: Dict[str, List[float]] = {f"{i}{j}": [] for i in components for j in components}

    for eps_str, eps_val in zip(epsilon_list_str, epsilon_floats):
        for i in components:
            def_row = press_def[eps_str][i]          # stresses under strain in direction i
            for j in components:
                try:
                    sij_def = float(def_row[j])
                    sij_ref = float(press_ref[j])
                except KeyError as e:
                    print(f"Error: missing stress component {e} in files for '{i}' or ref.", file=sys.stderr)
                    sys.exit(2)
                except ValueError as e:
                    print(f"Error: non-numeric stress in JSON: {e}", file=sys.stderr)
                    sys.exit(2)
                cij = -(sij_def - sij_ref) / eps_val
                stiffness_eps[f"{i}{j}"].append(cij)

    # Statistics: mean and SEM per component
    stiffness_mean: Dict[str, float] = {}
    stiffness_sem: Dict[str, float] = {}

    for comp, values in stiffness_eps.items():
        arr = np.asarray(values, dtype=float)
        if arr.size == 0:
            print(f"Warning: no samples for {comp}", file=sys.stderr)
            stiffness_mean[comp] = float("nan")
            stiffness_sem[comp] = float("nan")
            continue
        stiffness_mean[comp] = float(np.mean(arr))
        sem = float(np.std(arr, ddof=1) / np.sqrt(arr.size)) if arr.size > 1 else 0.0
        stiffness_sem[comp] = sem

    Path(f"{prefix_out}_stiffness.json").write_text(json.dumps(stiffness_mean, indent=2))

    info = {
        "components": components,
        "epsilons": epsilon_list_str,
        "stiffness": {
            "per_epsilon": stiffness_eps,
            "mean": stiffness_mean,
            "sem": stiffness_sem,
        },
    }
    Path(f"{prefix_out}_info.json").write_text(json.dumps(info, indent=2))

if __name__ == "__main__":
    main()

