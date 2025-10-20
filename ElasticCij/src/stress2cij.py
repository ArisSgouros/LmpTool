#!/usr/bin/env python3
import os
import json
import argparse
import numpy as np

def load_stress(prefix: str, index: int, sign: str) -> np.ndarray:
    """Load a 6-component stress vector for a given deformation index and sign."""
    fname = f"{prefix}{index}{sign}.json"
    if not os.path.exists(fname):
        print(f"Error: {fname} not found.")
        quit()
    with open(fname, "r") as f:
        data = json.load(f)
    return np.array([float(data[str(i)]) for i in range(1, 7)], dtype=float)

def build_stress_matrix(prefix: str, sign: str) -> np.ndarray:
    """Build a 6x6 matrix with stress vectors as columns for given deformation sign."""
    return np.column_stack([load_stress(prefix, i, sign) for i in range(1, 7)])

def main():
    parser = argparse.ArgumentParser(description="Build 6×6 stiffness matrix from ± stress JSON files.")
    parser.add_argument("--prefix", default="o.stress_", help="File prefix (default: o.stress_)")
    parser.add_argument("--epsilon", type=float, default=1e-6, help="Strain magnitude used in deformation (default: 1e-6)")
    args = parser.parse_args()

    stress_plus  = build_stress_matrix(args.prefix, "+")
    stress_minus = build_stress_matrix(args.prefix, "-")

    stiffness = {}
    for i in range(6):
        for j in range(6):
            key = f"{i+1}{j+1}"
            stiffness[key] = (stress_minus[i][j] - stress_plus[i][j]) / (2.0 * args.epsilon)

    file_stiffness = 'o.stiffness.json'
    with open(file_stiffness, 'w') as f: 
        json.dump(stiffness, f, indent=2)

    print(f"Stiffness matrix written to {file_stiffness} (ε = {args.epsilon:g})")

if __name__ == "__main__":
    main()
