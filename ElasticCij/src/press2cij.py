#!/usr/bin/env python3
import os
import json
import argparse
import numpy as np

def read_json(fname):
    with open(fname, "r") as f:
        data = json.load(f)
    return data

def main():
    parser = argparse.ArgumentParser(description="Build 6×6 stiffness matrix from ± press JSON files.")
    parser.add_argument("--prefix_in", default="", help="File prefix)")
    parser.add_argument("--prefix_out", default="press2cij", help="File prefix)")
    parser.add_argument("--components", type=str, default="xx,yy,zz,yz,xz,xy")
    parser.add_argument("--epsilon", type=str, default="")

    args = parser.parse_args()

    prefix_in = args.prefix_in
    prefix_out = args.prefix_out
    epsilon_list = args.epsilon.split(",")
    component_list = args.components.split(",")

    print(f"prefix_in: {prefix_in}")
    print(f"prefix_out: {prefix_out}")
    print(f"epsilon_list: {epsilon_list}")
    print(f"component_list: {component_list}")

    fname = f"{prefix_in}_ref.json"
    press_ref = read_json(fname)

    press_def = {}
    for epsilon in epsilon_list:
        press_def[epsilon] = {}
        for direction in component_list:
            fname = f"{prefix_in}_{direction}_{epsilon}.json"
            press_def[epsilon][direction] = read_json(fname)

    with open(f'{prefix_out}_press_all.json', 'w') as f: 
        json.dump(press_def, f, indent=2)

    stiffness_eps = {}
    for direction in component_list:
        for press in component_list:
            component = f"{direction}{press}"
            stiffness_eps[component] = []
            for epsilon in epsilon_list:
                epsilon_float = float(epsilon)
                stiffness_eps[component].append( -(press_def[epsilon][direction][press] - press_ref[press]) / epsilon_float )

    with open(f'{prefix_out}_stiffness_eps.json', 'w') as f: 
        json.dump(stiffness_eps, f, indent=2)

    # Calculate the average stiffness tensor
    stiffness_av = {}
    for component in stiffness_eps:
        stiffness_av[component] = np.mean(stiffness_eps[component])
    with open(f'{prefix_out}_stiffness_av.json', 'w') as f: 
        json.dump(stiffness_av, f, indent=2)

    stiffness_std = {}
    for component in stiffness_eps:
        values = np.array(stiffness_eps[component])
        stiffness_std[component] = np.std(values, ddof=1) / np.sqrt(len(values))
    with open(f'{prefix_out}_stiffness_std.json', 'w') as f: 
        json.dump(stiffness_std, f, indent=2)

if __name__ == "__main__":
    main()
