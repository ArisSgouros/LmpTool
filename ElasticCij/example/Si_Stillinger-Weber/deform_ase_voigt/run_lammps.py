#!/usr/bin/env python3
import os

def main():
    folder = "."  # current directory
    files = [f for f in os.listdir(folder) if f.endswith(".data")]

    first_parts = set()
    second_parts = set()

    for filename in files:
        if filename.count("_") < 2:
            continue  # skip files that don't have at least two underscores

        # Remove the ".data" suffix
        name = filename[:-5]

        # Split by "_"
        parts = name.split("_")

        component = parts[-2]   # the part before the last underscore
        epsilon = parts[-1]  # the part after the last underscore

        press_file = f"_{component}_{epsilon}.json"
        cmd = f"lmp_mpi -in in.emin -v datafile {filename} -v press_file {press_file} > log.lammps"
        print(cmd)
        os.system(cmd)

if __name__ == "__main__":
    main()

