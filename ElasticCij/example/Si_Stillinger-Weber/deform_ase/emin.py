#!/usr/bin/env python3
import os, sys

def main():
    for dir_idx in range(1, 7):
        for sign in (-1, +1):
            tag = f"{dir_idx}{'-' if sign < 0 else '+'}"
            os.system(f"lmp_mpi -in in.emin -v tag {tag} > o.log_{tag}")
            print(f"Run {tag}")

if __name__ == "__main__":
    main()

