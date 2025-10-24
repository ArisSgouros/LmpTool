#!/usr/bin/env python3
import argparse
from pathlib import Path
import sys

def main():
    ap = argparse.ArgumentParser(description="Convert NEP to ASE xyz files.")
    ap.add_argument("input", help="Input file")
    ap.add_argument("-o", "--output", help="Write result to this file")
    ap.add_argument("--inplace", action="store_true", help="Overwrite the input file")
    args = ap.parse_args()

    in_path = Path(args.input)
    text = in_path.read_text(encoding="utf-8")

    replacements = {
        "lattice=": "Lattice=",
        "properties=": "Properties=",
        "force:": "forces:",
    }
    for old, new in replacements.items():
        text = text.replace(old, new)

    if args.inplace:
        in_path.write_text(text, encoding="utf-8")
    elif args.output:
        Path(args.output).write_text(text, encoding="utf-8")
    else:
        sys.stdout.write(text)

if __name__ == "__main__":
    main()
