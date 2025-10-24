import os, sys, argparse
import json

# Add the parent directory to sys.path
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../ParseQeOutput/'))
sys.path.insert(0, parent_dir)
from qe_to_json import qe_out_to_json

parser = argparse.ArgumentParser(description="ASE deformation via deformation gradient F (LAMMPS-compatible).")
parser.add_argument("--out_qe", type=str, default=None, help="Path of relaxed QE output")
args = parser.parse_args()
out_qe = args.out_qe

print('Parsing file {:s}'.format(out_qe))
data = qe_out_to_json(out_qe)
with open(f"{out_qe}.json", "w") as f:
    json.dump(data, f, indent=4, sort_keys=True)
