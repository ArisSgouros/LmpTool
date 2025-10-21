import os, sys, argparse

# Add the parent directory to sys.path
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../ParseQeOutput/'))
sys.path.insert(0, parent_dir)
from parse_qe_output import parse_qe_output

# Add the path of extxyz
extxyz_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../ExtXyz/'))
sys.path.insert(0, extxyz_dir)
from extxyz import write_xyz_file

out_qe = "MoS2.scf_0.000.out"
xyz_ref = "in.ref.xyz"
case_ = "ref"

parser = argparse.ArgumentParser(description="ASE deformation via deformation gradient F (LAMMPS-compatible).")
parser.add_argument("--out_qe", type=str, default=None, help="Path of relaxed QE output")
parser.add_argument("--xyz_ref", type=str, default=None, help="Path of nre reference xyz file")
args = parser.parse_args()
out_qe = args.out_qe
xyz_ref = args.xyz_ref

print()
print('####################################')
print('Parsing file {:s}'.format(out_qe))
frame_def, status = parse_qe_output(out_qe, xyzfmt = 'ase', debug = False)

print()
print('Display parsed properties:')
# Check whether the properties were parsed
for key in frame_def['cell_properties']:
   print(key, frame_def['cell_properties'][key])

if ('forces' in frame_def['atoms'][0] or 'force' in frame_def['atoms'][0]):
   print('found forces')

if 'virial' in frame_def['cell_properties']:
   print('found virial')

# set the default weight
frame_def['cell_properties']['config_type'] = case_
write_xyz_file(xyz_ref, [frame_def], status='w', xyzfmt = 'ase')
