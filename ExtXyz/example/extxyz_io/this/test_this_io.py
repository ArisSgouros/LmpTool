import os
import sys

# Add the parent directory to sys.path
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../'))
sys.path.insert(0, parent_dir)

from extxyz import read_xyz_file, write_xyz_file

inxyz = "../i.bi2se3_nep.xyz"
parsed_data = read_xyz_file(inxyz)

print(parsed_data)

write_xyz_file("o.nep_nep.xyz", parsed_data, status='w', xyzfmt = 'nep')
write_xyz_file("o.nep_ase.xyz", parsed_data, status='w', xyzfmt = 'ase')


inxyz = "../i.bi2se3_ase.xyz"
parsed_data = read_xyz_file(inxyz)

print(parsed_data)

write_xyz_file("o.ase_nep.xyz", parsed_data, status='w', xyzfmt = 'nep')
write_xyz_file("o.ase_ase.xyz", parsed_data, status='w', xyzfmt = 'ase')
