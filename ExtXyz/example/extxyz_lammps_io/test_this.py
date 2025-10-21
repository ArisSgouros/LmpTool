import os
import sys
# Add the parent directory to sys.path
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../'))
sys.path.insert(0, parent_dir)

from ase_lammps_io import read_lammps, write_xyz, read_xyz, write_lammps

atoms = read_lammps("in.data", style="atomic")
write_xyz(atoms, "o.xyz")

atoms2 = read_xyz("o.xyz", frame=0)
write_lammps(atoms2, "o.data", style="atomic", units="metal")

