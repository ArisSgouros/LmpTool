#!/bin/bash
ase_io="../../ase_io.py"
python $ase_io --file_in i.ase.xyz --type_in xyz --file_out o.xyz.xyz --type_out xyz
python $ase_io --file_in o.xyz.xyz --type_in xyz --file_out o.xyz.lammps --type_out lammps
python $ase_io --file_in o.xyz.lammps --type_in lammps --file_out o.lammps.xyz --type_out xyz

