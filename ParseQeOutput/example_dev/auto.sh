#!/bin/bash

path_pwo="../../QuantumEspressoHelper/parse_pwo.py"
path_pwi="../../QuantumEspressoHelper/parse_pwi.py"
path_json2xyz="../../Convert/json2xyz_ase.py"
path_json2lammps="../../Convert/json2lammps_ase.py"
path_nep2ase_xyz="../../Convert/nep2ase_xyz.py"
path_ase2nep_xyz="../../Convert/ase2nep_xyz.py"

# parse single input files
python $path_pwi i.MoSSe_3x3_03_scf.in --pretty -o o.03_in.json

# parse single output files
python $path_pwo i.MoSSe_3x3_01_vcrelax.out --pretty -o o.01.json
python $path_pwo i.MoSSe_3x3_02_relax.out --pretty -o o.02.json
python $path_pwo i.MoSSe_3x3_03_scf.out --pretty -o o.03.json
python $path_pwo i.MoSSe_3x3_04_scf_no_force_no_stress.out --pretty -o o.04.json

# test addition of custom fields
python add_custom.py o.01.json "config_type" "01" "weight" "1.0" > o.01_add.json

# convert json to xyz
python $path_json2xyz o.01.json > o.01_ase.xyz
python $path_json2xyz o.02.json > o.02_ase.xyz
python $path_json2xyz o.03.json > o.03_ase.xyz
python $path_json2xyz o.03_in.json > o.03_in_ase.xyz
python $path_json2xyz o.04.json > o.04_ase.xyz

# convert ase to nep
python $path_ase2nep_xyz o.01_ase.xyz -o o.01_nep.xyz
python $path_ase2nep_xyz o.02_ase.xyz -o o.02_nep.xyz
python $path_ase2nep_xyz o.03_ase.xyz -o o.03_nep.xyz
python $path_ase2nep_xyz o.03_in_ase.xyz -o o.03_in_nep.xyz
python $path_ase2nep_xyz o.04_ase.xyz -o o.04_nep.xyz
