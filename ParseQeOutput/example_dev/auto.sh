#!/bin/bash

path_pwo="../../QuantumEspressoHelper/parse_pwo.py"

# parse single output files
python $path_pwo i.MoSSe_3x3_01_vcrelax.out --pretty -o o.01.json
python $path_pwo i.MoSSe_3x3_02_relax.out --pretty -o o.02.json
python $path_pwo i.MoSSe_3x3_03_scf.out --pretty -o o.03.json
python $path_pwo i.MoSSe_3x3_04_scf_no_force_no_stress.out --pretty -o o.04.json
