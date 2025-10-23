#!/bin/bash

# parse single input files
python ../parse_pwi.py MoSSe_3x3_scf.in --pretty -o MoSSe_3x3_scf.in.json

# parse several input files
python ../parse_pwi.py "./*.in" --pretty > MoSSe_3x3_scf.in.many.json

# parse single output files
python ../parse_pwo.py MoSSe_3x3_scf.out --pretty -o MoSSe_3x3_scf.out.json
python ../parse_pwo.py MoSSe_3x3_scf_no_force_stress.out --pretty -o MoSSe_3x3_scf_no_force_stress.out.json
python ../parse_pwo.py MoSSe_3x3_relax.out --pretty -o MoSSe_3x3_relax.out.json
python ../parse_pwo.py MoSSe_3x3_vcrelax.out --pretty -o MoSSe_3x3_vcrelax.out.json

# parse sevaral output files
python ../parse_pwo.py "./MoSSe_3x3_relax*.out" --pretty > MoSSe_3x3_scf.out.many.json
