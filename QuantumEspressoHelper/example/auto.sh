#!/bin/bash

pwi2json="../pwi2json.py"
pwo2json="../pwo2json.py"
json2xyz="../json2xyz.py"

# parse single input files
python $pwi2json MoSSe_3x3_scf.in --pretty -o MoSSe_3x3_scf.in.json

# parse several input files
python $pwi2json "./*.in" --pretty > MoSSe_3x3_scf.in.many.json

# parse single output files
python $pwo2json MoSSe_3x3_scf.out --pretty -o MoSSe_3x3_scf.out.json
python $pwo2json  MoSSe_3x3_scf_no_force_stress.out --pretty -o MoSSe_3x3_scf_no_force_stress.out.json
python $pwo2json  MoSSe_3x3_relax.out --pretty -o MoSSe_3x3_relax.out.json
python $pwo2json  MoSSe_3x3_vcrelax.out --pretty -o MoSSe_3x3_vcrelax.out.json

# parse sevaral output files
python $pwo2json  "./MoSSe_3x3_relax*.out" --pretty > MoSSe_3x3_scf.out.many.json

# convert to extxyz
python $json2xyz MoSSe_3x3_scf.in.json  > MoSSe_3x3_scf.in.xyz
python $json2xyz MoSSe_3x3_scf.out.json > MoSSe_3x3_scf.out.xyz
python $json2xyz MoSSe_3x3_scf.out.json --no-cell --no-forces --no-info  > MoSSe_3x3_scf_lite.out.xyz
python $json2xyz MoSSe_3x3_scf.out.json --info-keys energy virial > MoSSe_3x3_scf_custom-info.out.xyz
python $json2xyz MoSSe_3x3_scf_no_force_stress.out.json > MoSSe_3x3_scf_no_force_stress.out.xyz
