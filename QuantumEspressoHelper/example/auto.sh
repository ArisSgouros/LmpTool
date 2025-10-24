#!/bin/bash

pwi2json="../pwi2json.py"
pwo2json="../pwo2json.py"
json2xyz="../json2xyz.py"

for case_ in 01a_MoSSe.scf_lat-ang.in 01b_MoSSe.scf_ang-lat.in 02a_MoSSe.relax.in; do
    python $pwi2json $case_ --pretty -o $case_.json
done

for case_ in 01c_MoSSe.scf.out 02b_MoSSe.relax.out 03b_MoSSe.vcrelax_lat.out 03c_MoSSe.vcrelax_ang.out; do
    python $pwo2json $case_ --pretty -o $case_.json
done

python ../pwo2json.py "./03*.out" --pretty > 03bc.out.combined.json
