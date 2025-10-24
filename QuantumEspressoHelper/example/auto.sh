#!/bin/bash

pwi2json="../pwi2json.py"
pwo2json="../pwo2json.py"
json2xyz="../json2xyz.py"

# input files
for case_ in 01a_MoSSe.scf_lat-ang.in 01b_MoSSe.scf_ang-lat.in 02a_MoSSe.relax.in; do
    python $pwi2json $case_ --pretty -o $case_.json
    python $json2xyz $case_.json > $case_.xyz
done
python ../pwi2json.py "./01*.in" --pretty > 01ba_in.combined.json
python $json2xyz 01ba_in.combined.json > 01ba_in.combined.xyz

# output files
for case_ in 01c_MoSSe.scf.out 02b_MoSSe.relax.out 03b_MoSSe.vcrelax_lat.out 03c_MoSSe.vcrelax_ang.out; do
    python $pwo2json $case_ --pretty -o $case_.json
    python $json2xyz $case_.json > $case_.xyz
done
python ../pwo2json.py "./03*.out" --pretty > 03cb_out.combined.json
python $json2xyz 03cb_out.combined.json > 03cb_out.combined.xyz

# partial xyz files
python $json2xyz 02b_MoSSe.relax.out.json --info-keys energy job_done --no-forces > 02b_MoSSe.relax.out.partial.json
