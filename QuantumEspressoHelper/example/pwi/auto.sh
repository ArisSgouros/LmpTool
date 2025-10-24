#!/bin/bash

pwi2json="../../pwi2json.py"
json2xyz="../../json2xyz.py"
#01a_MoSSe.scf_lat-ang.in 01b_MoSSe.scf_lat-ang.in 01c_MoSSe.scf_ang-lat.in 02a_MoSSe.relax.in 03a_MoSSe.vcrelax.in

# input files
for case_ in 01a_MoSSe.scf_lat-ang.in 01b_MoSSe.scf_lat-ang.in 01c_MoSSe.scf_ang-lat.in 02a_MoSSe.relax.in 03a_MoSSe.vcrelax.in ; do
    python $pwi2json $case_ --pretty -o $case_.json
    python $json2xyz $case_.json > $case_.xyz
done
python $pwi2json "./01*.in" --pretty > 01d_combine_a-c.json
python $json2xyz 01d_combine_a-c.json > 01d_combine_a-c.xyz
