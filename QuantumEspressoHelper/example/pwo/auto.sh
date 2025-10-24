#!/bin/bash

pwo2json="../../pwo2json.py"
json2xyz="../../json2xyz.py"

# output files
for case_ in 01c_MoSSe.scf.out 02b_MoSSe.relax.out 03b_MoSSe.vcrelax_lat.out 03c_MoSSe.vcrelax_ang.out; do
    python $pwo2json $case_ --pretty -o $case_.json
    python $json2xyz $case_.json > $case_.xyz
done
python $pwo2json "./03*.out" --pretty > 03d_combine_a-c.json
python $json2xyz 03d_combine_a-c.json > 03d_combine_a-c.xyz

# partial xyz files
python $json2xyz 02b_MoSSe.relax.out.json --info-keys energy job_done --no-forces > 02b_MoSSe.relax.out.partial.json
