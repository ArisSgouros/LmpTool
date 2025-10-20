up="1.0e-6"
lmp_mpi -in in.elastic -v up $up > log.lammps
python parse_elastic.py --log "log.lammps" --epsilon $up
python ../../../src/stress2cij.py --epsilon $up
mv o.stiffness.json o.stiffness_this_$up.json
