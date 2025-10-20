lmp_mpi -in in.elastic > log.lammps
python parse_elastic.py
python ../../../src/stress2cij.py --epsilon 1e-6
