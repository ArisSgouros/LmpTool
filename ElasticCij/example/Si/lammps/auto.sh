lmp_mpi -in in.elastic -v up 1.0e-6 > log.lammps
python parse_elastic.py
python ../../../src/stress2cij.py --epsilon 1e-6
