# finite deformation (same for all directions)
epsilon="1.0e-6"

# run lammps example in.elastic
lmp_mpi -in in.elastic -v up $epsilon > log.lammps

# parse the stiffness tensor from the lammps log file
python parse_elastic.py --log "log.lammps" > o.log_parse_elastic.py
mv o.stiffness.json o.lammps_stiffness.json

# compute the stiffness matrix using the custom code
python ../../../src/press2cij.py --prefix_in "o.press" --prefix_out "o.press2cij" --epsilon "$epsilon,-$epsilon" --components "xx,yy,zz,yz,xz,xy" > o.log_press2cij

# remove redudant files
rm restart.equil log.lammps
