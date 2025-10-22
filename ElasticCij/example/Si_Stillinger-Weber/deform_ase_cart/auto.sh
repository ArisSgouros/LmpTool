#!/bin/bash

# finite deformation (same for all directions)
epsilon="1.0e-6"

# run deform script (cartesion notation)
python ../../../src/deform.py --data_in in.data --prefix_out o.pos --epsilon "$epsilon,-$epsilon" --components xx,yy,zz,yz,xz,xy --component_fmt cart

# minimize the potential energy with lammps and extract the pressure tensor
PARA="lmp_mpi -in in.emin "
prefix_in=o.press

# reference configuration
$PARA -v datafile o.pos_ref.data > log.lammps
mv press.json o.press_ref.json

# deformed configuration
for i in xx yy zz yz xz xy; do
    tag="$i"_-"$epsilon"
    $PARA -v datafile o.pos_$tag.data > log.lammps
    mv press.json "$prefix_in"_$tag.json
    tag="$i"_"$epsilon"
    $PARA -v datafile o.pos_$tag.data > log.lammps
    mv press.json "$prefix_in"_$tag.json
done

# compute the stiffness matrix
python ../../../src/press2cij.py --prefix_in $prefix_in --prefix_out "o.press2cij" --epsilon "$epsilon,-$epsilon" --components xx,yy,zz,yz,xz,xy > o.log_press2cij

# remove redudant files
rm log.lammps
