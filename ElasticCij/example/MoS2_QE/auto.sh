#!/bin/bash

# finite deformation (same for all directions)
epsilon="0.01"

case="o.demo"

# generate a reference extxyz file
python pymake_extxyz.py --out_qe=MoS2.scf_0.000.out --xyz_ref="$case"_ref.xyz

# run deform script (cartesion notation)
python ../../src/deform.py --data_in "$case"_ref.xyz --prefix_out $case --epsilon "$epsilon,-$epsilon" --components xx,yy,xy --component_fmt cart --in_filetype xyz --out_filetype xyz

# generate the quantum espresso scripts
cmd="mpirun -np 6 pw.x -npool 1"
python ../../src/toqe.py --header header.in --xyz "$case"_ref.xyz --cmd "$cmd"

# deformed configuration
#for i in xx yy xy; do
#    python qe.py --header header.in --xyz "$case"_"$i"_-"$epsilon".xyz
#    python qe.py --header header.in --xyz "$case"_"$i"_"$epsilon".xyz
#done

# compute the stiffness matrix
#python ../../../src/press2cij.py --prefix_in $prefix_in --prefix_out "o.press2cij" --epsilon "$epsilon,-$epsilon" --components xx,yy,zz,yz,xz,xy > o.log_press2cij

# remove redudant files
#rm log.lammps

#python 3.pyqe.py --case=$case --run
#python 4.pyparse_stress.py --case=$case
#python ../../src/stress2cij.py --epsilon $up --prefix=$case
#mv o.stiffness.json o.stiffness_this_$up.json
