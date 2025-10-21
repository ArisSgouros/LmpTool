#!/bin/bash
up="0.01"
case="o.demo_"
#python 1.pymake_extxyz.py --out_qe="MoS2.scf_0.000.out" --xyz_ref=$case"ref.xyz"
#python 2.pydeform.py --epsilon=$up --case=$case --xyz_ref=$case"ref.xyz"
#python 3.pyqe.py --case=$case --run
#python 4.pyparse_stress.py --case=$case
python ../../src/stress2cij.py --epsilon $up --prefix=$case
#mv o.stiffness.json o.stiffness_this_$up.json
