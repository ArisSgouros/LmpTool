#!/bin/bash

# Roll a 10x20x1 xy-sheet along y
python ../../sheet_to_cnt.py in.sq10x20x1.dat -data_file_out=o.sq10x20x1_xz.dat -dim_roll=x -dim_norm=z > o.sq10x20x1_xz.log

# Roll a 10x20x1 xy-sheet (with z = 15) along y
python ../../sheet_to_cnt.py in.sq10x20x1_zref15.dat -data_file_out=o.sq10x20x1_zref15_xz.dat -dim_roll=x -dim_norm=z > o.sq10x20x1_zref15_xz.log

# Roll a 10x20x1 xy-sheet along x
python ../../sheet_to_cnt.py in.sq10x20x1.dat -data_file_out=o.sq10x20x1_yz.dat -dim_roll=y -dim_norm=z > o.sq10x20x1_yz.log

# Roll a y-oriented CNT along x
python ../../sheet_to_cnt.py o.sq10x20x1_xz.dat -data_file_out=o.sq10x20x1_xz_yx.dat -dim_roll=y -dim_norm=x > o.sq10x20x1_xz_yx.log

# Roll a 1x10x20 yz-sheet along z
python ../../sheet_to_cnt.py in.sq1x10x20.dat -data_file_out=o.sq1x10x20_yx.dat -dim_roll=y -dim_norm=x > o.sq1x10x20_yx.log

# Roll a 1x10x20 yz-sheet along y
python ../../sheet_to_cnt.py in.sq1x10x20.dat -data_file_out=o.sq1x10x20_zx.dat -dim_roll=z -dim_norm=x > o.sq1x10x20_zx.log

# Roll a 200x4x1 xy-sheet along y
python ../../sheet_to_cnt.py in.200_4.dat -data_file_out=o.200_4_xz.dat -dim_roll=y -dim_norm=z > o.200_4_xz.log

# Roll a y-oriented CNT along z
python ../../sheet_to_cnt.py o.200_4_xz.dat -data_file_out=o.200_4_xz_yx.dat -dim_roll=x -dim_norm=y > o.200_4_xz_yx.log

