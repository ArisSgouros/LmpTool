#!/bin/bash
python ../../sheet_to_cnt.py in.sq10x20x1.dat o.sq10x20x1_xz.dat x z > o.sq10x20x1_xz.log

# displace along reference
python ../../sheet_to_cnt.py in.sq10x20x1_zref15.dat o.sq10x20x1_zref15_xz.dat x z > o.sq10x20x1_zref15_xz.log


python ../../sheet_to_cnt.py in.sq10x20x1.dat o.sq10x20x1_yz.dat y z > o.sq10x20x1_yz.log


#
python ../../sheet_to_cnt.py o.sq10x20x1_xz.dat o.sq10x20x1_xz_yx.dat y x > o.sq10x20x1_xz_yx.log

python ../../sheet_to_cnt.py in.sq1x10x20.dat o.sq1x10x20_yx.dat y x > o.sq1x10x20_yx.log
python ../../sheet_to_cnt.py in.sq1x10x20.dat o.sq1x10x20_zx.dat z x > o.sq1x10x20_zx.log

python ../../sheet_to_cnt.py in.200_4.dat o.200_4_xz.dat y z > o.200_4_xz.log
python ../../sheet_to_cnt.py o.200_4_xz.dat o.200_4_xz_yx.dat x y > o.200_4_xz_yx.log

