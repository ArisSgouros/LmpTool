#!/bin/bash
python ../rmv_dupl_coeff.py datafile.orig datafile.final_unique_pair_coeff
python ../rmv_dupl_coeff.py datafile.orig datafile.final_duplicate_pair_coeff -unique_pair_coeff=0
