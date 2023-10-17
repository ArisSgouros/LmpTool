#!/bin/bash
python ../dump_const_atom.py dump.lammpstrj 1 "-1 0.0 0.0 0.0" > dump_max_1_orig000.lammpstrj
python ../dump_const_atom.py dump.lammpstrj 6 "-1 0.0 0.0 0.0" > dump_max_6_orig000.lammpstrj
python ../dump_const_atom.py dump.lammpstrj 6 "-1 5.0 5.0 5.0" > dump_max_6_orig555.lammpstrj
