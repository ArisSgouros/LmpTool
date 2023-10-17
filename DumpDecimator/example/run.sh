#!/bin/bash
python ../dump_decimator.py dump.lammpstrj 1 > dump_ev1.lammpstrj
python ../dump_decimator.py dump.lammpstrj 2 > dump_ev2.lammpstrj
python ../dump_decimator.py dump.lammpstrj 5 > dump_ev5.lammpstrj
python ../dump_decimator.py dump.lammpstrj 10 > dump_ev10.lammpstrj
python ../dump_decimator.py dump.lammpstrj 11 > dump_ev11.lammpstrj
