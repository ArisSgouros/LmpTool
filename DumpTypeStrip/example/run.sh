#!/bin/bash
python ../type_strip.py dump.lammpstrj 2 -l "2,3" > dump_t1.lammpstrj
python ../type_strip.py dump.lammpstrj 2 -l "1"   > dump_t2_3.lammpstrj
