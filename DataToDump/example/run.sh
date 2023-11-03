#!/bin/bash
python ../data_to_dump.py in.pos.dat full --fmt="id,type,mol,x,y,z" -dump_file="o.dump.lammpstrj"
python ../data_to_dump.py in.pos_im.dat full --fmt="id,type,mol,x,y,z" -dump_file="o.dump_im.lammpstrj"
