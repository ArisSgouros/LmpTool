#!/bin/bash
python ../data_to_dump.py in.pos.dat full -file_type="xyz" --fmt="type,x,y,z" -dump_file="o.dump.xyz"
python ../data_to_dump.py in.pos.dat full -file_type="lammpstrj" --fmt="id,type,mol,x,y,z" -dump_file="o.dump.lammpstrj"
python ../data_to_dump.py in.pos_im.dat full -file_type="lammpstrj" --fmt="id,type,mol,x,y,z" -dump_file="o.dump_im.lammpstrj"
