###############################################################################
# MIT License
#
# Copyright (c) 2023 ArisSgouros
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
###############################################################################

import sys
import os
import argparse
from itertools import islice

parser = argparse.ArgumentParser()
parser.add_argument('dump_file', type=str, help='Path of the lammps dump file')
parser.add_argument('n_atom', type=int, help='const number of atoms')
parser.add_argument('dummy_attrib', type=str, help='attributes of dummy atoms')

if __name__ == "__main__":
   args = parser.parse_args()
   dump_file    = args.dump_file
   n_atom_const = args.n_atom
   dummy_str    = args.dummy_attrib

   print(dump_file)
   print(n_atom_const)
   print(dummy_str)


   with open(dump_file, 'r') as infile:
      while True:

         # Read the first header line
         line = infile.readline()
         if (line == ""):
           break;

         # Write the first header line
         print(line, end='') # ITEM: TIMESTEP

         # Read the remaining 8 lines from header
         lines_header   = []
         for line in islice(infile, 8):
            lines_header.append(line)
         n_atom_frame = int(lines_header[2])

         # Write the remaining header lines
         print(lines_header[0],end='')
         print(lines_header[1],end='') # ITEM: NUMBER OF ATOMS
         print(n_atom_const)
         print(lines_header[3],end='') # ITEM: BOX BOUNDS
         print(lines_header[4],end='')
         print(lines_header[5],end='')
         print(lines_header[6],end='')
         print(lines_header[7],end='') # ITEM: ATOMS

         lines_atom   = []
         for line in islice(infile, n_atom_frame):
            lines_atom.append(line)

         id = 1
         for ii in range(min(n_atom_frame, n_atom_const)):
            print(lines_atom[ii], end='')
            #tmp = lines.split()
            #tmp[0] = id
            #for attrib in tmp:
            #   print(attrib, end='')
            id += 1

         # Add extra atoms so the total is n_atom_const
         for id in range(n_atom_frame+1, n_atom_const+1):
            print(id, dummy_str)
