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

parser = argparse.ArgumentParser(description='Decimate lammps dump files')
parser.add_argument('dump_file', type=str, help='Path of the lammps dump file')
parser.add_argument('column', type=int, help='Column to be sorted')

if __name__ == "__main__":
   args = parser.parse_args()
   dump_file = args.dump_file
   column = args.column - 1

   #
   # Start reading the frames and stripping the chosen types!
   #
   data_file = open(dump_file,'rt')
   while True:
      current_line = data_file.readline() #ITEM: TIMESTEP
      # Check if we reached end of file
      if current_line == '': break
      print( current_line, end = "")

      print( data_file.readline(), end ="")
      print( data_file.readline(), end ="") #ITEM: NUMBER OF ATOMS
      n_atom = int(data_file.readline())
      print( n_atom)
      print( data_file.readline(), end ="") #ITEM: BOX BOUNDS pp pp pp
      print( data_file.readline(), end ="")
      print( data_file.readline(), end ="")
      print( data_file.readline(), end ="")
      print( data_file.readline(), end ="") #ATOM STYLE FORMAT

      lines_sorted = {}
      for ii in range(n_atom):
         current_line = data_file.readline()
         current_line_split = current_line.split()
         col_val = int(current_line_split[column])

         lines_sorted[col_val] = current_line

      for ii in sorted(lines_sorted):
         print( lines_sorted[ii], end ="")
