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

parser = argparse.ArgumentParser()
parser.add_argument('dump_file', type=str, help='Path of the lammps dump file')
parser.add_argument('type_col', type=int, help='column of the atom type')
parser.add_argument('-l', '--list', help='delimited list input', type=str)

if __name__ == "__main__":
   args = parser.parse_args()
   fileName = args.dump_file
   typeColumn = args.type_col
   strippedTypes = [int(item) for item in args.list.split(',')]

   newNumberOfAtoms = 0

   #
   #Read the first frame so an to compute the new number of atoms
   #
   data_file = open(fileName,'rt')

   data_file.readline() #ITEM: TIMESTEP
   data_file.readline()
   data_file.readline() #ITEM: NUMBER OF ATOMS
   oldNumberOfAtoms = int(data_file.readline())
   data_file.readline() #ITEM: BOX BOUNDS pp pp pp
   data_file.readline()
   data_file.readline()
   data_file.readline()
   data_file.readline() #ITEM: ATOMS id mol x y z
   for ii in range(0, oldNumberOfAtoms):
      currentLine = data_file.readline()
      currentLineSplit = currentLine.split()
      itype = int(currentLineSplit[typeColumn-1])
      if not itype in strippedTypes:
         newNumberOfAtoms += 1
   data_file.close()

   #
   # Start reading the frames and stripping the chosen types!
   #
   data_file = open(fileName,'rt')
   while True:
      current_line = data_file.readline() #ITEM: TIMESTEP
      # Check if we reached end of file
      if current_line == '': break
      print( current_line, end='')
      print( data_file.readline(), end='')
      print( data_file.readline(), end='') #ITEM: NUMBER OF ATOMS
      data_file.readline()
      print( newNumberOfAtoms)
      print( data_file.readline(), end='') #ITEM: BOX BOUNDS pp pp pp
      print( data_file.readline(), end='')
      print( data_file.readline(), end='')
      print( data_file.readline(), end='')
      print( data_file.readline(), end='') #ATOM STYLE FORMAT

      id = 0
      for ii in range(oldNumberOfAtoms):
         currentLine = data_file.readline()
         currentLineSplit = currentLine.split()
         itype = int(currentLineSplit[typeColumn-1])

         if not itype in strippedTypes:
            id += 1
            print( currentLine, end='')
   data_file.close()
