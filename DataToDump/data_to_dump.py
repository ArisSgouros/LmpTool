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
import copy as cp
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('data_file', type=str, help='Path of the lammps data file')
parser.add_argument('atom_style',  type=str, help='Lammps atom style')
parser.add_argument('-fmt', '--fmt', help='format of Lammps dump file', type=str)
parser.add_argument('-dump_file', type=str, default="dump.lammpstrj", help='Path of the Lammps dump file.')


if __name__ == "__main__":
   args = parser.parse_args()
   file_data = args.data_file
   atom_style = args.atom_style
   file_dump = args.dump_file
   fmt = [item for item in args.fmt.split(',')]

   print("data file  : ", file_data)
   print("atom_style : ", atom_style)
   print("dump_file  : ", file_dump)
   print("format     : ", fmt)

   atoms = {}

   n_atom = 0
   xlo = xhi = ylo = yhi = zlo = zhi = lx = ly = lz = -1

   col = {}
   if atom_style == "full":
      col["id"] = 0
      col["mol"] = 1
      col["type"] = 2
      col["q"] = 3
      col["x"] = 4
      col["y"] = 5
      col["z"] = 6
   elif atom_style == "atomic":
      col["id"] = 0
      col["type"] = 1
      col["x"] = 2
      col["y"] = 3
      col["z"] = 4
   else:
      print("Error: unsupported atom style ", atom_style)
      sys.exit()

   foo = open(file_data, 'r')
   lines = []
   line_num = 0
   while True:
      line = foo.readline()

      # Check if this is the end of file; if yes break the loop.
      if line == '':
         break

      line_split = line.split()
      if "atoms" in line:
         n_atom = int(line_split[0])
      if "xlo xhi" in line:
         xlo = line_split[0]
         xhi = line_split[1]
         lx = float(xhi) - float(xlo)
      if "ylo yhi" in line:
         ylo = line_split[0]
         yhi = line_split[1]
         ly = float(yhi) - float(ylo)
      if "zlo zhi" in line:
         zlo = line_split[0]
         zhi = line_split[1]
         lz = float(zhi) - float(zlo)
      if "Atoms" in line:
         Atoms_start = line_num + 2
      lines.append(line)
      line_num += 1
   foo.close()

   # Export the header of the dump file
   foo = open(file_dump, 'w')
   foo.write("ITEM: TIMESTEP\n")
   foo.write("0\n")
   foo.write("ITEM: NUMBER OF ATOMS\n")
   foo.write("%d\n" % (n_atom))
   foo.write("ITEM: BOX BOUNDS pp pp pp\n")
   foo.write("%s %s\n" % (xlo, xhi))
   foo.write("%s %s\n" % (ylo, yhi))
   foo.write("%s %s\n" % (zlo, zhi))
   foo.write("ITEM: ATOMS")
   for kind in fmt:
      foo.write(" %s" % (kind))
   foo.write("\n")

   fmt_no_coord = cp.copy(fmt)
   if "x" in fmt_no_coord:
      fmt_no_coord.remove("x")
   if "y" in fmt_no_coord:
      fmt_no_coord.remove("y")
   if "z" in fmt_no_coord:
      fmt_no_coord.remove("z")

   # Read the atom section of the data file while exporting the dump file
   for ii in range(n_atom):
      line = lines[Atoms_start+ii].split()
      atom = {}

      # deal with the coordinates
      atom["x"] = line[col["x"]]
      atom["y"] = line[col["y"]]
      atom["z"] = line[col["z"]]
      try:
         ix = int(line[col["x"] + 3])
         iy = int(line[col["y"] + 3])
         iz = int(line[col["z"] + 3])
         atom["x"] = str(float(atom["x"]) + ix*lx)
         atom["y"] = str(float(atom["y"]) + iy*ly)
         atom["z"] = str(float(atom["z"]) + iz*lz)
      except:
         pass
      # deal with the other attributes
      for kind in fmt_no_coord:
         atom[kind] = line[col[kind]]

      # export attributes to the dump file
      for kind in fmt:
         foo.write(atom[kind]+" ")
      foo.write("\n")
   foo.close()
