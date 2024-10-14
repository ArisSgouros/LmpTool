###############################################################################
# MIT License
#
# Copyright (c) 2024 ArisSgouros
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
import argparse

parser = argparse.ArgumentParser(description='Rescale lammps data file')
parser.add_argument('data_file', type=str, help='Path of the lammps data file')
parser.add_argument('-atomtype', type=str, default='full', help='specify the atomtype')
parser.add_argument('-xscale', type=float, default=1.0, help='xscale')
parser.add_argument('-yscale', type=float, default=1.0, help='yscale')
parser.add_argument('-zscale', type=float, default=1.0, help='zscale')

if __name__ == "__main__":
   args = parser.parse_args()
   data_file = args.data_file
   xscale = args.xscale
   yscale = args.yscale
   zscale = args.zscale
   atomtype = args.atomtype

   if atomtype == 'full':
      col_x = 4
      col_y = 5
      col_z = 6
   else:
      print('Atom type %s not supported..' % (atomtype))
      sys.exit()

   with open(data_file, 'r') as foo:
      lines = foo.readlines()

   line_box = -1
   line_atoms = -1
   line_Atoms = -1
   count = 0
   for line in lines:
      if 'atoms' in line:
        line_atoms = count
      if 'xlo' in line and 'xhi' in line:
        line_box = count
      if 'Atoms' in line:
        line_Atoms = count + 2
      count += 1

   atoms = int(lines[line_atoms].split()[0])

   box_lo = [0.0]*3
   box_hi = [0.0]*3
   # fetch the original box bounds
   box_lo[0], box_hi[0], aux, aux = lines[line_box + 0].split()
   box_lo[1], box_hi[1], aux, aux = lines[line_box + 1].split()
   box_lo[2], box_hi[2], aux, aux = lines[line_box + 2].split()

   box_lo[0] = float(box_lo[0])*xscale
   box_hi[0] = float(box_hi[0])*xscale
   box_lo[1] = float(box_lo[1])*yscale
   box_hi[1] = float(box_hi[1])*yscale
   box_lo[2] = float(box_lo[2])*zscale
   box_hi[2] = float(box_hi[2])*zscale

   # replace box bounds with the transposed ones
   lines[line_box + 0] = "%.15f %.15f xlo xhi\n" % (box_lo[0], box_hi[0])
   lines[line_box + 1] = "%.15f %.15f ylo yhi\n" % (box_lo[1], box_hi[1])
   lines[line_box + 2] = "%.15f %.15f zlo zhi\n" % (box_lo[2], box_hi[2])
         
   # replace atom coordinates with the transposed ones
   for ii in range(atoms):
      line = lines[line_Atoms + ii]
      data = line.split()
      rr = [data[col_x], data[col_y], data[col_z]]
      # transpose
      data[col_x] = str(float(data[col_x])*xscale)
      data[col_y] = str(float(data[col_y])*yscale)
      data[col_z] = str(float(data[col_z])*zscale)
      line_scaled = ' '.join(data) + '\n'
      lines[line_Atoms + ii] = line_scaled

for line in lines:
   print(line, end='')

