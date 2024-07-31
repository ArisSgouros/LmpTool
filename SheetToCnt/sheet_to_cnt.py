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
import os
import argparse
import numpy as np

kColX = 4
kColY = 5
kColZ = 6

parser = argparse.ArgumentParser(description='Fold sheet to cnt')
parser.add_argument('data_file_in', type=str, help='Path of the lammps data file')
parser.add_argument('data_file_out', type=str, help='Path of the lammps data file')
parser.add_argument('side', type=str, help='x, y or z')
parser.add_argument('normal', type=str, help='x, y or z')

dim2int = {"x":0, "y":1, "z":2}

if __name__ == "__main__":
  args = parser.parse_args()
  data_file_in = args.data_file_in
  data_file_out = args.data_file_out
  side = dim2int[args.side]
  normal = dim2int[args.normal]

  print("data_in: ", data_file_in)
  print("data_out: ", data_file_out)
  print("side: ", side)
  print("normal: ", normal)

  # Read data file
  lines = []
  LL = [0.0, 0.0, 0.0]
  Llo = [0.0, 0.0, 0.0]
  Lhi = [0.0, 0.0, 0.0]
  line_atom = 0
  line_box = 0
  line_i = 0
  n_atom = 0
  with open(data_file_in, 'r') as g:
    while True:
      line = g.readline()
      # Check if this is the end of file; if yes break the loop.
      if line == '':
         break
 
      line_split = line.split()

      if "atoms" in line:
         n_atom = int(line_split[0])
      if "Atoms" in line:
         line_atom = line_i + 2

      if "xlo xhi" in line:
        Llo[0] = float(line_split[0])
        Lhi[0] = float(line_split[1])
        LL[0] = Lhi[0] - Llo[0]
      if "ylo yhi" in line:
        Llo[1] = float(line_split[0])
        Lhi[1] = float(line_split[1])
        LL[1] = Lhi[1] - Llo[1]
      if "zlo zhi" in line:
        Llo[2] = float(line_split[0])
        Lhi[2] = float(line_split[1])
        LL[2] = Lhi[2] - Llo[2]

      lines.append(line)
      line_i += 1

  # Read the atom coordinates
  rr_list = []
  for ii in range(n_atom):
     lsplit = lines[line_atom+ii].split()
     xx = float(lsplit[kColX])
     yy = float(lsplit[kColY])
     zz = float(lsplit[kColZ])
     
     rr_list.append([xx, yy, zz])

  # calculate the reference normal position of the sheet
  norm_ref = 0.0
  for rr in rr_list:
    norm_ref += rr[normal]
  norm_ref /= n_atom
  print("norm_ref: ", norm_ref)

  # calculate the length of the sheet along the folding direction
  length = LL[side]
  print("length: ", length)

  radius = length / (2.0*np.pi)

  rr_cnt_list = []
  for rr in rr_list:
    theta = rr[side]/radius
    rr_cnt = rr
    rr_cnt[normal] = radius*np.cos(theta) + norm_ref
    rr_cnt[side] = radius*np.sin(theta)
    rr_cnt_list.append(rr_cnt)

  # export data out
  ii = 0
  for rr in rr_cnt_list:
    line = lines[line_atom+ii]
    lsplit = line.split()
    lsplit[kColX] = str(rr[0])
    lsplit[kColY] = str(rr[1])
    lsplit[kColZ] = str(rr[2])
    jline =  " ".join(lsplit)
    print(jline)
    lines[line_atom+ii] = jline + "\n"
    ii += 1
    
  with open(data_file_out, 'w') as g:
    for line in lines:
       g.write(line)

  # export VMD
  foo = open(data_file_out+".lammpstrj", "w")
  foo.write("ITEM: TIMESTEP\n")
  foo.write("0\n")
  foo.write("ITEM: NUMBER OF ATOMS\n")
  foo.write("%d\n" %(n_atom))
  foo.write("ITEM: BOX BOUNDS pp pp pp\n")
  foo.write("%f %f\n" % (Llo[0], Lhi[0]))
  foo.write("%f %f\n" % (Llo[1], Lhi[1]))
  foo.write("%f %f\n" % (Llo[2], Lhi[2]))
  foo.write("ITEM: ATOMS id type xu yu zu\n")

  ii = 1
  type_ = 1
  for rr in rr_cnt_list:
    foo.write("%d %d %f %f %f\n" % (ii, type_, rr[0], rr[1], rr[2]))
    ii += 1
  foo.close()

