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

parser = argparse.ArgumentParser(description='Fold sheet to cnt')
parser.add_argument('data_file_in', type=str, help='Path of the lammps data file')
parser.add_argument('-data_file_out', default='o.pos.data', type=str, help='Path of the lammps data file')
parser.add_argument('-dim_roll', type=str, default='x', help='x, y or z')
parser.add_argument('-dim_norm', type=str, default='z', help='x, y or z')
parser.add_argument('-atomtype', type=str, default='full', help='Set the atomtype in the Lammps data file (e.g., full, atomic, etc.)')
parser.add_argument('-adjust_norm', type=float, default='30.0', help='Adjust the normal direction of the box to avoid overlap')

dim2int = {"x":0, "y":1, "z":2}

if __name__ == "__main__":
  args = parser.parse_args()
  data_file_in = args.data_file_in
  data_file_out = args.data_file_out
  dim_roll = dim2int[args.dim_roll]
  dim_norm = dim2int[args.dim_norm]
  atomtype = args.atomtype
  adjust_norm = args.adjust_norm

  # Print
  print("data_in     : ", data_file_in)
  print("data_out    : ", data_file_out)
  print("dim_roll    : ", dim_roll)
  print("dim_norm    : ", dim_norm)
  print("atomtype    : ", atomtype)
  print("adjust_norm : 2*r + ", adjust_norm)

  # Deal with the format of the Lammps data file
  if atomtype == 'full':
    kType = 2
    kColX = 4
    kColY = 5
    kColZ = 6
  else:
    print('unsupported format of Lammps data files for the atom type ' + atomtype)
    sys.exit()

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
        line_box = line_i
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
  type_list = []
  for ii in range(n_atom):
     lsplit = lines[line_atom+ii].split()
     type_ = int(lsplit[kType])
     xx = float(lsplit[kColX])
     yy = float(lsplit[kColY])
     zz = float(lsplit[kColZ])
     rr_list.append([xx, yy, zz])
     type_list.append(type_)

  # calculate the reference normal position of the sheet
  norm_ref = 0.0
  for rr in rr_list:
    norm_ref += rr[dim_norm]
  norm_ref /= n_atom
  print("norm_ref: ", norm_ref)

  # calculate the length of the sheet along the folding direction
  # and the radius of the NP
  length = LL[dim_roll]
  print("length: ", length)
  radius = length / (2.0*np.pi)

  # apply the deformation
  rr_cnt_list = []
  for rr in rr_list:
    theta = rr[dim_roll]/radius
    dr_ref = rr[dim_norm]-norm_ref
    rr_cnt = rr
    radius_ref = radius + dr_ref
    rr_cnt[dim_norm] = radius_ref*np.cos(theta) + norm_ref
    rr_cnt[dim_roll] = radius_ref*np.sin(theta)
    rr_cnt_list.append(rr_cnt)

  if adjust_norm > 1e-5:
    box_dim_norm_new = 2.0*radius + adjust_norm
    scale_norm = box_dim_norm_new/LL[dim_norm]
    LL[dim_norm] *= scale_norm
    Llo[dim_norm] *= scale_norm
    Lhi[dim_norm] *= scale_norm


  # update data file lines

  lines[line_box + 0] = "%.9f %.9f xlo xhi\n" % (Llo[0], Lhi[0])
  lines[line_box + 1] = "%.9f %.9f ylo yhi\n" % (Llo[1], Lhi[1])
  lines[line_box + 2] = "%.9f %.9f zlo zhi\n" % (Llo[2], Lhi[2])

  ii = 0
  for rr in rr_cnt_list:
    line = lines[line_atom+ii]
    lsplit = line.split()
    lsplit[kColX] = str(rr[0])
    lsplit[kColY] = str(rr[1])
    lsplit[kColZ] = str(rr[2])
    jline =  " ".join(lsplit)
    lines[line_atom+ii] = jline + "\n"
    ii += 1

  # export data file
  with open(data_file_out, 'w') as g:
    for line in lines:
       g.write(line)

  # export .lammpstrj
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
  for ii in range(len(rr_cnt_list)):
    rr = rr_cnt_list[ii]
    type_ = type_list[ii]
    foo.write("%d %d %f %f %f\n" % (ii, type_, rr[0], rr[1], rr[2]))
    ii += 1
  foo.close()
