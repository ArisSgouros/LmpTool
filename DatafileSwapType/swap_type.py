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

parser = argparse.ArgumentParser(description='Swap coefficients in lammps data files')
parser.add_argument('data_file_orig', type=str, help='Path of the lammps data file')
parser.add_argument('data_file_new', type=str, help='Path of the output data file')
parser.add_argument('swap_type_file', type=str, help='Path of the swap type file')

def swap_coeffs(swap_file_lines, type_start_swap):
   type_swap = {}
   if type_start_swap:
      for ii in range(sys.maxsize):
         if type_start_swap+1+ii > len(swap_file_lines):
            break
         line = swap_file_lines[ type_start_swap + ii ]
         if line == "\n":
            break
         ii = line.split()[0]
         jj = line.split()[1]
         type_swap[ii] = jj
   return type_swap

def ReplaceNthSubstring(string, substring, nn):
   tmp = string.split(" ")
   tmp[nn] = substring
   return " ".join(tmp)

if __name__ == "__main__":
   args = parser.parse_args()
   data_file = args.data_file_orig
   data_file_swap = args.data_file_new
   swap_type_file = args.swap_type_file

   g = open(data_file, 'r')

   atoms = None
   bonds = None
   angles = None
   dihedrals = None
   impropers = None
   atom_types = None
   bond_types = None
   angle_types = None
   dihedral_types = None
   improper_types = None
   Pair_coeffs_start = None
   Bond_coeffs_start = None
   Angle_coeffs_start = None
   Dihedral_coeffs_start = None
   Improper_coeffs_start = None
   Masses_start = None
   Atoms_start = None
   Bonds_start = None
   Angles_start = None
   Dihedrals_start = None
   Impropers_start = None

   datafile_lines = []
   line_num = 0
   while True:
      line = g.readline()
      # Check if this is the end of file; if yes break the loop.
      if line == '':
         break

      line_split = line.split()

      if "atoms" in line:
         atoms = int(line_split[0])
      if "bonds" in line:
         bonds = int(line_split[0])
      if "angles" in line:
         angles = int(line_split[0])
      if "dihedrals" in line:
         dihedrals = int(line_split[0])
      if "impropers" in line:
         impropers = int(line_split[0])
      if "atom types" in line:
         atom_types = int(line_split[0])
      if "bond types" in line:
         bond_types = int(line_split[0])
      if "angle types" in line:
         angle_types = int(line_split[0])
      if "dihedral types" in line:
         dihedral_types = int(line_split[0])
      if "improper types" in line:
         improper_types = int(line_split[0])

      if "xlo xhi" in line:
         xlo = line_split[0]
         xhi = line_split[1]
      if "ylo yhi" in line:
         ylo = line_split[0]
         yhi = line_split[1]
      if "zlo zhi" in line:
         zlo = line_split[0]
         zhi = line_split[1]

      if "Pair Coeffs" in line:
         Pair_coeffs_start = line_num + 2
      if "Bond Coeffs" in line:
         Bond_coeffs_start = line_num + 2
      if "Angle Coeffs" in line:
         Angle_coeffs_start = line_num + 2
      if "Dihedral Coeffs" in line:
         Dihedral_coeffs_start = line_num + 2
      if "Improper Coeffs" in line:
         Improper_coeffs_start = line_num + 2

      if "Masses" in line:
         Masses_start = line_num + 2
      if "Atoms" in line:
         Atoms_start = line_num + 2
      if "Bonds" in line:
         Bonds_start = line_num + 2
      if "Angles" in line:
         Angles_start = line_num + 2
      if "Dihedrals" in line:
         Dihedrals_start = line_num + 2
      if "Impropers" in line:
         Impropers_start = line_num + 2

      datafile_lines.append(line)
      line_num += 1

   g.close()
   #
   # Read the coefficients to be replaced
   #
   g = open(swap_type_file, 'r')

   Masses_start_swap = 0
   Bond_coeffs_start_swap = 0
   Angle_coeffs_start_swap = 0
   Dihedral_coeffs_start_swap = 0
   Improper_coeffs_start_swap = 0

   swap_file_lines = []
   line_num = 0
   while True:
      line = g.readline()
      # Check if this is the end of file; if yes break the loop.
      if line == '':
         break

      line_split = line.split()

      if "Masses" in line:
         Masses_start_swap = line_num + 2
      if "Bond Coeffs" in line:
         Bond_coeffs_start_swap = line_num + 2
      if "Angle Coeffs" in line:
         Angle_coeffs_start_swap = line_num + 2
      if "Dihedral Coeffs" in line:
         Dihedral_coeffs_start_swap = line_num + 2
      if "Improper Coeffs" in line:
         Improper_coeffs_start_swap = line_num + 2

      swap_file_lines.append(line)
      line_num += 1
   g.close()

   #pair_coeffs_swap = swap_coeffs(swap_file_lines, Pair_coeffs_start_swap)
   ptype_swap = swap_coeffs(swap_file_lines, Masses_start_swap)
   btype_swap = swap_coeffs(swap_file_lines, Bond_coeffs_start_swap)
   atype_swap = swap_coeffs(swap_file_lines, Angle_coeffs_start_swap)
   dtype_swap = swap_coeffs(swap_file_lines, Dihedral_coeffs_start_swap)
   itype_swap = swap_coeffs(swap_file_lines, Improper_coeffs_start_swap)

   #
   # Generate the new datafile
   #
   f = open(data_file_swap, 'w')

   f.write(datafile_lines[0])
   f.write(datafile_lines[1])
   f.write('%d atoms\n' % atoms)
   if bonds:
      f.write('%d bonds\n' % bonds)
   if angles:
      f.write('%d angles\n' % angles)
   if dihedrals:
      f.write('%d dihedrals\n' % dihedrals)
   if impropers:
      f.write('%d impropers\n' % impropers)
   f.write('\n')
   if atom_types:
      f.write('%d atom types\n' % atom_types)
   if bond_types:
      f.write('%d bond types\n' % bond_types)
   if angle_types:
      f.write('%d angle types\n' % angle_types)
   if dihedral_types:
      f.write('%d dihedral types\n' % dihedral_types)
   if improper_types:
      f.write('%d improper types\n' % improper_types)
   f.write('\n')
   f.write('%s %s xlo xhi\n' % (xlo, xhi))
   f.write('%s %s ylo yhi\n' % (ylo, yhi))
   f.write('%s %s zlo zhi\n' % (zlo, zhi))

   if Masses_start:
      f.write('\n')
      f.write('Masses\n')
      f.write('\n')
      for ii in range(atom_types):
         line = datafile_lines[Masses_start + ii]
         type_ = line.split()[0]
         if type_ in ptype_swap.keys():
            line = ReplaceNthSubstring(line, ptype_swap[type_], 0)
         f.write(line)
   if Atoms_start:
      f.write('\n')
      f.write('Atoms\n')
      f.write('\n')
      for ii in range(atoms):
         line = datafile_lines[Atoms_start + ii]
         nn = 2
         type_ = line.split()[nn]
         if type_ in ptype_swap.keys():
            line = ReplaceNthSubstring(line, ptype_swap[type_], nn)
         f.write(line)
   if Bonds_start:
      f.write('\n')
      f.write('Bonds\n')
      f.write('\n')
      for ii in range(bonds):
         line = datafile_lines[Bonds_start + ii]
         nn = 1
         type_ = line.split()[nn]
         if type_ in btype_swap.keys():
            line = ReplaceNthSubstring(line, btype_swap[type_], nn)
         f.write(line)
   if Angles_start:
      f.write('\n')
      f.write('Angles\n')
      f.write('\n')
      for ii in range(angles):
         line = datafile_lines[Angles_start + ii]
         nn = 1
         type_ = line.split()[nn]
         if type_ in atype_swap.keys():
            line = ReplaceNthSubstring(line, atype_swap[type_], nn)
         f.write(line)
   if Dihedrals_start:
      f.write('\n')
      f.write('Dihedrals\n')
      f.write('\n')
      for ii in range(dihedrals):
         line = datafile_lines[Dihedrals_start + ii]
         nn = 1
         type_ = line.split()[nn]
         if type_ in dtype_swap.keys():
            line = ReplaceNthSubstring(line, dtype_swap[type_], nn)
         f.write(line)
   if Impropers_start:
      f.write('\n')
      f.write('Impropers\n')
      f.write('\n')
      for ii in range(impropers):
         line = datafile_lines[Impropers_start + ii]
         nn = 1
         type_ = line.split()[nn]
         if type_ in itype_swap.keys():
            line = ReplaceNthSubstring(line, itype_swap[type_], nn)
         f.write(line)
   if Pair_coeffs_start:
      f.write('\n')
      f.write('Pair Coeffs\n')
      f.write('\n')
      for ii in range(atom_types):
         line = datafile_lines[Pair_coeffs_start + ii]
         nn = 0
         type_ = line.split()[nn]
         if type_ in ptype_swap.keys():
            line = ReplaceNthSubstring(line, ptype_swap[type_], nn)
         f.write(line)
   if Bond_coeffs_start:
      f.write('\n')
      f.write('Bond Coeffs\n')
      f.write('\n')
      for ii in range(bond_types):
         line = datafile_lines[Bond_coeffs_start + ii]
         nn = 0
         type_ = line.split()[nn]
         if type_ in btype_swap.keys():
            line = ReplaceNthSubstring(line, btype_swap[type_], nn)
         f.write(line)
   if Angle_coeffs_start:
      f.write('\n')
      f.write('Angle Coeffs\n')
      f.write('\n')
      for ii in range(angle_types):
         line = datafile_lines[Angle_coeffs_start + ii]
         nn = 0
         type_ = line.split()[nn]
         if type_ in atype_swap.keys():
            line = ReplaceNthSubstring(line, atype_swap[type_], nn)
         f.write(line)
   if Dihedral_coeffs_start:
      f.write('\n')
      f.write('Dihedral Coeffs\n')
      f.write('\n')
      for ii in range(dihedral_types):
         line = datafile_lines[Dihedral_coeffs_start + ii]
         nn = 0
         type_ = line.split()[nn]
         if type_ in dtype_swap.keys():
            line = ReplaceNthSubstring(line, dtype_swap[type_], nn)
         f.write(line)
   if Improper_coeffs_start:
      f.write('\n')
      f.write('Improper Coeffs\n')
      f.write('\n')
      for ii in range(improper_types):
         line = datafile_lines[Improper_coeffs_start + ii]
         nn = 0
         type_ = line.split()[nn]
         if type_ in itype_swap.keys():
            line = ReplaceNthSubstring(line, itype_swap[type_], nn)
         f.write(line)

   f.close()
