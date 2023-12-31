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
import argparse

parser = argparse.ArgumentParser(description='Swap coefficients in lammps data files')
parser.add_argument('data_file_orig', type=str, help='Path of the lammps data file')
parser.add_argument('data_file_new', type=str, help='Path of the output data file')
parser.add_argument('swap_coeff_file', type=str, help='Path of the coeff file')

def swap_coeffs(swap_file_lines, type_start_swap):
   type_swap = {}
   if type_start_swap:
      for ii in range(sys.maxsize):
         if type_start_swap+1+ii > len(swap_file_lines):
            break
         line = swap_file_lines[ type_start_swap + ii ]
         if line == "\n":
            break
         Id = line.split()[0]
         new_coeffs = line
         type_swap[Id] = line
   return type_swap

if __name__ == "__main__":
   args = parser.parse_args()
   data_file = args.data_file_orig
   data_file_swap = args.data_file_new
   swap_coeff_file = args.swap_coeff_file

   g = open(data_file, 'r')

   atoms = 0
   bonds = 0
   angles = 0
   dihedrals = 0
   impropers = 0

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
   # Read the coefficients to be replaces
   #
   g = open(swap_coeff_file, 'r')

   Masses_start_swap = 0
   Pair_coeffs_start_swap = 0
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
      if "Pair Coeffs" in line:
         Pair_coeffs_start_swap = line_num + 2
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

   pair_coeffs_swap = swap_coeffs(swap_file_lines, Pair_coeffs_start_swap)
   masses_swap = swap_coeffs(swap_file_lines, Masses_start_swap)
   bond_coeffs_swap = swap_coeffs(swap_file_lines, Bond_coeffs_start_swap)
   angle_coeffs_swap = swap_coeffs(swap_file_lines, Angle_coeffs_start_swap)
   dihedral_coeffs_swap = swap_coeffs(swap_file_lines, Dihedral_coeffs_start_swap)
   improper_coeffs_swap = swap_coeffs(swap_file_lines, Improper_coeffs_start_swap)
   #
   # Generate the new datafile
   #
   f = open(data_file_swap, 'w')

   f.write(datafile_lines[0])
   f.write(datafile_lines[1])
   f.write('%d atoms\n' % atoms)
   f.write('%d bonds\n' % bonds)
   f.write('%d angles\n' % angles)
   f.write('%d dihedrals\n' % dihedrals)
   f.write('%d impropers\n' % impropers)
   f.write('\n')
   f.write('%d atom types\n' % atom_types)
   f.write('%d bond types\n' % bond_types)
   f.write('%d angle types\n' % angle_types)
   f.write('%d dihedral types\n' % dihedral_types)
   f.write('%d improper types\n' % improper_types)
   f.write('\n')
   f.write('%s %s xlo xhi\n' % (xlo, xhi))
   f.write('%s %s ylo yhi\n' % (ylo, yhi))
   f.write('%s %s zlo zhi\n' % (zlo, zhi))

   f.write('\n')
   f.write('Masses\n')
   f.write('\n')
   for ii in range(atom_types):
      line = datafile_lines[Masses_start + ii]
      Id = line.split()[0]
      if Id in masses_swap.keys():
         f.write(masses_swap[Id])
      else:
         f.write(line)

   f.write('\n')
   f.write('Atoms\n')
   f.write('\n')
   for ii in range(atoms):
      f.write(datafile_lines[Atoms_start + ii])
   f.write('\n')
   f.write('Bonds\n')
   f.write('\n')
   for ii in range(bonds):
      f.write(datafile_lines[Bonds_start + ii])
   f.write('\n')
   f.write('Angles\n')
   f.write('\n')
   for ii in range(angles):
      f.write(datafile_lines[Angles_start + ii])
   f.write('\n')
   f.write('Dihedrals\n')
   f.write('\n')
   for ii in range(dihedrals):
      f.write(datafile_lines[Dihedrals_start + ii])
   f.write('\n')
   f.write('Impropers\n')
   f.write('\n')
   for ii in range(impropers):
      f.write(datafile_lines[Impropers_start + ii])

   f.write('\n')
   f.write('Pair Coeffs\n')
   f.write('\n')
   for ii in range(atom_types):
      line = datafile_lines[Pair_coeffs_start + ii]
      Id = line.split()[0]
      if Id in pair_coeffs_swap.keys():
         f.write(pair_coeffs_swap[Id])
      else:
         f.write(line)

   if bonds:
      f.write('\n')
      f.write('Bond Coeffs\n')
      f.write('\n')
      for ii in range(bond_types):
         line = datafile_lines[Bond_coeffs_start + ii]
         Id = line.split()[0]
         if Id in bond_coeffs_swap.keys():
               f.write(bond_coeffs_swap[Id])
         else:
            f.write(line)

   if angles:
      f.write('\n')
      f.write('Angle Coeffs\n')
      f.write('\n')
      for ii in range(angle_types):
         line = datafile_lines[Angle_coeffs_start + ii]
         Id = line.split()[0]
         if Id in angle_coeffs_swap.keys():
            f.write(angle_coeffs_swap[Id])
         else:
            f.write(line)

   if dihedrals:
      f.write('\n')
      f.write('Dihedral Coeffs\n')
      f.write('\n')
      for ii in range(dihedral_types):
         line = datafile_lines[Dihedral_coeffs_start + ii]
         Id = line.split()[0]
         if Id in dihedral_coeffs_swap.keys():
            f.write(dihedral_coeffs_swap[Id])
         else:
            f.write(line)

   if impropers:
      f.write('\n')
      f.write('Improper Coeffs\n')
      f.write('\n')
      for ii in range(improper_types):
         line = datafile_lines[Improper_coeffs_start + ii]
         Id = line.split()[0]
         if Id in improper_coeffs_swap.keys():
            f.write(improper_coeffs_swap[Id])
         else:
            f.write(line)



   f.close()
