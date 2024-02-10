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

parser = argparse.ArgumentParser(description='Decimate lammps dump files')
parser.add_argument('data_file_orig', type=str, help='Path of the original datafile')
parser.add_argument('data_file_final', type=str, help='Path of the final datafile')
parser.add_argument('-unique_pair_coeff', type=int, default=1, help='Generate the new lists containing the unique pair coefficients')

args = parser.parse_args()
data_file_orig  = args.data_file_orig
data_file_final = args.data_file_final
unique_pair_coeff    = args.unique_pair_coeff

#
# A class containing the id of the type and coeffs.
#
class c_coeff:
   def __init__(self, Id, coeffs):
      self.Id = Id
      self.coeffs = coeffs
      self.mass = ""
#
# This function removes the duplicate coefficients and returns
# a list with the unique coefficients plus a dictionary which
# returns the new unique id having as an input the old type id.
#
def rmv_duplicate_coeffs(type_coeffs_all):
   newId = 1
   type_coeffs_short = []
   type_coeffs_swap_id = {}
   for old in type_coeffs_all:
      firstTime = True
      for new in type_coeffs_short:
         if str(old.coeffs) == str(new.coeffs):
            firstTime = False
            break
      if firstTime:
         type_coeffs_short.append( c_coeff(newId, old.coeffs) )
         type_coeffs_swap_id[old.Id] = newId
         newId += 1
      else:
         type_coeffs_swap_id[old.Id] = new.Id

   return [type_coeffs_short, type_coeffs_swap_id]

g = open(data_file_orig, 'r')

atoms = 0
bonds = 0
angles = 0
dihedrals = 0
impropers = 0

atom_types = 0
bond_types = 0
angle_types = 0
dihedral_types = 0
improper_types = 0

lines = []
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

   lines.append(line)
   line_num += 1

g.close()

#
# Read the bond coeffs
#
bond_coeffs_all = [ 0 for ii in range(bond_types)]
for ii in range(bond_types):
   line = lines[Bond_coeffs_start+ii].split()
   Id = int(line[0])
   coeffs = []
   for jj in range(1,len(line)):
      icoeff = line[jj]
      if icoeff == '#':
         break
      else:
         coeffs.append(icoeff)
   bond_coeffs_all[ii] = c_coeff(Id, coeffs)
#
# Read the angle coeffs
#
angle_coeffs_all = [ 0 for ii in range(angle_types)]
for ii in range(angle_types):
   line = lines[Angle_coeffs_start+ii].split()
   Id = int(line[0])
   coeffs = []
   for jj in range(1,len(line)):
      icoeff = line[jj]
      if icoeff == '#':
         break
      else:
         coeffs.append(icoeff)
   angle_coeffs_all[ii] = c_coeff(Id, coeffs)
#
# Read the dihedral coeffs
#
dihedral_coeffs_all = [ 0 for ii in range(dihedral_types)]
for ii in range(dihedral_types):
   line = lines[Dihedral_coeffs_start+ii].split()
   Id = int(line[0])
   coeffs = []
   for jj in range(1,len(line)):
      icoeff = line[jj]
      if icoeff == '#':
         break
      else:
         coeffs.append(icoeff)
   dihedral_coeffs_all[ii] = c_coeff(Id, coeffs)
#
# Read the improper coeffs
#
improper_coeffs_all = [ 0 for ii in range(improper_types)]
for ii in range(improper_types):
   line = lines[Improper_coeffs_start+ii].split()
   Id = int(line[0])
   coeffs = []
   for jj in range(1,len(line)):
      icoeff = line[jj]
      if icoeff == '#':
         break
      else:
         coeffs.append(icoeff)
   improper_coeffs_all[ii] = c_coeff(Id, coeffs)
#
# Read the pair coeffs and masses
#
pair_coeffs_all = [ 0 for ii in range(atom_types)]
for ii in range(atom_types):
   line = lines[Pair_coeffs_start+ii].split()
   Id = int(line[0])
   coeffs = []
   for jj in range(1,len(line)):
      icoeff = line[jj]
      if icoeff == '#':
         break
      else:
         coeffs.append(icoeff)
   pair_coeffs_all[Id-1] = c_coeff(Id, coeffs)
for ii in range(atom_types):
   line = lines[Masses_start+ii].split()
   Id = int(line[0])
   mass = line[1]
   pair_coeffs_all[Id-1].mass = mass
#
# Generate the new lists containing the unique pair coefficients
#
if unique_pair_coeff:
   type_coeffs_all = pair_coeffs_all
   newId = 1
   type_coeffs_short = []
   type_coeffs_swap_id = {}
   for old in type_coeffs_all:
      firstTime = True
      for new in type_coeffs_short:
         if str(old.coeffs) == str(new.coeffs):
            firstTime = False
            if str(old.mass) != str(new.mass):
               print("ERROR: Atoms with same coeffs have different masses!")
               print("Please rerun with the unique_pair_coeff flag set to False..")
               print("Exiting..")
               sys.exit()
            break
      if firstTime:
         icoeff = c_coeff(newId, old.coeffs)
         icoeff.mass = old.mass
         type_coeffs_short.append( icoeff )
         type_coeffs_swap_id[old.Id] = newId
         newId += 1
      else:
         type_coeffs_swap_id[old.Id] = new.Id

   pair_coeffs_short = type_coeffs_short
   pair_coeffs_swap_id = type_coeffs_swap_id
else:
  type_coeffs_swap_id = {}
  for id in range(1,atom_types+1):
     type_coeffs_swap_id[id] = id
  [pair_coeffs_short, pair_coeffs_swap_id] = [pair_coeffs_all, type_coeffs_swap_id]

[bond_coeffs_short, bond_coeffs_swap_id] = rmv_duplicate_coeffs(bond_coeffs_all)
[angle_coeffs_short, angle_coeffs_swap_id] = rmv_duplicate_coeffs(angle_coeffs_all)
[dihedral_coeffs_short, dihedral_coeffs_swap_id] = rmv_duplicate_coeffs(dihedral_coeffs_all)
[improper_coeffs_short, improper_coeffs_swap_id] = rmv_duplicate_coeffs(improper_coeffs_all)

atom_types = len(pair_coeffs_short)
bond_types = len(bond_coeffs_short)
angle_types = len(angle_coeffs_short)
dihedral_types = len(dihedral_coeffs_short)
improper_types = len(improper_coeffs_short)

#for kk in pair_coeffs_short: print(kk.Id, kk.coeffs)
#for kk in bond_coeffs_short: print(kk.Id, kk.coeffs)
#for kk in angle_coeffs_short: print(kk.Id, kk.coeffs)
#for kk in dihedral_coeffs_short: print(kk.Id, kk.coeffs)
#for kk in improper_coeffs_short: print(kk.Id, kk.coeffs)

#
# Generate the new datafile
#
f = open(data_file_final, 'w')

f.write(lines[0])
f.write(lines[1])
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
for ii in pair_coeffs_short:
   f.write("%d %s\n" % (ii.Id, ii.mass))
#
#for ii in range(atom_types):
#   line = lines[Masses_start + ii]
f.write('\n')
f.write('Atoms\n')
f.write('\n')
for ii in range(atoms):
   line_split = lines[Atoms_start + ii].split()
   id_mol = ' '.join(line_split[0:2])
   old_type = int(line_split[2])
   new_type = pair_coeffs_swap_id[old_type]

   rest = ' '.join(line_split[3::])
   f.write("%s %d %s\n" % ( id_mol, new_type, rest ))
if bonds:
   f.write('\n')
   f.write('Bonds\n')
   f.write('\n')
   for ii in range(bonds):
      line_split = lines[Bonds_start + ii].split()
      id = line_split[0]
      old_type = int(line_split[1])
      new_type = bond_coeffs_swap_id[old_type]
      coeffs = ' '.join(line_split[2::])
      f.write("%s %d %s\n" % ( id, new_type, coeffs ))
if angles:
   f.write('\n')
   f.write('Angles\n')
   f.write('\n')
   for ii in range(angles):
      line_split = lines[Angles_start + ii].split()
      id = line_split[0]
      old_type = int(line_split[1])
      new_type = angle_coeffs_swap_id[old_type]
      coeffs = ' '.join(line_split[2::])
      f.write("%s %d %s\n" % ( id, new_type, coeffs ))
if dihedrals:
   f.write('\n')
   f.write('Dihedrals\n')
   f.write('\n')
   for ii in range(dihedrals):
      line_split = lines[Dihedrals_start + ii].split()
      id = line_split[0]
      old_type = int(line_split[1])
      new_type = dihedral_coeffs_swap_id[old_type]
      coeffs = ' '.join(line_split[2::])
      f.write("%s %d %s\n" % ( id, new_type, coeffs ))
if impropers:
   f.write('\n')
   f.write('Impropers\n')
   f.write('\n')
   for ii in range(impropers):
      line_split = lines[Impropers_start + ii].split()
      id = line_split[0]
      old_type = int(line_split[1])
      new_type = improper_coeffs_swap_id[old_type]
      coeffs = ' '.join(line_split[2::])
      f.write("%s %d %s\n" % ( id, new_type, coeffs ))
f.write('\n')
f.write('Pair Coeffs\n')
f.write('\n')
for ii in pair_coeffs_short:
   id = ii.Id
   coeffs = ' '.join(ii.coeffs)
   f.write("%d %s\n" % (id, coeffs))
if bonds:
   f.write('\n')
   f.write('Bond Coeffs\n')
   f.write('\n')
   for ii in bond_coeffs_short:
      id = ii.Id
      coeffs = ' '.join(ii.coeffs)
      f.write("%d %s\n" % (id, coeffs))
if angles:
   f.write('\n')
   f.write('Angle Coeffs\n')
   f.write('\n')
   for ii in angle_coeffs_short:
      id = ii.Id
      coeffs = ' '.join(ii.coeffs)
      f.write("%d %s\n" % (id, coeffs))
if dihedrals:
   f.write('\n')
   f.write('Dihedral Coeffs\n')
   f.write('\n')
   for ii in dihedral_coeffs_short:
      id = ii.Id
      coeffs = ' '.join(ii.coeffs)
      f.write("%d %s\n" % (id, coeffs))
if impropers:
   f.write('\n')
   f.write('Improper Coeffs\n')
   f.write('\n')
   for ii in improper_coeffs_short:
      id = ii.Id
      coeffs = ' '.join(ii.coeffs)
      f.write("%d %s\n" % (id, coeffs))
f.close()
