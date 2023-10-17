import sys
import copy as cp
import math as m

SET_LO_AT_ZERO = True
#
# A class containing the id of the type and coeffs.
#
class c_atom:
   def __init__(self, Id, molId, type, q, x, y, z):
      self.Id = Id
      self.molId = molId
      self.type = type
      self.q = q
      self.x = x
      self.y = y
      self.z = z

class c_bond:
   def __init__(self, Id, type, i, j):
      self.Id = Id
      self.type = type
      self.i = i
      self.j = j

class c_angle:
   def __init__(self, Id, type, i, j, k):
      self.Id = Id
      self.type = type
      self.i = i
      self.j = j
      self.k = k

class c_dihedral:
   def __init__(self, Id, type, i, j, k, l):
      self.Id = Id
      self.type = type
      self.i = i
      self.j = j
      self.k = k
      self.l = l
#
# This section parses the command line arguments
#
nargs = len(sys.argv) - 1
print("nargs = ", nargs)
ptype_offset = 0
btype_offset = 0
atype_offset = 0
dtype_offset = 0
iype_offset = 0
if nargs > 3:
   data_file_A = sys.argv[1]
   print("First data file:", data_file_A)
   data_file_B = sys.argv[2]
   print("Second data file:", data_file_B)
   data_file_merged = sys.argv[3]
   print("Merged data file:", data_file_merged)
   side = sys.argv[4]
   if side in ["-x", "+x", "-y", "+y", "-z", "+z", "0"]:
      print("Side:", side)
   else:
      print("wrong side value.", side)
      print("Choose between -x, +x, -y, +y, -z, +z, and 0")
      quit()
   if nargs > 4: ptype_offset = int(sys.argv[5])
   print("atom type offset: ",ptype_offset)
   if nargs > 5: btype_offset = int(sys.argv[6])
   print("bond type offset: ",btype_offset)
   if nargs > 6: atype_offset = int(sys.argv[7])
   print("angle type offset: ",atype_offset)
   if nargs > 7: dtype_offset = int(sys.argv[8])
   print("dihedral type offset: ",dtype_offset)
   if nargs > 8: iype_offset = int(sys.argv[9])
   print("improper type offset: ",iype_offset)
else:
   print("Please set the name of the data file")
   print("by issuing the command:\n")
   print("python lmp_merge_datafile.py \"datafile_A\" \"datafile_B\" \"newdatafile\"\n")
   print("exiting...")
   exit()

atoms = {}
bonds = {}
angles = {}
dihedrals = {}
impropers = {}
bond_coeffs = {}
angle_coeffs = {}
dihedral_coeffs = {}
improper_coeffs = {}
pair_coeffs = {}
masses = {}

#
# Deal with the first data file
#
print("Reading the datafile:", data_file_A)

n_atoms_A = 0

n_atoms = 0
n_bonds = 0
n_angles = 0
n_dihedrals = 0
n_impropers = 0

n_atom_types = 0
n_bond_types = 0
n_angle_types = 0
n_dihedral_types = 0
n_improper_types = 0

Pair_coeffs_start = 0
Bond_coeffs_start = 0
Angle_coeffs_start = 0
Dihedral_coeffs_start = 0
Improper_coeffs_start = 0

Llo = [0] * 3
Lhi = [0] * 3
LL  = [0] * 3

# Open the first file
g = open(data_file_A, 'r')

lines = []
line_num = 0
while True:
   line = g.readline()
   # Check if this is the end of file; if yes break the loop.
   if line == '':
      break

   line_split = line.split()

   if "atoms" in line:
      n_atoms = int(line_split[0])
   if "bonds" in line:
      n_bonds = int(line_split[0])
   if "angles" in line:
      n_angles = int(line_split[0])
   if "dihedrals" in line:
      n_dihedrals = int(line_split[0])
   if "impropers" in line:
      n_impropers = int(line_split[0])
   if "atom types" in line:
      n_atom_types = int(line_split[0])
   if "bond types" in line:
      n_bond_types = int(line_split[0])
   if "angle types" in line:
      n_angle_types = int(line_split[0])
   if "dihedral types" in line:
      n_dihedral_types = int(line_split[0])
   if "improper types" in line:
      n_improper_types = int(line_split[0])

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
try:
   print("bond coeffs:", n_bond_types)
   for ii in range(n_bond_types):
      line = lines[Bond_coeffs_start+ii].split()
      Id = int(line[0])
      coeffs = " ".join(line[1:])
      bond_coeffs[Id] = coeffs
except:
   print("Problem with reading bond coeffs section/not existing")
#
# Read the angle coeffs
try:
   print("angle coeffs:", n_angle_types)
   for ii in range(n_angle_types):
      line = lines[Angle_coeffs_start+ii].split()
      Id = int(line[0])
      coeffs = " ".join(line[1:])
      angle_coeffs[Id] = coeffs
except:
   print("Problem with reading angle coeffs section/not existing")
#
# Read the dihedral coeffs
try:
   print("dihedral coeffs:", n_dihedral_types)
   for ii in range(n_dihedral_types):
      line = lines[Dihedral_coeffs_start+ii].split()
      Id = int(line[0])
      coeffs = " ".join(line[1:])
      dihedral_coeffs[Id] = coeffs
except:
   print("Problem with reading dihedral coeffs section/not existing")
#
# Read the improper coeffs
try:
   print("improper coeffs:", n_improper_types)
   for ii in range(n_improper_types):
      line = lines[Improper_coeffs_start+ii].split()
      Id = int(line[0])
      coeffs = " ".join(line[1:])
      improper_coeffs[Id] = coeffs
except:
   print("Problem with reading improper coeffs section/not existing")
#
# Read the pair coeffs
try:
   print("pair coeffs:", n_atom_types)
   for ii in range(n_atom_types):
      line = lines[Pair_coeffs_start+ii].split()
      Id = int(line[0])
      coeffs = " ".join(line[1:])
      pair_coeffs[Id] = coeffs
except:
   print("Problem with reading pair coeffs section/not existing")
#
# Read the Masses section
try:
   print("Masses section")
   for ii in range(n_atom_types):
      line = lines[Masses_start+ii].split()
      Id = int(line[0])
      mass = line[1]
      masses[Id] = mass
except:
   print("Problem with reading masses section/not existing")
#
# Read the Atoms section
molId_offset = 0
try:
   print("atom section:", n_atoms)
   for ii in range(n_atoms):
      line = lines[Atoms_start+ii].split()
      Id = int(line[0])
      molId = int(line[1])
      molId_offset = max(molId_offset, molId)
      type = int(line[2])
      q = float(line[3])
      x = float(line[4])
      y = float(line[5])
      z = float(line[6])
      ix, iy, iz = [0, 0, 0]
      try:
         ix = int(line[7])
         iy = int(line[8])
         iz = int(line[9])
         x += ix * LL[0]
         y += iy * LL[1]
         z += iz * LL[2]
      except:
         pass

      iatom = c_atom(Id, molId, type, q, x, y, z)
      atoms[Id] = iatom
except:
   print("Problem with reading atom section")
#
# Read the Bonds section
try:
   print("bonds:", n_bonds)
   for ii in range(n_bonds):
      line = lines[Bonds_start+ii].split()
      Id = int(line[0])
      type = int(line[1])
      i = float(line[2])
      j = float(line[3])

      ibond = c_bond(Id, type, i, j)
      bonds[Id] = ibond
except:
   print("Problem with reading bond section")
#
# Read the Angles section
try:
   print("angles:", n_angles)
   for ii in range(n_angles):
      line = lines[Angles_start+ii].split()
      Id = int(line[0])
      type = int(line[1])
      i = float(line[2])
      j = float(line[3])
      k = float(line[4])

      iangle = c_angle(Id, type, i, j, k)
      angles[Id] = iangle
except:
   print("Problem with reading angle section")
#
# Read the Dihedrals section
try:
   print("dihedrals:", n_impropers)
   for ii in range(n_dihedrals):
      line = lines[Dehidrals_start+ii].split()
      Id = int(line[0])
      type = int(line[1])
      i = float(line[2])
      j = float(line[3])
      k = float(line[4])
      l = float(line[5])

      idihedral = c_dihedral(Id, type, i, j, k, l)
      dihedralss[Id] = idihedral
except:
   print("Problem with reading dihedrals section")
#
# Read the Impropers section
try:
   print("impropers", n_impropers)
   for ii in range(n_impropers):
      line = lines[Impropers_start+ii].split()
      Id = int(line[0])
      type = int(line[1])
      i = float(line[2])
      j = float(line[3])
      k = float(line[4])
      l = float(line[5])

      iimpriper = c_dihedral(Id, type, i, j, k, l)
      impropers[Id] = iimproper
except:
   print("Problem with reading impropers section")
#
# Compute the offset of the Ids and the dimensions of the box
atId_offset = n_atoms
bId_offset = n_bonds
aId_offset = n_angles
dId_offset = n_dihedrals
iId_offset = n_impropers

Llo_A = cp.deepcopy(Llo)
Lhi_A = cp.deepcopy(Lhi)
LL_A = [Lhi_A[d] - Llo_A[d] for d in range(3)]

print("box A dimensions:", Llo_A[0], Lhi_A[0], "xlo xhi")
print("                 ", Llo_A[1], Lhi_A[1], "ylo yhi")
print("                 ", Llo_A[2], Lhi_A[2], "zlo zhi")

#
# Deal with the second data file
#
print("Reading the datafile:", data_file_B)

n_atoms = 0
n_bonds = 0
n_angles = 0
n_dihedrals = 0
n_impropers = 0

n_atom_types = 0
n_bond_types = 0
n_angle_types = 0
n_dihedral_types = 0
n_improper_types = 0

Pair_coeffs_start = 0
Bond_coeffs_start = 0
Angle_coeffs_start = 0
Dihedral_coeffs_start = 0
Improper_coeffs_start = 0

Llo = [0] * 3
Lhi = [0] * 3
LL  = [0] * 3

# Open the first file
g = open(data_file_B, 'r')

lines = []
line_num = 0
while True:
   line = g.readline()
   # Check if this is the end of file; if yes break the loop.
   if line == '':
      break

   line_split = line.split()

   if "atoms" in line:
      n_atoms = int(line_split[0])
   if "bonds" in line:
      n_bonds = int(line_split[0])
   if "angles" in line:
      n_angles = int(line_split[0])
   if "dihedrals" in line:
      n_dihedrals = int(line_split[0])
   if "impropers" in line:
      n_impropers = int(line_split[0])
   if "atom types" in line:
      n_atom_types = int(line_split[0])
   if "bond types" in line:
      n_bond_types = int(line_split[0])
   if "angle types" in line:
      n_angle_types = int(line_split[0])
   if "dihedral types" in line:
      n_dihedral_types = int(line_split[0])
   if "improper types" in line:
      n_improper_types = int(line_split[0])

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
# Compute the offset of the Ids and the dimensions of the box
Llo_B = cp.deepcopy(Llo)
Lhi_B = cp.deepcopy(Lhi)
LL_B = [Lhi_B[d] - Llo_B[d] for d in range(3)]

print("box B dimensions:", Llo_B[0], Lhi_B[0], "xlo xhi")
print("                 ", Llo_B[1], Lhi_B[1], "ylo yhi")
print("                 ", Llo_B[2], Lhi_B[2], "zlo zhi")

coord_offset = [0.0, 0.0, 0.0]

Llo_merged = cp.deepcopy(Llo_A)
Lhi_merged = cp.deepcopy(Lhi_A)

if side == "-x":
   coord_offset[0] = Llo_A[0] - Lhi_B[0]
   Llo_merged[0] = Llo_B[0] + coord_offset[0]
elif side == "+x":
   coord_offset[0] = Lhi_A[0] - Llo_B[0]
   Lhi_merged[0] = Lhi_B[0] + coord_offset[0]
elif side == "-y":
   coord_offset[1] = Llo_A[1] - Lhi_B[1]
   Llo_merged[1] = Llo_B[1] + coord_offset[1]
elif side == "+y":
   coord_offset[1] = Lhi_A[1] - Llo_B[1]
   Lhi_merged[1] = Lhi_B[1] + coord_offset[1]
elif side == "-z":
   coord_offset[2] = Llo_A[2] - Lhi_B[2]
   Llo_merged[2] = Llo_B[2] + coord_offset[2]
elif side == "+z":
   coord_offset[2] = Lhi_A[2] - Llo_B[2]
   Lhi_merged[2] = Lhi_B[2] + coord_offset[2]
elif side == "0":
   print("no coord offset")

LL_final = [Lhi_merged[d] - Llo_merged[d] for d in range(3)]

#
# Read the bond coeffs
try:
   print("bond coeffs:", n_bond_types)
   for ii in range(n_bond_types):
      line = lines[Bond_coeffs_start+ii].split()
      Id = int(line[0]) + btype_offset
      coeffs = " ".join(line[1:])
      bond_coeffs[Id] = coeffs
except:
   print("Problem with reading bond coeffs section/not existing")
#
# Read the angle coeffs
try:
   print("angle coeffs:", n_angle_types)
   for ii in range(n_angle_types):
      line = lines[Angle_coeffs_start+ii].split()
      Id = int(line[0]) + atype_offset
      coeffs = " ".join(line[1:])
      angle_coeffs[Id] = coeffs
except:
   print("Problem with reading angle coeffs section/not existing")
#
# Read the dihedral coeffs
try:
   print("dihedral coeffs:", n_dihedral_types)
   for ii in range(n_dihedral_types):
      line = lines[Dihedral_coeffs_start+ii].split()
      Id = int(line[0]) + dtype_offset
      coeffs = " ".join(line[1:])
      dihedral_coeffs[Id] = coeffs
except:
   print("Problem with reading dihedral coeffs section/not existing")
#
# Read the improper coeffs
try:
   print("improper coeffs:", n_improper_types)
   for ii in range(n_improper_types):
      line = lines[Improper_coeffs_start+ii].split()
      Id = int(line[0]) + itype_offset
      coeffs = " ".join(line[1:])
      improper_coeffs[Id] = coeffs
except:
   print("Problem with reading improper coeffs section/not existing")
#
# Read the pair coeffs
try:
   print("pair coeffs:", n_atom_types)
   for ii in range(n_atom_types):
      line = lines[Pair_coeffs_start+ii].split()
      Id = int(line[0]) + ptype_offset
      coeffs = " ".join(line[1:])
      pair_coeffs[Id] = coeffs
except:
   print("Problem with reading pair coeffs section/not existing")
#
# Read the Masses section
try:
   print("Masses section")
   for ii in range(n_atom_types):
      line = lines[Masses_start+ii].split()
      Id = int(line[0]) + ptype_offset
      mass = line[1]
      masses[Id] = mass
except:
   print("Problem with reading masses section/not existing")
#
# Read the Atoms section
try:
   print("atom section:", n_atoms)
   for ii in range(n_atoms):
      line = lines[Atoms_start+ii].split()
      Id = int(line[0]) + atId_offset
      molId = int(line[1]) + molId_offset
      type = int(line[2]) + ptype_offset
      q = float(line[3])
      x = float(line[4])
      y = float(line[5])
      z = float(line[6])
      ix, iy, iz = [0, 0, 0]
      try:
         ix = int(line[7])
         iy = int(line[8])
         iz = int(line[9])
         x += ix * LL[0]
         y += iy * LL[1]
         z += iz * LL[2]
      except:
         pass

      x += coord_offset[0]
      y += coord_offset[1]
      z += coord_offset[2]

      iatom = c_atom(Id, molId, type, q, x, y, z)
      atoms[Id] = iatom
except:
   print("Problem with reading atom section")
#
# Read the Bonds section
try:
   print("bonds:", n_bonds)
   for ii in range(n_bonds):
      line = lines[Bonds_start+ii].split()
      Id = int(line[0]) + bId_offset
      type = int(line[1]) + btype_offset
      i = float(line[2]) + atId_offset
      j = float(line[3]) + atId_offset

      ibond = c_bond(Id, type, i, j)
      bonds[Id] = ibond
except:
   print("Problem with reading bond section")
#
# Read the Angles section
try:
   print("angles:", n_angles)
   for ii in range(n_angles):
      line = lines[Angles_start+ii].split()
      Id = int(line[0]) + aId_offset
      type = int(line[1]) + atype_offset
      i = float(line[2]) + atId_offset
      j = float(line[3]) + atId_offset
      k = float(line[4]) + atId_offset

      iangle = c_angle(Id, type, i, j, k)
      angles[Id] = iangle
except:
   print("Problem with reading angle section")
#
# Read the Dihedrals section
try:
   print("dihedrals:", n_impropers)
   for ii in range(n_dihedrals):
      line = lines[Dehidrals_start+ii].split()
      Id = int(line[0]) + dId_offset
      type = int(line[1]) + dtype_offset
      i = float(line[2]) + atId_offset
      j = float(line[3]) + atId_offset
      k = float(line[4]) + atId_offset
      l = float(line[5]) + atId_offset

      idihedral = c_dihedral(Id, type, i, j, k, l)
      dihedralss[Id] = idihedral
except:
   print("Problem with reading dihedrals section")
#
# Read the Impropers section
try:
   print("impropers", n_impropers)
   for ii in range(n_impropers):
      line = lines[Impropers_start+ii].split()
      Id = int(line[0]) + iId_offset
      type = int(line[1]) + itype_offset
      i = float(line[2]) + atId_offset
      j = float(line[3]) + atId_offset
      k = float(line[4]) + atId_offset
      l = float(line[5]) + atId_offset

      iimpriper = c_dihedral(Id, type, i, j, k, l)
      impropers[Id] = iimproper
except:
   print("Problem with reading impropers section")
#
# Write the merged data file
#

print("Generate the merged datafile:", data_file_merged)

#
# Generate the new datafile
#
tol = 1e-5
if (side in ["-x", "+x"] and ( abs(LL_A[1]-LL_B[1]) > tol or abs(LL_A[2]-LL_B[2]) > tol)) or \
   (side in ["-y", "+y"] and ( abs(LL_A[0]-LL_B[0]) > tol or abs(LL_A[2]-LL_B[2]) > tol)) or \
   (side in ["-z", "+z"] and ( abs(LL_A[0]-LL_B[0]) > tol or abs(LL_A[1]-LL_B[1]) > tol)):
   print("\nWARNING: the box crossections along the normal direction are not equal!\n")

if SET_LO_AT_ZERO:
   shift_orig = cp.deepcopy(Llo_merged)
   Llo_merged = [Llo_merged[ii]-shift_orig[ii] for ii in range(3)]
   Lhi_merged = [Lhi_merged[ii]-shift_orig[ii] for ii in range(3)]
   for iat in atoms.values():
      iat.x -= shift_orig[0]
      iat.y -= shift_orig[1]
      iat.z -= shift_orig[2]


print("merged box dimensions:")
print("                      ", Llo_merged[0], Lhi_merged[0], "xlo xhi")
print("                      ", Llo_merged[1], Lhi_merged[1], "ylo yhi")
print("                      ", Llo_merged[2], Lhi_merged[2], "zlo zhi")



f = open(data_file_merged, 'w')

f.write("# A data file merged from: %s and %s\n" % (data_file_A, data_file_B))
f.write("#\n")
f.write('%d atoms\n' % len(atoms))
f.write('%d bonds\n' % len(bonds))
f.write('%d angles\n' % len(angles))
f.write('%d dihedrals\n' % len(dihedrals))
f.write('%d impropers\n' % len(impropers))
f.write('\n')
f.write('%d atom types\n' % len(masses))
f.write('%d bond types\n' % len(bond_coeffs))
f.write('%d angle types\n' % len(angle_coeffs))
f.write('%d dihedral types\n' % len(dihedral_coeffs))
f.write('%d improper types\n' % len(improper_coeffs))
f.write('\n')
f.write('%f %f xlo xhi\n' % (Llo_merged[0], Lhi_merged[0]))
f.write('%f %f ylo yhi\n' % (Llo_merged[1], Lhi_merged[1]))
f.write('%f %f zlo zhi\n' % (Llo_merged[2], Lhi_merged[2]))
f.write('\n')

try:
   f.write('Masses\n')
   f.write('\n')
   for key in masses:
      f.write("%d %s\n" % (key, masses[key]))
except:
   print("Problem with writting mass section")
#
#for ii in range(atom_types):
#   line = lines[Masses_start + ii]
f.write('\n')
f.write('Atoms\n')
f.write('\n')
for ii in atoms.values():
   f.write("%d %d %d %f %f %f %f\n" % ( ii.Id, ii.molId, ii.type, ii.q, ii.x, ii.y, ii.z ))

try:
   if bonds:
      f.write('\n')
      f.write('Bonds\n')
      f.write('\n')
      for ii in bonds.values():
         f.write("%d %d %d %d\n" % ( ii.Id, ii.type, ii.i, ii.j ))
except:
   print("Problem with writting bonds/not existing")

try:
   if angles:
      f.write('\n')
      f.write('Angles\n')
      f.write('\n')
      for ii in angles.values():
         f.write("%d %d %d %d %d\n" % ( ii.Id, ii.type, ii.i, ii.j, ii.k ))
except:
   print("Problem with writting angles/not existing")

try:
   if dihedrals:
      f.write('\n')
      f.write('Dihedrals\n')
      f.write('\n')
      for ii in dihedrals.values():
         f.write("%d %d %d %d %d %d\n" % ( ii.Id, ii.type, ii.i, ii.j, ii.k, ii.l ))
except:
   print("Problem with writting dihedrals/not existing")

try:
   if impropers:
      f.write('\n')
      f.write('Impropers\n')
      f.write('\n')
      for ii in impropers.values():
         f.write("%d %d %d %d %d %d\n" % ( ii.Id, ii.type, ii.i, ii.j, ii.k, ii.l ))
except:
   print("Problem with writting impropers/not existing")

try:
   if pair_coeffs:
      f.write('\n')
      f.write('Pair Coeffs\n')
      f.write('\n')
      for key in pair_coeffs:
         f.write("%d %s\n" % (key, pair_coeffs[key]))
except:
   print("Problem with writting pair coeffs/not existing")

try:
   if bond_coeffs:
      f.write('\n')
      f.write('Bond Coeffs\n')
      f.write('\n')
      for key in bond_coeffs:
         f.write("%d %s\n" % (key, bond_coeffs[key]))
except:
   print("Problem with writting bond coeffs/not existing")

try:
   if angle_coeffs:
      f.write('\n')
      f.write('Angle Coeffs\n')
      f.write('\n')
      for key in angle_coeffs:
         f.write("%d %s\n" % (key, angle_coeffs[key]))
except:
   print("Problem with writting angle coeffs/not existing")

try:
   if dihedral_coeffs:
      f.write('\n')
      f.write('Dihedral Coeffs\n')
      f.write('\n')
      for key in dihedral_coeffs:
         f.write("%d %s\n" % (key, dihedral_coeffs[key]))
except:
   print("Problem with writting dihedral coeffs/not existing")

try:
   if improper_coeffs:
      f.write('\n')
      f.write('Improper Coeffs\n')
      f.write('\n')
      for key in improper_coeffs:
         f.write("%d %s\n" % (key, improper_coeffs[key]))
except:
   print("Problem with writting improper coeffs/not existing")

f.close()


#
# Write a LAMMPS dump file
#

f = open(data_file_merged+".lammpstrj", 'w')
f.write("ITEM: TIMESTEP\n")
f.write("0\n")
f.write("ITEM: NUMBER OF ATOMS\n")
f.write("%d\n" %(len(atoms)))
f.write("ITEM: BOX BOUNDS pp pp pp\n")
f.write('%f %f\n' % (Llo_merged[0], Lhi_merged[0]))
f.write('%f %f\n' % (Llo_merged[1], Lhi_merged[1]))
f.write('%f %f\n' % (Llo_merged[2], Lhi_merged[2]))
f.write("ITEM: ATOMS id mol type x y z\n")
for ii in atoms.values():
   f.write("%d %d %d %f %f %f\n" % ( ii.Id, ii.molId, ii.type, ii.x, ii.y, ii.z ))
f.close()
