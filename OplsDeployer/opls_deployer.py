import sys
from xml.dom import minidom
import xml.etree.ElementTree as ET
#
# A class containing the id of the type and coeffs.
#
class c_atom:
   def __init__(self):
      self.Id = -1
      self.molId = -1
      self.type = -1
      self.opls_name = ""
      self.opls_class = ""
      self.opls_element = ""
      self.opls_mass = -1.0
      self.opls_charge = ""
      self.default_charge = ""
      self.x = ""
      self.y = ""
      self.z = ""

class c_bond:
   def __init__(self):
      self.Id = -1
      self.type = -1
      self.i = -1
      self.j = -1

class c_angle:
   def __init__(self):
      self.Id = -1
      self.type = -1
      self.i = -1
      self.j = -1
      self.k = -1

class c_torsion:
   def __init__(self):
      self.Id = -1
      self.type = -1
      self.i = -1
      self.j = -1
      self.k = -1
      self.l = -1
##
## This section parses the command line arguments
##
if len(sys.argv) == 5:
   data_file = sys.argv[1]
   SMART_MOL_file = sys.argv[2]
   XMLFileName = sys.argv[3]
   new_data_file = sys.argv[4]
   print(data_file)
else:
   print("The format of the command arguments is the following:\n")
   print("python lmp_apply_OPLS_v2.py \"data_file\" \"SMART_MOL_file\" \"OPLS_PARAM_XML_FILE\" \"OUTPUT_NAME\" \n")
   print("example\n")
   print("python lmp_apply_OPLS_v2.py POS.SYSTEM.dat SMART_SYSTEM.dat oplsaa_from_foyer_group.xml out.data\n")
   print("exiting...")
   exit()

use_lmp_charges = False
generate_impropers = False
display_warnings = True

if use_lmp_charges:
   print("The charges from the LAMMPS datafile will be used..")
else:
   print("The charges from the OPLS force field will be used..")

if generate_impropers:
   print("Impropers will be generated..")
else:
   print("Impropers will not be generated..")
#
# Read the lammps data file
#

n_atoms = 0
n_bonds = 0
n_angles = 0
n_dihedrals = 0
n_impropers = 0

print("Reading LAMMPS data file \"" + data_file + "\"..")
g = open(data_file, 'r')

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

   if "xlo xhi" in line:
      xlo = line_split[0]
      xhi = line_split[1]
   if "ylo yhi" in line:
      ylo = line_split[0]
      yhi = line_split[1]
   if "zlo zhi" in line:
      zlo = line_split[0]
      zhi = line_split[1]

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

# reassign n_impropers to zero in case they are not needed
if not generate_impropers:
   n_impropers = 0

g.close()

atom_list = {}
bond_list = {}
angle_list = {}
dihedral_list = {}
improper_list = {}
#
# Read the atom section
print("   atoms section..")
for ii in range(n_atoms):
   line_split = lines[Atoms_start + ii].split()
   iat = c_atom()
   iat.Id = int(line_split[0])
   iat.molId = int(line_split[1])
   iat.default_charge = line_split[3]
   iat.x = line_split[4]
   iat.y = line_split[5]
   iat.z = line_split[6]
   atom_list[iat.Id] = iat
#
# Read the bond section
print("   bonds section..")
for ii in range(n_bonds):
   line_split = lines[Bonds_start + ii].split()
   ibond = c_bond()
   ibond.Id = int(line_split[0])
   ibond.i = float(line_split[2])
   ibond.j = float(line_split[3])
   bond_list[ibond.Id] = ibond
#
# Read the angle section
print("   angles section..")
for ii in range(n_angles):
   line_split = lines[Angles_start + ii].split()
   iangle = c_angle()
   iangle.Id = int(line_split[0])
   iangle.i = float(line_split[2])
   iangle.j = float(line_split[3])
   iangle.k = float(line_split[4])
   angle_list[iangle.Id] = iangle
#
# Read the dihedral section
print("   dihedrals section..")
for ii in range(n_dihedrals):
   line_split = lines[Dihedrals_start + ii].split()
   idihedral = c_torsion()
   idihedral.Id = int(line_split[0])
   idihedral.i = float(line_split[2])
   idihedral.j = float(line_split[3])
   idihedral.k = float(line_split[4])
   idihedral.l = float(line_split[5])
   dihedral_list[idihedral.Id] = idihedral
#
# Read the improper section
print("   impropers section..")
for ii in range(n_impropers):
   line_split = lines[Impropers_start + ii].split()
   iimproper = c_torsion()
   iimproper.Id = int(line_split[0])
   iimproper.i = float(line_split[2])
   iimproper.j = float(line_split[3])
   iimproper.k = float(line_split[4])
   iimproper.l = float(line_split[5])
   improper_list[iimproper.Id] = iimproper

#
# TEMPORAL SECTION.
# MUST BE MODIFIED IN ORDER TO WORK WITH MORE COMPLEX TOPOLOGIES
#
# Maybe we need a separate code that generates the OPLS parameters of the system
#
print("Reading the SMART topology files from", SMART_MOL_file, "..")

g = open(SMART_MOL_file, 'r')

for Id in range(1,n_atoms+1):
   line_split = g.readline().split()
   atom_list[Id].opls_name = line_split[3]
   atom_list[Id].opls_class = line_split[4]
   #print(Id, atom_list[Id].opls_name)

g.close()
#
# END OF TEMPORAL SECTION
#

#
# Read the OPLS atom types from the force field file
#
#XMLFileName = "oplsaa_from_foyer_group_mod.xml"
tree = ET.parse(XMLFileName)
root = tree.getroot()

child_AtomTypes = 0
child_HarmonicBondForce = 1
child_HarmonicAngleForce = 2
child_RBTorsionForce = 3
child_NonbondedForce = 4

#
# Read the atom attributes from the opls force field file
print("Matching the atom types from the \"" + XMLFileName + "\" force field file..")
def read_opls_atom_type(opls_name):
   for child in root[child_AtomTypes]:
      if child.attrib['name'] == opls_name:
         opls_element = child.attrib['element']
         opls_mass = child.attrib['mass']
#         opls_class = child.attrib['class']
#         return [opls_class, opls_mass, opls_element]
         print(opls_name, opls_mass, opls_element) # TEMP
         return [opls_mass, opls_element]
         break
   if display_warnings:
      print("   #Warning: type " + opls_name + " not found!!!")
   return [opls_name+"_mass", opls_name+"_element"]
for atom in atom_list.values():
   opls_name = atom.opls_name
   Id = atom.Id
   [atom_list[Id].opls_mass, atom_list[Id].opls_element] = read_opls_atom_type(opls_name)
#
# Read the opls coefficients
#
print("Retriving the opls coefficients..")
#
# Read the nonbonded coeffs from the opls force field file
print("   nonbonded coeffs..")
# Read the pair coeffs for each pair
def get_opls_atom_types(type):
   for child in root[child_NonbondedForce]:
      if child.attrib['type'] == type:
          return [child.attrib['epsilon'], child.attrib['sigma'], child.attrib['charge']]
   if display_warnings:
      print("      #Warning: Coeff for type " + type +" not found!")
   return [ type+"_eps", type+"_sig", 0.0 ]
#
atom_types = {}
n_atom_types = 0
type_Id = 1
for atom in atom_list.values():
   Id = atom.Id
   opls_name = atom_list[Id].opls_name
   [epsilon, sigma, opls_charge] = get_opls_atom_types(opls_name)
   coeffs = [epsilon, sigma, atom.opls_mass, atom.opls_element]
   # Check if the coefficient is unique.
   # If not, assign a previously defined type
   isCoeffUnique = True
   for itype in atom_types:
      if str(coeffs) == str(atom_types[itype]):
         type = itype
         isCoeffUnique = False
         break
   # If yes, define a new type
   if isCoeffUnique:
      type = type_Id
      atom_types[type_Id] = coeffs
      n_atom_types += 1
      type_Id += 1
   atom_list[Id].type = type
   atom_list[Id].opls_charge = opls_charge

# Read the bond coeffs for each bond
print("   bond coeffs..")
def get_opls_bond_coeffs(class1, class2):
   coeffs = []
   for child in root[child_HarmonicBondForce]:
      if (child.attrib['class1'] == class1 and child.attrib['class2'] == class2) or \
         (child.attrib['class1'] == class2 and child.attrib['class2'] == class1):
          coeffs.append([ child.attrib['k'], child.attrib['length'] ])
   n_coeffs = len(coeffs)
   if n_coeffs == 0 and display_warnings:
      coeffs.append([ class1+"-"+class2 ])
      print("      #Warning: Bond " + class1 +' '+ class2 + " not found!")
   elif n_coeffs == 1:
      coeffs = coeffs[0]
   elif n_coeffs > 1:
      coeffs.append([ class1+"-"+class2 ])
      print("      #Warning: Multiple (" + str(n_coeffs) + ") coeffs were found for " + class1 + " " + class2 + "!")

   return coeffs
#
bond_coeffs = {}
n_bond_types = 0
type_Id = 1
for bond in bond_list.values():
   Id = bond.Id
   class1 = atom_list[bond.i].opls_class
   class2 = atom_list[bond.j].opls_class
   coeff = get_opls_bond_coeffs(class1,class2)
   # Check if the coefficient is unique.
   # If not, assign a previously defined type
   isCoeffUnique = True
   for itype in bond_coeffs:
      if str(coeff) == str(bond_coeffs[itype]):
         type = itype
         isCoeffUnique = False
         break
   # If yes, define a new type
   if isCoeffUnique:
      type = type_Id
      bond_coeffs[type_Id] = coeff
      n_bond_types += 1
      type_Id += 1
   bond_list[Id].type = type
#
# Read the angle coeffs for each angle
print("   angle coeffs..")
def get_opls_angle_coeffs(class1, class2, class3):
   coeffs = []
   for child in root[child_HarmonicAngleForce]:
      if (child.attrib['class1'] == class1 and child.attrib['class2'] == class2 and child.attrib['class3'] == class3) or \
         (child.attrib['class1'] == class3 and child.attrib['class2'] == class2 and child.attrib['class3'] == class1):
          coeffs.append([ child.attrib['k'], child.attrib['angle'] ])
   n_coeffs = len(coeffs)
   if n_coeffs == 0 and display_warnings:
      coeffs.append([ class1+"-"+class2+"-"+class3 ])
      print("      #Warning: Angle " + class1 +'-'+ class2 + '-' + class3 + " not found!")
   elif n_coeffs == 1:
      coeffs = coeffs[0]
   elif n_coeffs > 1:
      coeffs.append([ class1+"-"+class2+"-"+class3 ])
      print("      #Warning: Multiple (" + str(n_coeffs) + ") coeffs were found for " + class1 + " " + class2 + " " + class3 + "!")
   return coeffs
#
angle_coeffs = {}
n_angle_types = 0
type_Id = 1
for angle in angle_list.values():
   Id = angle.Id
   class1 = atom_list[angle.i].opls_class
   class2 = atom_list[angle.j].opls_class
   class3 = atom_list[angle.k].opls_class
   coeff = get_opls_angle_coeffs(class1,class2,class3)
   # Check if the coefficient is unique.
   # If not, assign a previously defined type
   isCoeffUnique = True
   for itype in angle_coeffs:
      if str(coeff) == str(angle_coeffs[itype]):
         type = itype
         isCoeffUnique = False
         break
   # If yes, define a new type
   if isCoeffUnique:
      type = type_Id
      angle_coeffs[type_Id] = coeff
      n_angle_types += 1
      type_Id += 1
   angle_list[Id].type = type
#
# Read the dihedral coeffs for each dihedral
print("   dihedral coeffs..")
def get_opls_dihedral_coeffs(class1, class2, class3, class4):
   coeffs = []
   for child in root[child_RBTorsionForce]:
      if (child.attrib['class1'] == class1 and child.attrib['class2'] == class2 and child.attrib['class3'] == class3 and child.attrib['class4'] == class4) or \
         (child.attrib['class1'] == class4 and child.attrib['class2'] == class3 and child.attrib['class3'] == class2 and child.attrib['class4'] == class1):
          coeffs.append([  child.attrib['c0'], child.attrib['c1'], child.attrib['c2'], child.attrib['c3'], child.attrib['c4'], child.attrib['c5']  ])
   n_coeffs = len(coeffs)
   if n_coeffs == 0 and display_warnings:
      coeffs.append([ class1+"-"+class2+"-"+class3+"-"+class4 ])
      print("      #Warning: Dihedral " + class1+'-'+class2+"-"+class3+"-"+class4+ " not found!")
   elif n_coeffs == 1:
      coeffs = coeffs[0]
   elif n_coeffs > 1:
      coeffs.append([ class1+"-"+class2+"-"+class3+"-"+class4 ])
      print("      #Warning: Multiple (" + str(n_coeffs) + ") coeffs were found for " + class1+"-"+class2+"-"+class3+"-"+class4 + "!")
   return coeffs
#
dihedral_coeffs = {}
n_dihedral_types = 0
type_Id = 1
for dihedral in dihedral_list.values():
   Id = dihedral.Id
   class1 = atom_list[dihedral.i].opls_class
   class2 = atom_list[dihedral.j].opls_class
   class3 = atom_list[dihedral.k].opls_class
   class4 = atom_list[dihedral.l].opls_class
   coeff = get_opls_dihedral_coeffs(class1,class2,class3,class4)
   # Check if the coefficient is unique.
   # If not, assign a previously defined type
   isCoeffUnique = True
   for itype in dihedral_coeffs:
      if str(coeff) == str(dihedral_coeffs[itype]):
         type = itype
         isCoeffUnique = False
         break
   # If yes, define a new type
   if isCoeffUnique:
      type = type_Id
      dihedral_coeffs[type_Id] = coeff
      n_dihedral_types += 1
      type_Id += 1
   dihedral_list[Id].type = type
#
# Read the improper coeffs for each improper
print("   improper coeffs..")
def get_opls_improper_coeffs(class1, class2, class3, class4):
   coeffs = []
   for child in root[child_RBTorsionForce]:
      if (child.attrib['class1'] == class1 and child.attrib['class2'] == class2 and child.attrib['class3'] == class3 and child.attrib['class4'] == class4) or \
         (child.attrib['class1'] == class4 and child.attrib['class2'] == class3 and child.attrib['class3'] == class2 and child.attrib['class4'] == class1):
          coeffs.append([  child.attrib['c0'], child.attrib['c1'], child.attrib['c2'], child.attrib['c3'], child.attrib['c4'], child.attrib['c5']  ])
   n_coeffs = len(coeffs)
   if n_coeffs == 0 and display_warnings:
      coeffs.append([ class1+"-"+class2+"-"+class3+"-"+class4 ])
      print("      #Warning: Improper " + class1+'-'+class2+"-"+class3+"-"+class4+ " not found!")
   elif n_coeffs == 1:
      coeffs = coeffs[0]
   elif n_coeffs > 1:
      coeffs.append([ class1+"-"+class2+"-"+class3+"-"+class4 ])
      print("      #Warning: Multiple (" + str(n_coeffs) + ") coeffs were found for " + class1+"-"+class2+"-"+class3+"-"+class4 + "!")
   return coeffs
#
improper_coeffs = {}
n_improper_types = 0
type_Id = 1
for improper in improper_list.values():
   Id = improper.Id
   class1 = atom_list[improper.i].opls_class
   class2 = atom_list[improper.j].opls_class
   class3 = atom_list[improper.k].opls_class
   class4 = atom_list[improper.l].opls_class
   coeff = get_opls_improper_coeffs(class1,class2,class3,class4)
   # Check if the coefficient is unique.
   # If not, assign a previously defined type
   isCoeffUnique = True
   for itype in improper_coeffs:
      if str(coeff) == str(improper_coeffs[itype]):
         type = itype
         isCoeffUnique = False
         break
   # If yes, define a new type
   if isCoeffUnique:
      type = type_Id
      improper_coeffs[type_Id] = coeff
      n_improper_types += 1
      type_Id += 1
   improper_list[Id].type = type
#
# Generate the new data_file
#
print("Generating the new lammps datafile with name " + new_data_file + "..")
f = open(new_data_file, 'w')

f.write(lines[0])
f.write(lines[1])
f.write('%d atoms\n' % n_atoms)
f.write('%d bonds\n' % n_bonds)
f.write('%d angles\n' % n_angles)
f.write('%d dihedrals\n' % n_dihedrals)
f.write('%d impropers\n' % n_impropers)
f.write('\n')
f.write('%d atom types\n' % n_atom_types)
f.write('%d bond types\n' % n_bond_types)
f.write('%d angle types\n' % n_angle_types)
f.write('%d dihedral types\n' % n_dihedral_types)
f.write('%d improper types\n' % n_improper_types)
f.write('\n')
f.write('%s %s xlo xhi\n' % (xlo, xhi))
f.write('%s %s ylo yhi\n' % (ylo, yhi))
f.write('%s %s zlo zhi\n' % (zlo, zhi))
f.write('\n')
f.write('Atoms\n')
f.write('\n')
for atom in atom_list.values():
   if use_lmp_charges:
      f.write("%d %d %d %s %s %s %s\n" % ( atom.Id, atom.molId, atom.type, atom.default_charge, atom.x, atom.y, atom.z ))
   else:
      f.write("%d %d %d %s %s %s %s\n" % ( atom.Id, atom.molId, atom.type, atom.opls_charge, atom.x, atom.y, atom.z ))

if n_bonds:
   f.write('\n')
   f.write('Bonds\n')
   f.write('\n')
   for bond in bond_list.values():
      f.write("%d %d %d %d\n" % ( bond.Id, bond.type, bond.i, bond.j ))
if n_angles:
   f.write('\n')
   f.write('Angles\n')
   f.write('\n')
   for angle in angle_list.values():
      f.write("%d %d %d %d %d\n" % ( angle.Id, angle.type, angle.i, angle.j, angle.k ))
if n_dihedrals:
   f.write('\n')
   f.write('Dihedrals\n')
   f.write('\n')
   for dihedral in dihedral_list.values():
      f.write("%d %d %d %d %d %d\n" % ( dihedral.Id, dihedral.type, dihedral.i, dihedral.j, dihedral.k, dihedral.l ))
if n_impropers:
   f.write('\n')
   f.write('Impropers\n')
   f.write('\n')
   for improper in improper_list.values():
      f.write("%d %d %d %d %d %d\n" % ( improper.Id, improper.type, improper.i, improper.j, improper.k, improper.l ))

f.write('\n')
f.write('Masses\n')
f.write('\n')
for type in atom_types:
   f.write("%s %s # %s\n" % (type, atom_types[type][2], atom_types[type][3]))
#
# Coefficients section and conversions from gromacs to lammps
#
# Regarding harmonic bonds, gmx: 2k(x - x0)^2,  while lammps: k(x-x0)^2
gmx2lmp_b_k = 0.001195 # 0.5 kj/mol/nm2 to kcal/mol/A2.
gmx2lmp_b_len = 10.0 # nm to A.
# Eg. CT-CT coeffs should change from 0.1529 nm and 224262.4 j/mol/nm^2 to 1.529 A
# and 268 kcal/mol/A^2 (see SI in OPLS paper)

gmx2lmp_a_k = 0.1195# 0.5 kj/mol to kcal/mol
gmx2lmp_a_ang = 57.29577951# rad to deg
# Eg. CT-CT-CT coeffs should change from th = 1.966986067 rad and k = 488.273" kj/mol/rad2,
# to 112.7 deg and 58.35 kcal/mol/rad2. (see SI in OPLS paper)

gmx2lmp_d_k = 0.239# kj/mol to kcal/mol

gmx2lmp_sig = 10.0 # nm to A.
gmx2lmp_eps = 0.239 # kj/mol to kcal/mol.
# Eg. CT-CT-CT-N coeffs should change from c0=5.48732, c1=0.02719, c2=0.0, c3=-5.51451, c4=0.0, c5=0.0
# (in kj/mol units) to the coeffs of the Ryckaert-Bellemans function (the procedure is discussed in gromacs
# manual-5.5.1.pdf, p. 81,82): F1=1.964, F2 = 0.0, F3=0.659, F4=0.0 (see 10.1021/ja9621760, p.10)
# Also, note that there are discrepancies regarding the C-C-C-C, H-C-C-H for hydrocarbons between gromax
# database and the original papers [10.1021/ja9621760].
# This is still a mystery to me and others, but since the parameters remain unchanged over the last 15 years
# they have to be correct...
# https://mailman-1.sys.kth.se/pipermail/gromacs.org_gmx-users/2012-March/068930.html
# https://mailman-1.sys.kth.se/pipermail/gromacs.org_gmx-users/2010-May/050780.html
f.write('\n')
f.write('Pair Coeffs\n')
f.write('\n')
for type in atom_types:
   try:
      lmp_eps = float(atom_types[type][0]) * gmx2lmp_eps
      lmp_sig = float(atom_types[type][1]) * gmx2lmp_sig
      f.write("%s %f %f\n" % (type, lmp_eps, lmp_sig))
   except:
      f.write("%s\n" % (type))
if n_bonds:
   f.write('\n')
   f.write('Bond Coeffs\n')
   f.write('\n')
   for type in bond_coeffs:
      try:
         k = float(bond_coeffs[type][0]) * gmx2lmp_b_k
         len = float(bond_coeffs[type][1]) * gmx2lmp_b_len
         f.write("%d %f %f\n" % (type, k, len))
      except:
         f.write("%d %s\n" % (type, bond_coeffs[type]))
if n_angles:
   f.write('\n')
   f.write('Angle Coeffs\n')
   f.write('\n')
   for type in angle_coeffs:
      try:
         k = float(angle_coeffs[type][0]) * gmx2lmp_a_k
         ang = float(angle_coeffs[type][1]) * gmx2lmp_a_ang
         f.write("%d %f %f\n" % (type, k, ang))
      except:
         f.write("%d %s\n" % (type, angle_coeffs[type]))

icoeff = [0.225, 0.275, 1.5, -0.4, -1.6, 0.0]

def dih2RB(c):
   return [-1.5*c[3]-2*c[1], -c[4]-c[2], -0.5*c[3], -0.25*c[4]]

if n_dihedrals:
   f.write('\n')
   f.write('Dihedral Coeffs\n')
   f.write('\n')
   for type in dihedral_coeffs:
      try:
         C = dihedral_coeffs[type]
         C = [float(i) for i in C]
         V = dih2RB(C)
         V = [i * gmx2lmp_d_k for i in V]
         f.write("%d %f %f %f %f\n" % (type, V[0], V[1], V[2], V[3]))
      except:
         f.write("%d %s\n" % (type, dihedral_coeffs[type]))
if n_impropers:
   f.write('\n')
   f.write('Improper Coeffs\n')
   f.write('\n')
   for type in improper_coeffs:
      try:
         C = improper_coeffs[type]
         C = [float(i) for i in C]
         V = dih2RB(C)
         V = [i * gmx2lmp_d_k for i in V]
         f.write("%d %f %f %f %f\n" % (type, V[0], V[1], V[2], V[3]))
      except:
         f.write("%d %s\n" % (type, dihedral_coeffs[type]))
f.close()

print("Performing a charge summation test..")
default_charge = opls_charge = 0.0
for atom in atom_list.values():
   default_charge += float(atom.default_charge)
   opls_charge += float(atom.opls_charge)
print("   default = " + str(default_charge))
print("   opls = " + str(opls_charge))
