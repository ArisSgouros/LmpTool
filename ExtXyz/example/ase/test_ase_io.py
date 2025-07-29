from ase import io
from ase import __version__

print("ase version: {:s}".format(__version__))

#
# Testing ASE format
#

inxyz = "../i.bi2se3_ase.xyz"
outxyz = "o.bi2se3_ase.xyz"

## --- Load structure ---
print('\nLoading structure {:s}'.format(inxyz))

structure = io.read(inxyz)

print('\n  Lattice vectors:')
try:
   print(structure.get_cell())         # 3x3 lattice vectors
except:
   print("Warning: could not parse lattice vectors..")

print('\n  positions:')
try:
   print(structure.get_positions())
except:
   print("Warning: could not parse positions..")

print('\n  forces:')
try:
   print(structure.get_forces())       # Nx3 forces
except:
   print("Warning: could not parse forces..")

print('\n  Energy:')
try:
   print(structure.get_potential_energy())
except:
   print("Warning: could not parse energy..")

print('\n  chemical_symbols:')
try:
   print(structure.get_chemical_symbols())  # ['Mo', 'Se', 'Se', 'Mo']
except:
   print("Warning: could not parse chemical_symbols..")

print('\n  info:')
try:
   print(structure.info)
except:
   print("Warning: could not parse info")

#print("All attributes and methods:", structure.__dir__())

print('\nExporting structure to {:s}'.format(outxyz))
io.write(outxyz, structure)



#
# Testing NEP format
#

inxyz = "../i.bi2se3_nep.xyz"
outxyz = "o.bi2se3_nep.xyz"

## --- Load structure ---
print('\nLoading structure {:s}'.format(inxyz))

structure = io.read(inxyz)

print('\n  Lattice vectors:')
try:
   print(structure.get_cell())         # 3x3 lattice vectors
except:
   print("Warning: could not parse lattice vectors..")

print('\n  positions:')
try:
   print(structure.get_positions())
except:
   print("Warning: could not parse positions..")

print('\n  forces:')
try:
   print(structure.get_forces())       # Nx3 forces
except:
   print("Warning: could not parse forces..")

print('\n  Energy:')
try:
   print(structure.get_potential_energy())
except:
   print("Warning: could not parse energy..")

print('\n  chemical_symbols:')
try:
   print(structure.get_chemical_symbols())  # ['Mo', 'Se', 'Se', 'Mo']
except:
   print("Warning: could not parse chemical_symbols..")

print('\n  info:')
try:
   print(structure.info)
except:
   print("Warning: could not parse info")

#print("All attributes and methods:", structure.__dir__())

print('\nExporting structure to {:s}'.format(outxyz))
io.write(outxyz, structure)
