from ase import io
from ase import __version__
from numpy import linalg as LA

print("ase version: {:s}".format(__version__))

#
# Testing ASE format
#

for [inxyz, outxyz] in [["../i.bi2se3_ase.xyz", "o.bi2se3_ase.xyz"], \
                        ["../i.bi2se3_nep.xyz", "o.bi2se3_nep.xyz"]]:

   ## --- Load structure ---
   print()
   print('*************************************')
   print('  Loading structure {:s}'.format(inxyz))
   print('*************************************')

   structure = io.read(inxyz)

   print('\n  Lattice vectors:')
   try:
      cell = structure.get_cell()
      print(cell[0])
      print(cell[1])
      print(cell[2])

      lx = LA.linalg.norm(cell[0])
      ly = LA.linalg.norm(cell[1])
      lz = LA.linalg.norm(cell[2])
      if lx == 0.0 or ly == 0.0 or lz == 0.0:
          print(  "  Warning: a cell dimension is zero")
   except:
      print("  Warning: could not parse lattice vectors..")

   print('\n  positions:')
   try:
      print(structure.get_positions())
   except:
      print("  Warning: could not parse positions..")

   print('\n  forces:')
   try:
      print(structure.get_forces())       # Nx3 forces
   except:
      print("  Warning: could not parse forces..")

   print('\n  Energy:')
   try:
      print(structure.get_potential_energy())
   except:
      print("  Warning: could not parse energy..")

   print('\n  chemical_symbols:')
   try:
      print(structure.get_chemical_symbols())  # ['Mo', 'Se', 'Se', 'Mo']
   except:
      print("  Warning: could not parse chemical_symbols..")

   print('\n  info:')
   try:
      print(structure.info)
   except:
      print("  Warning: could not parse info")

   #print("All attributes and methods:", structure.__dir__())

   print('\n  Exporting structure to {:s}'.format(outxyz))
   io.write(outxyz, structure)
