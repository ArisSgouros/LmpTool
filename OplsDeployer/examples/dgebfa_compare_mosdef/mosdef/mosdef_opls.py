import sys
import parmed as pmd
import mbuild as mb

from foyer import Forcefield
from foyer.tests.utils import get_fn

if __name__ == "__main__":
    name = sys.argv[1]
    ff_path = sys.argv[2]
    mol2_path = name + ".mol2"

    untyped_group = pmd.load_file(mol2_path, structure=True)
    oplsaa_aps = Forcefield(forcefield_files=ff_path)
    group = oplsaa_aps.apply(untyped_group, references_file=name+".bib", verbose=True)

    verbose = True
    if verbose:
        print("Atoms:")
        for atom in group.atoms:
            print("Atom {} is typed as {}".format(atom, atom.type))

        print("Bonds:")
        for bond in group.bonds:
            print("{} ".format(bond))

        print("Angles:")
        for angle in group.angles:
            print("{} ".format(angle))

        print("Dihedrals:")
        for dihedral in group.dihedrals:
            print("{} ".format(dihedral))

    # Save to Lammps [https://mbuild.mosdef.org/en/stable/getting_started/writers/LAMMPS_file_writers.html]
    # bug https://github.com/ParmEd/ParmEd/issues/1029
    boxlo = [0.0, 0.0, 0.0]
    boxhi = [0.1*group.box[0], 0.1*group.box[1], 0.1*group.box[2]]
    mb.formats.lammpsdata.write_lammpsdata(group, name+"_mosdef.data", atom_style='full', mins=boxlo, maxs=boxhi)

    # Save to GROMACS
    #group.save(name+".gro")
    #group.save(name+".top")

    # Within the `Forcefield.apply` method, an intermediate OpenMM system is
    # created. If you wish to use OpenMM, e.g. for use of potential forms not
    # yet supported by ParmEd, you can simply stop the conversion process
    # after the OpenMM System creation by directly invoking the internal
    # method calls:
    #from foyer.forcefield import generate_topology
    #
    #omm_topology, positions = generate_topology(untyped_group)
    #omm_system = oplsaa_aps.createSystem(topology=omm_topology)
