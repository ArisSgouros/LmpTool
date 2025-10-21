import sys
import re
import json

kRyToEv = 13.6057039763
kRyToJoule = 2.1798741E-18
kBohrToAngstrom = 0.529177
kEvToJoule = 1.60218e-19
kPaToAtm = 1.0 / 101325
kAngstromToMeter = 1e-10
kRy_Bohr3ToPa = kRyToJoule/(kBohrToAngstrom*kAngstromToMeter)**3
kRy_Bohr3ToAtm = kRy_Bohr3ToPa*kPaToAtm

def parse_qe_output(file_path, xyzfmt = 'ase', debug = False):

    if xyzfmt == 'ase':
       kForceFmt = "forces"
    elif xyzfmt == 'nep':
       kForceFmt = "force"

    line_stress = None
    line_forces = None
    line_energy = None
    line_atomic_positions = None
    line_cartesian_axes = None
    line_cell_parameters = None
    line_crystal_axes = None
    line_unit_cell_volume = None
    line_alat = None

    is_job_done = False

    lines = []
    try:
        with open(file_path) as f:
            lines = [line.rstrip('\n') for line in f]

        ii = 0
        for line in lines:
            if "CELL_PARAMETERS" in line:
                line_cell_parameters = ii + 1
                if debug: print(ii, "CELL_PARAMETERS")
            if "crystal axes:" in line:
                line_crystal_axes = ii + 1
                if debug: print(ii, "crystal axes:")
            if "ATOMIC_POSITIONS" in line:
                line_atomic_positions = ii + 1
                if debug: print(ii, "ATOMIC_POSITIONS")
            if "Cartesian axes" in line:
                line_cartesian_axes = ii + 3
                if debug: print(ii, "Cartesian axes")
            if "total" in line and "stress" in line and "Ry/bohr**3" in line:
                line_stress = ii + 1
                if debug: print(ii, "total", "stress", "Ry/bohr**3")
            if "unit-cell volume" in line:
                line_unit_cell_volume = ii
                if debug: print(ii, "unit-cell volume")
            if "!    total energy" in line:
                line_energy = ii
                if debug: print(ii, "total energy")
            if "Forces acting on atoms" in line:
                line_forces = ii + 2
                if debug: print(ii, "Forces acting on atoms")
            if "lattice parameter (alat)" in line:
                line_alat = ii
                if debug: print(ii, "lattice parameter (alat)")
            if "JOB DONE" in line:
                is_job_done = True
            ii += 1

        if debug:
            if not line_energy:
                print(f"Warning: issue with parsing energy")
            if not line_atomic_positions and not line_cartesian_axes:
                print(f"Warning: issue with parsing atomic positions")
            if not line_stress:
                print(f"Warning: issue with parsing stress")
            if not line_unit_cell_volume and not line_crystal_axes:
                print(f"Warning: issue with parsing cell volume")
            if not line_cell_parameters:
                print(f"Warning: issue with parsing cell parameters")
            if not line_forces:
                print(f"Warning: issue with parsing forces")

    except FileNotFoundError:
        print(f"Error: File not found at '{file_path}'")
        return None
    except Exception as e:
        print(f"An error occurred while parsing the file: {e}")
        return None

    parsed_frame = {'num_atoms':None, 'cell_properties':{}, 'atoms': []}
    parsed_frame['cell_properties']['atom_section_format'] = []

    # Parse alat
    if line_alat:
        try:
            lsplit = lines[line_alat].split()
            alat_bohr = float(lsplit[4])
            print('alat (bohr): ', alat_bohr)
        except:
            print('Error with parsing the unit cell volume')

    # Parse energy
    if line_energy:
        try:
            lsplit = lines[line_energy].split()
            energy = float(lsplit[4])
            energy_to_ev = float(lsplit[4]) * kRyToEv
            parsed_frame['cell_properties']['energy'] = energy_to_ev
        except:
            print('Error with parsing the energy')

    # Parse cell_volume
    unit_cell_volume_bohr3 = None
    if line_unit_cell_volume:
        try:
            lsplit = lines[line_unit_cell_volume].split()
            offset = 3
            if 'new' in lsplit:
                offset = 4
            unit_cell_volume_bohr3 = float(lsplit[offset])
            print('unit_cell_volume: ',unit_cell_volume_bohr3)
        except:
            print('Error with parsing the unit cell volume')

    # Parse CELL_PARAMETERS
    if line_cell_parameters:
        lsplit = lines[line_cell_parameters + 0].split()
        ax = float(lsplit[0])
        ay = float(lsplit[1])
        az = float(lsplit[2])
        lsplit = lines[line_cell_parameters + 1].split()
        bx = float(lsplit[0])
        by = float(lsplit[1])
        bz = float(lsplit[2])
        lsplit = lines[line_cell_parameters + 2].split()
        cx = float(lsplit[0])
        cy = float(lsplit[1])
        cz = float(lsplit[2])
        lattice_vectors = [ax, ay, az, bx, by, bz, cx, cy, cz]
    elif line_crystal_axes:
        print('Warning: could not find CELL_PARAMETERS; will use init crystal axes')
        offset = 3
        lsplit = lines[line_crystal_axes + 0].split()
        ax = float(lsplit[offset + 0])
        ay = float(lsplit[offset + 1])
        az = float(lsplit[offset + 2])
        lsplit = lines[line_crystal_axes + 1].split()
        bx = float(lsplit[offset + 0])
        by = float(lsplit[offset + 1])
        bz = float(lsplit[offset + 2])
        lsplit = lines[line_crystal_axes + 2].split()
        cx = float(lsplit[offset + 0])
        cy = float(lsplit[offset + 1])
        cz = float(lsplit[offset + 2])
        lattice_vectors_bohr = [ax, ay, az, bx, by, bz, cx, cy, cz]
        lattice_vectors = [ij*alat_bohr*kBohrToAngstrom for ij in lattice_vectors_bohr]

    parsed_frame['cell_properties']['lattice_vectors'] = lattice_vectors


    # Parse stress
    if line_stress and line_unit_cell_volume:
        lsplit = lines[line_stress + 0].split()
        sxx = float(lsplit[0])
        sxy = float(lsplit[1])
        sxz = float(lsplit[2])
        lsplit = lines[line_stress + 1].split()
        syx = float(lsplit[0])
        syy = float(lsplit[1])
        syz = float(lsplit[2])
        lsplit = lines[line_stress + 2].split()
        szx = float(lsplit[0])
        szy = float(lsplit[1])
        szz = float(lsplit[2])

        stress_tensor_Ry_bohr3 = [sxx, sxy, sxz, syx, syy, syz, szx, szy, szz]
        virial_tensor_ev = [sij*unit_cell_volume_bohr3*kRyToEv for sij in stress_tensor_Ry_bohr3]
        stress_tensor_atm = [sij*kRy_Bohr3ToAtm for sij in stress_tensor_Ry_bohr3]

        parsed_frame['cell_properties']['virial'] = virial_tensor_ev
        parsed_frame['cell_properties']['stress'] = stress_tensor_atm

    # Parse ATOMIC_POSITIONS
    atoms_data = []
    if line_atomic_positions:
        for ii in range(100000):
           line = lines[line_atomic_positions + ii]
           if not line.strip(): # Check if empty line
              print("Warning: these are not final coordinates. Ignore if this is input file")
              break
           if "End final coordinates" in line:
              break

           lsplit = line.split()

           species = lsplit[0]
           x = float(lsplit[1])
           y = float(lsplit[2])
           z = float(lsplit[3])

           atom_info = {'species': species, 'pos': [x, y, z]}
           parsed_frame['atoms'].append(atom_info)

    elif line_cartesian_axes:
        for ii in range(100000):
           line = lines[line_cartesian_axes + ii]

           lsplit = line.split()
           if not lsplit:
              break

           species = lsplit[1]
           offset = 6
           x = float(lsplit[offset + 0])
           y = float(lsplit[offset + 1])
           z = float(lsplit[offset + 2])

           pos_alat = [x, y, z]
           pos_angstrom = [aa*alat_bohr*kBohrToAngstrom for aa in pos_alat]

           atom_info = {'species': species, 'pos': pos_angstrom}
           parsed_frame['atoms'].append(atom_info)

    # Specify the species/pos format
    fmt = {'name': 'species', 'type': 'S', 'size': 1}
    parsed_frame['cell_properties']['atom_section_format'].append(fmt)

    fmt = {'name': 'pos', 'type': 'R', 'size': 3}
    parsed_frame['cell_properties']['atom_section_format'].append(fmt)

    # Number of atoms
    nat = len(parsed_frame['atoms'])
    parsed_frame['num_atoms'] = nat

    # Parse forces
    if line_forces and (line_atomic_positions or line_cartesian_axes):
        for ii in range(nat):
           line = lines[line_forces + ii]
           lsplit = line.split()

           fx = float(lsplit[6])
           fy = float(lsplit[7])
           fz = float(lsplit[8])

           ff_Ry_au = [fx, fy, fz]
           ff_eV_angstrom = [fi*kRyToEv/kBohrToAngstrom for fi in ff_Ry_au]

           parsed_frame['atoms'][ii][kForceFmt] = ff_eV_angstrom

        fmt = {'name': kForceFmt, 'type': 'R', 'size': 3}
        parsed_frame['cell_properties']['atom_section_format'].append(fmt)

    status = {}
    status['job_done'] = is_job_done

    return parsed_frame, status

