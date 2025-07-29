import sys
import re
import json

kRyToEv = 13.6057039763
kBohrToAngstrom = 0.529177
#debug = False

def read_json(file):
   with open(file, 'r') as f:
      pp = json.load(f)
      return pp

def write_json(file, pp):
   with open(file, 'w') as f:
      json.dump(pp, f)

def parse_xyz_frame(frame_lines):
    """
    Parses a single XYZ-like frame from a list of strings.

    Args:
        frame_lines (list): A list of strings representing a single frame.

    Returns:
        dict: A dictionary containing the parsed data for the frame.
    """
    if not frame_lines or len(frame_lines) < 2:
        return None

    num_atoms = int(frame_lines[0].strip())
    properties_line = frame_lines[1].strip()

    # Parse cell properties
    energy_match = re.search(r'energy=(-?\d+\.?\d*)', properties_line)
    energy = float(energy_match.group(1)) if energy_match else None

    config_type_match = re.search(r'config_type=([a-zA-Z0-9._]+)', properties_line)
    config_type = str(config_type_match.group(1)) if config_type_match else None

    weight_match = re.search(r'weight=(-?\d+\.?\d*)', properties_line)
    weight = float(weight_match.group(1)) if weight_match else None

    lattice_match = re.search(r'lattice="([^"]*)"', properties_line, re.IGNORECASE)
    lattice_str = lattice_match.group(1) if lattice_match else ""
    lattice_vectors = [float(x) for x in lattice_str.split()]

    virial_match = re.search(r'virial="([^"]*)"', properties_line)
    virial_str = virial_match.group(1) if virial_match else ""
    virial = [float(x) for x in virial_str.split()]

    # Corrected parsing for properties
    format_match = re.search(r'Properties=([a-zA-Z0-9:]+)', properties_line, re.IGNORECASE)
    format_str = format_match.group(1) if format_match else ""
    atom_format = []
    format_parts = format_str.split(':')
    for i in range(0, len(format_parts), 3):
        if i + 2 < len(format_parts):
            atom_format.append({
                'name': format_parts[i],
                'type': format_parts[i+1],
                'size': int(format_parts[i+2])
            })

    atoms_data = []
    for line in frame_lines[2:]:
        parts = line.split()
        if not parts:
            continue

        atom_info = {}
        idx = 0
        for prop_def in atom_format:
            prop_name = prop_def['name']
            prop_type = prop_def['type']
            prop_size = prop_def['size']

            if prop_type == 'S':
                atom_info[prop_name] = parts[idx]
                idx += 1
            elif prop_type == 'R':
                if idx + prop_size <= len(parts):
                    atom_info[prop_name] = [float(x) for x in parts[idx:idx + prop_size]]
                    idx += prop_size
                else:
                    print(f"Warning: Not enough data for {prop_name} in line: {line.strip()}")
                    break
        atoms_data.append(atom_info)

    return {
        'num_atoms': num_atoms,
        'cell_properties': {
            'energy': energy,
            'lattice_vectors': lattice_vectors,
            'virial': virial,
            'atom_section_format': atom_format,
            'config_type': config_type,
            'weight': weight,
        },
        'atoms': atoms_data
    }

def read_xyz_file(file_path):
    """
    Reads an XYZ-like file containing multiple frames.

    Args:
        file_path (str): The path to the XYZ-like file.

    Returns:
        list: A list of dictionaries, where each dictionary represents a parsed frame.
    """
    all_frames_data = []
    with open(file_path, 'r') as f:
        lines = f.readlines()

    current_frame_lines = []
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if not line:
            i += 1
            continue

        # Check if the line indicates the start of a new frame (number of atoms)
        if line.isdigit():
            if current_frame_lines: # If there was a previous frame being collected
                parsed_frame = parse_xyz_frame(current_frame_lines)
                if parsed_frame:
                    all_frames_data.append(parsed_frame)
                current_frame_lines = [] # Reset for the new frame

            num_atoms = int(line)
          
            current_frame_lines.append(line)
            # Add the next (num_atoms + 1) lines to current_frame_lines
            # +1 for the properties line
            for j in range(1, num_atoms + 2):
                if i + j < len(lines):
                    current_frame_lines.append(lines[i+j].strip())
            i += num_atoms + 2 # Move index past the current frame
        else:
            i += 1 # Move to the next line if it's not a start of a frame

    # Add the last frame if any lines were collected
    if current_frame_lines:
        parsed_frame = parse_xyz_frame(current_frame_lines)
        if parsed_frame:
            all_frames_data.append(parsed_frame)

    return all_frames_data

def export_xyz_frame(frame_data, xyzfmt = 'ase'):

    """
    Exports a single XYZ-like frame from a dictionary to a list of strings.

    Args:
        frame_data (dict): A dictionary containing the frame data,
                           matching the format returned by parse_xyz_frame.
        xyzfmt (str): select format of extended xyz files

    Returns:
        list: A list of strings representing the XYZ frame, or None if invalid data.
    """

    if xyzfmt == 'ase':
       kLatticeFmt = "Lattice"
       kPropertyFmt = "Properties"
       kForceFmt = "forces"
    elif xyzfmt == 'nep':
       kLatticeFmt = "lattice"
       kPropertyFmt = "properties"
       kForceFmt = "force"
    else:
        print(f"Error: Unknown format {xyzfmt}")
        return None

    if not isinstance(frame_data, dict):
        print("Error: Input frame_data must be a dictionary.")
        return None

    num_atoms = frame_data.get('num_atoms')
    cell_properties = frame_data.get('cell_properties', {})
    atoms_data = frame_data.get('atoms', [])

    if num_atoms is None or not isinstance(num_atoms, int) or num_atoms < 0:
        print("Error: 'num_atoms' missing or invalid in frame_data.")
        return None
    if not isinstance(atoms_data, list) or len(atoms_data) != num_atoms:
        print(f"Error: Number of atoms in 'atoms' list ({len(atoms_data)}) does not match 'num_atoms' ({num_atoms}).")
        return None

    # --- Construct the properties line ---
    properties_parts = []

    # Energy
    energy = cell_properties.get('energy')
    if energy is not None:
        properties_parts.append(f"energy={energy:.9f}") # Format to 9 decimal places

    # lattice
    lattice_vectors = cell_properties.get('lattice_vectors')
    if lattice_vectors and isinstance(lattice_vectors, list) and all(isinstance(x, (int, float)) for x in lattice_vectors):
        lattice_str = " ".join([f"{x:.9f}" for x in lattice_vectors])
        properties_parts.append(f'{kLatticeFmt}="{lattice_str}"')

    # virial
    virial = cell_properties.get('virial')
    if virial and isinstance(virial, list) and all(isinstance(x, (int, float)) for x in virial):
        virial_str = " ".join([f"{x:.9f}" for x in virial])
        properties_parts.append(f'virial="{virial_str}"')

    # Atom Section Format
    atom_format = cell_properties.get('atom_section_format')
    if atom_format and isinstance(atom_format, list):
        format_str_parts = []
        for prop_def in atom_format:
            name = prop_def.get('name')
            prop_type = prop_def.get('type')
            size = prop_def.get('size')
            if name in ['force', 'forces']:
                name = kForceFmt
            if name and prop_type and isinstance(size, int):
                format_str_parts.append(f"{name}:{prop_type}:{size}")
        if format_str_parts:
            properties_parts.append(f'{kPropertyFmt}=' + ':'.join(format_str_parts))

    # Config type
    config_type = cell_properties.get('config_type')
    if config_type is not None:
        properties_parts.append(f"config_type={config_type:s}") # Format to 9 decimal places

    # Weight
    weight = cell_properties.get('weight')
    if weight is not None:
        properties_parts.append(f"weight={weight:.1f}") # Format to 9 decimal places

    properties_line = " ".join(properties_parts)

    # --- Construct atom data lines ---
    atom_lines = []
    for atom_info in atoms_data:
        line_parts = []
        for prop_def in atom_format: # Use the same format to ensure order
            prop_name = prop_def['name']
            prop_type = prop_def['type']
            prop_size = prop_def['size']
            value = atom_info.get(prop_name)

            if value is None:
                # Handle missing properties gracefully, maybe skip or add placeholder
                # For now, print a warning and return None for malformed data
                print(f"Warning: Missing property '{prop_name}' for an atom. Cannot export frame.")
                return None

            if prop_type == 'S':
                line_parts.append(str(value))
            elif prop_type == 'R':
                if isinstance(value, list) and len(value) == prop_size:
                    line_parts.extend([f"{x:.9f}" for x in value]) # Format floats for consistency
                else:
                    print(f"Error: Invalid data for '{prop_name}' (expected list of size {prop_size}): {value}")
                    return None
            # Add other types if your format supports them (e.g., 'I' for integer)
            else:
                print(f"Warning: Unknown property type '{prop_type}' for '{prop_name}'. Skipping.")

        atom_lines.append(" ".join(line_parts))

    # Combine all parts
    output_lines = [
        str(num_atoms),
        properties_line
    ]
    output_lines.extend(atom_lines)

    return output_lines


def write_xyz_file(file_path, frames_data, status, xyzfmt):
    """
    Writes a list of parsed XYZ frame dictionaries to an XYZ-like file.

    Args:
        file_path (str): The path to the output XYZ file.
        frames_data (list): A list of dictionaries, where each dictionary
                            represents a frame to be written.
    """
    if not isinstance(frames_data, list):
        print("Error: Input frames_data must be a list of dictionaries.")
        return

    try:
        with open(file_path, status) as f:
            for i, frame_data in enumerate(frames_data):
                frame_lines = export_xyz_frame(frame_data, xyzfmt)
                if frame_lines:
                    f.write("\n".join(frame_lines))
                    f.write("\n") # Add a newline between frames
                else:
                    print(f"Warning: Skipping frame {i} due to invalid data.")
        print(f"Successfully exported data to '{file_path}'.")
    except Exception as e:
        print(f"Error writing XYZ file '{file_path}': {e}")


def replace_string_in_file(file_path, old_string, new_string):
    """
    Replaces all occurrences of old_string with new_string in a file.

    Args:
        file_path (str): The path to the file.
        old_string (str): The string to be replaced.
        new_string (str, int, float): The string to replace with.
    """
    try:
        # Make sure that the input value is a string
        new_string = str(new_string)

        # Read the entire file content into memory
        with open(file_path, 'r') as file:
            file_content = file.read()

        # Perform the replacement
        new_content = file_content.replace(old_string, new_string)

        # Write the modified content back to the file
        with open(file_path, 'w') as file:
            file.write(new_content)

        print(f"Successfully replaced '{old_string}' with '{new_string}' in '{file_path}'.")

    except FileNotFoundError:
        print(f"Error: File not found at '{file_path}'.")
    except Exception as e:
        print(f"An error occurred: {e}")

def append_line_to_file(file_path, string_to_append):
    """
    Appends a string to the end of a file.

    Args:
        file_path (str): The path to the file.
        string_to_append (str): The string to be appended.
    """
    try:
        # Open the file in append mode ('a')
        # The 'with' statement ensures the file is automatically closed
        with open(file_path, 'a') as file:
            file.write("%s\n" % string_to_append)
            print(f"Successfully appended '{string_to_append}' to '{file_path}'.")
    except FileNotFoundError:
        print(f"Error: File not found at '{file_path}'.")
    except Exception as e:
        print(f"An error occurred: {e}")

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
        virial_ev = [sij*unit_cell_volume_bohr3*kRyToEv for sij in stress_tensor_Ry_bohr3]

        parsed_frame['cell_properties']['virial'] = virial_ev

    # Parse ATOMIC_POSITIONS
    atoms_data = []
    if line_atomic_positions:
        for ii in range(100000):
           line = lines[line_atomic_positions + ii]
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

