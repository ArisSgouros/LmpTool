###############################################################################
# MIT License
#
# Copyright (c) 2025 ArisSgouros
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
import re
import json

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

