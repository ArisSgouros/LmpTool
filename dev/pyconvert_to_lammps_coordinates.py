import numpy as np

def transform_to_upper_triangular_lattice(a, b, c, atoms):
    """
    Transforms lattice vectors and atomic coordinates to an upper triangular form:
    a' = (a'x, a'y, a'z)
    b' = (0, b'y, b'z)
    c' = (0, 0, c'z)

    Parameters:
        a, b, c: lists or arrays of shape (3,) representing lattice vectors
        atoms: array of shape (N, 3) of atomic coordinates in Cartesian units

    Returns:
        a_new, b_new, c_new: new lattice vectors
        atoms_new: transformed atomic coordinates
    """
    a = np.array(a)
    b = np.array(b)
    c = np.array(c)
    atoms = np.array(atoms)

    # Step 1: Normalize a to get e1 (x-axis)
    e1 = a / np.linalg.norm(a)

    # Step 2: Remove projection of b along a to get e2 (in b–a plane, orthogonal to a)
    b_proj = b - np.dot(b, e1) * e1
    e2 = b_proj / np.linalg.norm(b_proj)

    # Step 3: Make e3 orthogonal to both e1 and e2
    e3 = np.cross(e1, e2)

    # Compose transformation matrix: columns are the new basis vectors
    R = np.vstack([e1, e2, e3])

    # Transform lattice vectors and atoms to new coordinate system
    a_new = R @ a
    b_new = R @ b
    c_new = R @ c
    atoms_new = (R @ atoms.T).T

    return a_new, b_new, c_new, atoms_new

# Example usage
if __name__ == "__main__":
    # Define lattice vectors
    a = [1.0, 0.0, 0.0]
    b = [0.5, 1.0, 0.0]
    c = [0.2, 0.3, 1.0]

    # Atomic coordinates (Cartesian)
    atoms = [
        [0.1, 0.1, 0.1],
        [0.5, 0.5, 0.5]
    ]

    a_new, b_new, c_new, atoms_new = transform_to_upper_triangular_lattice(a, b, c, atoms)

    print("Transformed Lattice Vectors:")
    print("a' =", a_new)
    print("b' =", b_new)
    print("c' =", c_new)

    print("\nTransformed Atomic Coordinates:")
    for atom in atoms_new:
        print(atom)

