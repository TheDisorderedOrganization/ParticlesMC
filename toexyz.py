import numpy as np

def lammps_to_extended_xyz(lammps_data):
    """Converts LAMMPS data to extended XYZ format.

    Args:
        lammps_data: A string containing the LAMMPS data.

    Returns:
        A string containing the extended XYZ data, or None if an error occurs.
    """

    try:
        lines = lammps_data.strip().split('\n')

        num_atoms_line = next(i for i, line in enumerate(lines) if "NUMBER OF ATOMS" in line)
        num_atoms = int(lines[num_atoms_line + 1].strip())

        box_bounds_line = next(i for i, line in enumerate(lines) if "BOX BOUNDS" in line)
        box_bounds = []
        for i in range(1, 4):
            bounds = list(map(float, lines[box_bounds_line + i].split()))
            box_bounds.append(bounds)
        
        x_lo, x_hi = box_bounds[0]
        y_lo, y_hi = box_bounds[1]
        z_lo, z_hi = box_bounds[2]
        lx = x_hi - x_lo
        ly = y_hi - y_lo
        lz = z_hi - z_lo

        atoms_start_line = next(i for i, line in enumerate(lines) if "ATOMS type" in line) + 1
        atoms_end_line = atoms_start_line + num_atoms

        extended_xyz_data = f"{num_atoms}\n"
        extended_xyz_data += f"Lattice=\"{lx} 0.0 0.0 0.0 {ly} 0.0 0.0 0.0 {lz}\" Properties=species:S:1:pos:R:3\n"

        for line in lines[atoms_start_line:atoms_end_line]:

            parts = line.split()
            atom_type = int(parts[0])
            x = float(parts[1])
            y = float(parts[2])
            extended_xyz_data += f"{atom_type} {x} {y}\n"

        return extended_xyz_data

    except (ValueError, StopIteration, IndexError) as e:
        print(f"Error converting LAMMPS data: {e}")
        return None

with open("test/config_0.lmp", "r") as g:
    extended_xyz = lammps_to_extended_xyz(g.read())
    if extended_xyz:
        print(extended_xyz)
        # You can save this to a file:
        with open("config_0.exyz", "w") as f:
            f.write(extended_xyz)
