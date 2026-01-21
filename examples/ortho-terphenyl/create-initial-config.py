# /// script
# requires-python = ">=3.13"
# dependencies = [
#     "numpy>=2.4.1",
# ]
# ///
import itertools

import numpy as np

M = 1000
N = M * 3
number_density = 1.2  # This is per particle, not per molecule
# rho = N/V
V = N / number_density
L = V ** (1 / 3)

sigmas = [0.9, 1, 1.1]


def main(output_file_name: str) -> None:
    """Create a regular lattice of molecules. These planar molecules are in the x-y plane."""

    number_in_each_direction = round(M ** (1 / 3))
    dxdydz = L / number_in_each_direction
    assert number_in_each_direction**3 == M, "M is not an int to power 3"

    r_ab = (sigmas[0] + sigmas[1]) / 2
    r_ac = (sigmas[0] + sigmas[2]) / 2
    cos_alpha = np.cos(60 / 180)
    sin_alpha = np.sin(60 / 180)

    f = open(output_file_name, "w")

    f.write(f"{N}\n")
    f.write(
        f"columns:molecule,species,position cell:{L},{L},{L},rho:{number_density}\n"
    )

    counter = 1
    for i, j, k in itertools.product(
        range(number_in_each_direction),
        range(number_in_each_direction),
        range(number_in_each_direction),
    ):
        r_a = (i * dxdydz, j * dxdydz, k * dxdydz)
        r_b = (i * dxdydz, j * dxdydz + r_ab, k * dxdydz)
        r_c = (i * dxdydz + r_ac * cos_alpha, j * dxdydz + r_ac * sin_alpha, k * dxdydz)

        f.write(f"{counter} 1 {r_a[0]} {r_a[1]} {r_a[2]}\n")
        f.write(f"{counter} 2 {r_b[0]} {r_b[1]} {r_b[2]}\n")
        f.write(f"{counter} 3 {r_c[0]} {r_c[1]} {r_c[2]}\n")

        counter += 1

    number_of_bonds = M * 3
    f.write(f"{number_of_bonds}\n")
    f.write("columns:bond\n")
    for i in range(M):
        index_of_a = 1 + i * 3
        # AB
        f.write(f"{index_of_a} {index_of_a + 1}\n")
        # BC
        f.write(f"{index_of_a + 1} {index_of_a + 2}\n")
        # AC
        f.write(f"{index_of_a} {index_of_a + 2}\n")


if __name__ == "__main__":
    output_file_name = "inputframe.xyz"
    main(output_file_name)
