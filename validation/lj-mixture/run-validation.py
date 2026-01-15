# /// script
# requires-python = ">=3.13"
# dependencies = [
#     "matplotlib>=3.10.8",
#     "numpy>=2.4.1",
#     "pandas>=2.3.3",
#     "seaborn>=0.13.2",
# ]
# ///

# Comes from Monte Carlo Simulations of Binary Lennard–Jones Mixtures: A Test of the van der Waals One-Fluid Model, https://doi.org/10.1023/A:1022614200488

import itertools
import json
import math
import os
import subprocess

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# Different from publication, because we start from rectangular lattice instead of fcc
# It shouldn't affect results, as we're not looking at crystallisation
N = 1000
path_to_exe = "particlesmc"

epsilon_1 = 1
epsilon_12 = 1.1523
epsilon_2 = 1.3702
sigma_1 = 1
sigma_12 = 1.0339
sigma_2 = 1.0640


def create_config(n1: int, n2: int, box_length: float, filepath: str) -> None:
    # We do a square lattice
    # Species are put randomly amongst the sites

    with open(filepath, "w") as f:
        f.write(f"{N:d}\n")
        f.write(
            f'Lattice="{box_length:.4f} 0.0 0.0 0.0 {box_length:.4f} 0.0 0.0 0.0 {box_length:.4f}" Properties=:species:S:1:pos:R:3\n'
        )

        species_indices = [1] * n1 + [2] * n2
        np.random.shuffle(species_indices)

        number_of_particle_in_each_direction = round(N ** (1 / 3))
        if number_of_particle_in_each_direction**3 != N:
            raise RuntimeError("N is not x^3")
        dxdydz = box_length / number_of_particle_in_each_direction
        counter = 0
        for i, j, k in itertools.product(
            range(number_of_particle_in_each_direction),
            range(number_of_particle_in_each_direction),
            range(number_of_particle_in_each_direction),
        ):
            x = (i + 0.5) * dxdydz - box_length / 2
            y = (j + 0.5) * dxdydz - box_length / 2
            z = (k + 0.5) * dxdydz - box_length / 2
            f.write(f"{species_indices[counter]} {x} {y} {z}\n")
            counter += 1


def create_params(parameters: dict, path_to_params: str) -> None:
    with open(path_to_params, "w") as f:
        f.write(f"""
        [system]
        config = "{parameters["config"]}"
        temperature = {parameters["temperature"]}
        density = {parameters["density"]}
        list_type = "LinkedList"

        [model]
        [model."1-1"]
        name = "LennardJones"
        epsilon = {parameters["epsilon_1"]:.2f}
        sigma = {parameters["sigma_1"]:.2f}
        rcut = {parameters["rcut"]:.2f}
        shift_potential = false

        [model."1-2"]
        name = "LennardJones"
        epsilon = {parameters["epsilon_12"]}
        sigma = {parameters["sigma_12"]}
        rcut = {parameters["rcut"]:.2f}
        shift_potential = false

        [model."2-2"]
        name = "LennardJones"
        epsilon = {parameters["epsilon_2"]}
        sigma = {parameters["sigma_2"]}
        rcut = {parameters["rcut"]:.2f}
        shift_potential = false

        [simulation]
        type = "Metropolis"
        steps = {int(parameters["steps"]):d}
        seed = 42
        parallel = false
        output_path = "./"

        [[simulation.move]]
        action = "Displacement"
        probability = 1.0
        policy = "SimpleGaussian"
        parameters = {{sigma = 0.05}}

        [[simulation.output]]
        algorithm = "StoreCallbacks"
        callbacks = ["energy", "acceptance"]
        scheduler_params = {{linear_interval = 100}}

        [[simulation.output]]
        algorithm = "StoreLastFrames"
        scheduler_params = {{linear_interval = 1000}}
        fmt = "EXYZ"

        """)


def run_simulations(output_path: str) -> None:
    df = pd.read_csv("reference-data.csv")

    path_to_config = "config.exyz"
    path_to_params = "params.toml"
    data = []
    for i, row in df.iterrows():
        workdir = f"./tmp/{i}"
        os.makedirs(workdir, exist_ok=True)

        print(row["t"], row["x"], row["density"])

        # density = N / box_length**3
        box_length = (N / row["density"]) ** (1 / 3)

        n2 = round(N * row["x"])
        n1 = N - n2
        create_config(n1, n2, box_length, f"{workdir}/{path_to_config}")

        # Generate input, and run
        rc = 4 * sigma_1
        parameters = {
            "config": path_to_config,
            "temperature": row["t"],
            "epsilon_1": epsilon_1,
            "epsilon_12": epsilon_12,
            "epsilon_2": epsilon_2,
            "sigma_1": sigma_1,
            "sigma_12": sigma_12,
            "sigma_2": sigma_2,
            "steps": 1e3,  # 1e4 equilibration in paper
            "rcut": rc,
            "density": N / box_length**3,
        }
        create_params(parameters, f"{workdir}/{path_to_params}")

        subprocess.run(
            [path_to_exe, path_to_params], cwd=workdir, stdout=subprocess.DEVNULL
        )

        # Post-process the energies
        energies = pd.read_csv(f"{workdir}/energy.dat", sep="\\s+", names=["i", "e"])[
            "e"
        ]
        # Remove the first half as equilibration, just to be sure
        energies = energies[int(len(energies) / 2) :]

        acceptance_rate = np.array(
            pd.read_csv(f"{workdir}/acceptance.dat", sep="\\s+", names=["i", "a"])["a"]
        )[-1]

        # Compute long-range corrections from the cutoff
        # Formula from Gromacs https://manual.gromacs.org/current/reference-manual/functions/long-range-vdw.html
        lr_correction = 0

        c6_11 = 4 * epsilon_1 * sigma_1**6
        rho_11 = n1 / box_length**3
        lr_correction += -2 / 3 * math.pi * n1 * rho_11 * c6_11 / rc**3

        c6_22 = 4 * epsilon_2 * sigma_2**6
        rho_22 = n2 / box_length**3
        lr_correction += -2 / 3 * math.pi * n2 * rho_22 * c6_22 / rc**3

        c6_12 = 4 * epsilon_12 * sigma_12**6
        lr_correction += -2 / 3 * math.pi * (n1 * rho_22 + n2 * rho_11) * c6_12 / rc**3

        data.append(
            {
                "t": row["t"],
                "x": row["x"],
                "density": row["density"],
                "energy": np.mean(energies),
                "energy_err": np.std(energies) / np.sqrt(len(energies)),
                "lr_correction": lr_correction,
                "acceptance_rate": json.loads(acceptance_rate)[0],
            }
        )

        df.merge(pd.DataFrame(data)).to_csv(output_path)


def plot_results(path_to_energies: str) -> None:
    df = pd.read_csv(path_to_energies)

    # u (from the paper) is actually total energy / epsilon_11 N
    # epsilon_11 = 1, so it's the energy per particule
    # Add long range corrections to the MC results, which is also in energy / particle
    df["energy_lr"] = df["energy"] + df["lr_correction"] / N

    sns.scatterplot(data=df, x="u", y="energy_lr", hue="density", style="x")
    plt.axline((-5, -5), slope=1)
    plt.xlabel("Published energy")
    plt.ylabel("Re-computed from ParticleMC")
    plt.savefig("correlation-plot.jpeg")
    plt.show()


if __name__ == "__main__":
    path_to_energies = "calculated-energies.csv"
    run_simulations(path_to_energies)
    plot_results(path_to_energies)
