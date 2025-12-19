# Running simulations

To perform a MC simulation of atomic or molecular systems, the procedure outlined in [Arianna](https://thedisorderedorganization.github.io/Arianna.jl/stable/#Usage) should be followed.
The first step involves initializing the system, represented by either an **Atoms** or **Molecules** structure. This requires specifying 
an initial configuration (including bonding information in the molecular case), the interaction potentials between different particle types, 
as well as the system’s density and temperature.
Subsequently, the desired *Arianna.jl* MCMC algorithm and the corresponding set of Monte Carlo moves are selected.
The next step consists of defining the simulation parameters, such as the total number of steps, the burn-in period, and the output directory.
Once these elements have been specified, the simulation can be executed.

This procedure can be carried out in two ways: either through the command-line interface (CLI) or programmatically within a Julia script, 
following the approach illustrated in *Arianna.jl* [examples](https://thedisorderedorganization.github.io/Arianna.jl/stable/#Usage). We detail here the CLI method, which is more user-friendly.

## Using the command-line interface

The CLI allows you to run simulations directly from the terminal. The basic command structure is as follows
```bash
particlesmc params.toml
```
where **params.toml** is a parameter file written in `TOML` format that specifies the simulation settings.
The `TOML` file must provide the following information:

- The *[system]* section specifies the input configuration file, the simulation temperature and density, and the type of neighbour list.
-  The *[model]* section defines the interaction potential between particles. 
-  The *[simulation]* section sets the simulation type to Metropolis, the number of Monte Carlo steps, the random seed, whether to run in parallel, and the output path.
-  The *[[simulation.move]]* section describes the Monte Carlo move to use.
- The *[[simulation.output]]* section configures the output (trajectories and callbacks).

## Example: Simulating a Lennard-Jones fluid
We perform a Monte Carlo simulation of a simple Lennard-Jones fluid using the package's CLI.
As explained above, we need to create a parameter file in TOML format that specifies the simulation settings.
An example of the parameter file is presented below for simulating a Lennard-Jones fluid at temperature $T=1.0$ and number density $\rho=1.0$.

```toml
[system]
config = "config.xyz"
temperature = 1.0
density = 1.0
list_type = "LinkedList"

[model]
[model."1-1"]
name = "LennardJones"
epsilon = 1.0
sigma = 1.0
rcut = 2.5

[simulation]
type = "Metropolis"
steps = 500   
output_path = "./"

[[simulation.move]]
action = "Displacement"
probability = 1.0
policy = "SimpleGaussian"
parameters = {sigma = 0.05}

[[simulation.output]]
algorithm = "StoreCallbacks"
callbacks = ["energy", "acceptance"]
scheduler_params = {linear_interval = 100}

[[simulation.output]]
algorithm = "StoreTrajectories"
scheduler_params = {linear_interval = 100}
fmt = "XYZ"
```

This example defines a minimal Monte Carlo simulation setup:
- The *[system]* The input configuration file is `config.xyz', the temperature is 1.0 and the density is 1.0. The neighbor list type is a linked list.
- The *[model]* The interaction potential is a Lennard-Jones potential is used for species pair `1--1` with specified parameters (`epsilon`, `sigma`, and cutoff `rcut`).
- The *[simulation]* The *Arianna.jl* MCMC method use is the Metropolis algorithm, the number of Monte Carlo steps is set to 500, and the output path is the simulation directory
- The *[[simulation.move]]* The Monte Carlo move is a displacement move with probability 1.0, guided by a simple Gaussian policy with a standard deviation (`sigma`) of 0.05.
- The *[[simulation.output]]* Two output algorithms are specified: one to store callbacks (energy and acceptance rate) every 100 steps, and another to store trajectories in XYZ format every 100 steps.

Executing the `particlesmc params.toml` command will run a basic Metropolis Monte Carlo simulation of particles 
interacting via the Lennard-Jones potential, using displacement moves, and periodically saving the system's trajectory. 
When the simulation ends, a file containing a summary of the simulation will be created in the output path.

This is a basic example of how to set up and run a Monte Carlo simulation using the package's CLI. Users can modify and extend this example to 
simulate a variety of atomic and molecular systems by changing the interactions, moves, and simtulation parameters.