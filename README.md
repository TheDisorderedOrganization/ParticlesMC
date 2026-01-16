<h1 align="center">
  <img src="https://raw.githubusercontent.com/TheDisorderedOrganization/ParticlesMC/main/logo.png" width="500"/>
</h1>

<div align="center">

  [![docs](https://img.shields.io/badge/docs-online-blue.svg)](https://TheDisorderedOrganization.github.io/ParticlesMC)
  [![license](https://img.shields.io/badge/license-GPL%203.0-red.svg)](https://github.com/TheDisorderedOrganization/ParticlesMC/blob/main/LICENSE)
  [![ci](https://github.com/TheDisorderedOrganization/ParticlesMC/actions/workflows/ci.yml/badge.svg)](https://github.com/TheDisorderedOrganization/ParticlesMC/actions/workflows/ci.yml)
  [![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
  [![codecov](https://codecov.io/gh/TheDisorderedOrganization/ParticlesMC/graph/badge.svg?token=URGL1HJOOI)](https://codecov.io/gh/TheDisorderedOrganization/ParticlesMC)

</div>

<p align="center">
ParticlesMC is a Julia package for performing atomic and molecular Monte Carlo simulations. It is designed to be efficient and user-friendly, making it suitable for both research and educational purposes. Built on top of the <a href="https://github.com/TheDisorderedOrganization/Arianna.jl"> Arianna </a> module, it leverages Arianna’s Monte Carlo framework.
</p>

<p align="center">
  <img src="https://raw.githubusercontent.com/TheDisorderedOrganization/ParticlesMC/main/examples/movie/movie.gif" alt="MC simulation of a 2D liquid" width="400">
</p>

<p align="center">
  MC simulation of a 2D liquid. This example can be reproduced by running <code>particlesmc params.toml</code> in the <code>examples/movie/</code> folder. Movie generated with <a href="https://www.ovito.org/"> ovito</a>.
</p>

## Features

- **Flexible execution**: Simulations can be run from standalone scripts or through a command-line interface (CLI), enabling easy integration into workflows.
- **Interaction potentials**: Supports a broad range of interatomic and intermolecular interaction potentials.
- **Monte Carlo moves**: Implements state-of-the-art Monte Carlo moves for both atomic and molecular simulations.
- **Computational efficiency**: Designed with performance in mind to enable fast simulations.
- **Arianna framework integration**: Leverages the [Arianna](https://github.com/TheDisorderedOrganization/Arianna.jl) Monte Carlo framework, benefiting from advanced techniques such as Policy-Guided Monte Carlo (PGMC) and parallel tempering (soon).


## Installation

### Requirements
- Julia version 1.9 or higher

### Installing ParticlesMC
You can install ParticlesMC using the Julia package manager in one of two ways:

1. Using the package mode (press `]` in the Julia REPL):
```julia
add https://github.com/TheDisorderedOrganization/ParticlesMC.git
```

2. Using the Pkg API:
```julia
using Pkg
Pkg.add(url="https://github.com/TheDisorderedOrganization/ParticlesMC.git")
```

### Building ParticlesMC

The build should be automatic when installing `ParticlesMC`. If it hasn't, you can manually build the package, by entering the package mode (press `]` in the Julia REPL) and by typing:

```julia
build
```
This will build the `particlesmc` executable at `~/.julia/bin` (please add this path to your PATH).

## Usage

### Running a Monte Carlo Simulation

To run a Monte Carlo simulation, you need an input atomic or molecular configuration file (typically with a `.xyz` extension) and a parameter file (in `TOML` format). The parameter file specifies both the system details (such as temperature, density, and interaction model) and the simulation details (such as simulation type, number of steps, Monte Carlo moves, and outputs). A minimal example is presented below. More detailed explanations can be found in the documentation.

**config.xyz**
```
3
Lattice="1.7321 0.0 0.0 0.0 1.7321 0.0 0.0 0.0 0.0" Properties=type:I:1:pos:R:2
1 0.1585 0.4965
1 1.7215 0.7468
1 0.7606 1.1439
```

**params.toml**
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
seed = 10
parallel = false
output_path = "./"

[[simulation.move]]
action = "Displacement"
probability = 1.0
policy = "SimpleGaussian"
parameters = {sigma = 0.05}

[[simulation.output]]
algorithm = "StoreTrajectories"
scheduler_params = {linear_interval = 50}
fmt = "XYZ"
```

**Explanation of the example:**

This example defines a minimal Monte Carlo simulation setup:

- The `[system]` section specifies the input configuration file (`config.xyz`), the simulation temperature and density, and the use of a linked list for neighbor searching.
- The `[model]` section defines the interaction potential between particles. Here, a Lennard-Jones potential is used for species pair "1-1" with specified parameters (`epsilon`, `sigma`, and cutoff `rcut`).
- The `[simulation]` section sets the simulation type to Metropolis, the number of Monte Carlo steps to 500, the random seed, whether to run in parallel, and the output path.
- The `[[simulation.move]]` section describes the Monte Carlo move to use: a displacement move with probability 1.0, guided by a simple Gaussian policy with a standard deviation (`sigma`) of 0.05.
- The `[[simulation.output]]` section configures the output: trajectories will be stored every 50 steps in the XYZ format.

By executing `particlesmc params.toml** you will run a basic Metropolis Monte Carlo simulation of particles interacting via the Lennard-Jones potential, using displacement moves, and periodically saving the system's trajectory.

**Extending beyond this simple example:**

More models and moves are available.

For instance, in a simulation that contains multiple species, __DiscreteSwap__ moves can be added to swap species between to random particles:
```toml
[[simulation.move]]
action = "DiscreteSwap"
probability = 0.1
policy = "DoubleUniform"
parameters = {species = [1, 2]}
```
is a move that will be attempted 10% of the time (probability = 0.1). It will select a random particle _i_ of species _1_, a random particle _j_ of species _2_, and will attempt to transform _i_ into species _2_, and _j_ into species _1_.

## Contributing

We welcome contributions from the community. If you have a new system or feature to add, please fork the repository, make your changes, and submit a pull request.

## Citing

If you use Arianna in your research, please cite it! You can find the citation information in the [CITATION](https://github.com/TheDisorderedOrganization/ParticlesMC/blob/main/CITATION.bib) file or directly through GitHub's "Cite this repository" button.

## License

This project is licensed under the GNU General Public License v3.0.  License. See the [LICENSE](https://github.com/TheDisorderedOrganization/ParticlesMC/blob/main/LICENSE) file for details.

## Contact

For any questions or issues, please open an issue on the GitHub repository or contact the maintainers.
