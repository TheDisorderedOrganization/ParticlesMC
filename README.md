<h1 align="center">
  <img src="particlesmc_logo.png" width="500"/>
</h1>

<div align="center">

  [![license](https://img.shields.io/badge/license-GPL%203.0-red.svg)](https://github.com/TheDisorderedOrganization/ParticlesMC/blob/main/LICENSE)
  [![ci](https://github.com/TheDisorderedOrganization/ParticlesMC/actions/workflows/ci.yml/badge.svg)](https://github.com/TheDisorderedOrganization/ParticlesMC/actions/workflows/ci.yml)
  [![codecov](https://codecov.io/gh/TheDisorderedOrganization/ParticlesMC/graph/badge.svg?token=URGL1HJOOI)](https://codecov.io/gh/TheDisorderedOrganization/ParticlesMC)
</div>

<p align="center">
ParticlesMC is a Julia package for performing atomic and molecular Monte Carlo simulations. It is designed to be efficient and user-friendly, making it suitable for both research and educational purposes. Built on top of the <a href="https://github.com/TheDisorderedOrganization/Arianna.jl"> Arianna </a> module, it leverages Arianna’s Monte Carlo framework.
</p>

<p align="center">
  <img src="https://github.com/TheDisorderedOrganization/ParticlesMC/blob/main/data/movie/movie.gif" alt="MC simulation of a 2D liquid" width="600">
  This example can be reproduced by running `particlesmc params.toml` in the `data/movie/` folder . Movie generated with [OVITO](https://www.ovito.org/)
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
After the 1st release, you will be able to install ParticlesMC using the Julia package manager in one of two ways:

1. Using the package mode (press `]` in the Julia REPL):
```julia
add ParticlesMC
```

2. Using the Pkg API:
```julia
using Pkg
Pkg.add("ParticlesMC")
```

### Building ParticlesMC

I you want to build `ParticlesMC`, enter the package mode (press `]` in the Julia REPL) and type:
```julia
build
```
This will build the `particlesmc` executable at `~/.julia/bin` (please add this path to your PATH).

## Usage

### Running a Monte Carlo Simulation

To run a Monte Carlo simulation, you need to define the systems and moves, and then use the provided functions to execute the simulation. Here is a basic example:

Add example code
```julia
```

## Contributing

We welcome contributions from the community. If you have a new system or feature to add, please fork the repository, make your changes, and submit a pull request.

## Citing

If you use Arianna in your research, please cite it! You can find the citation information in the [CITATION](https://github.com/TheDisorderedOrganization/ParticlesMC/blob/main/CITATION.bib) file or directly through GitHub's "Cite this repository" button.

## License

This project is licensed under the GNU General Public License v3.0.  License. See the [LICENSE](https://github.com/TheDisorderedOrganization/ParticlesMC/blob/main/LICENSE) file for details.

## Contact

For any questions or issues, please open an issue on the GitHub repository or contact the maintainers.
