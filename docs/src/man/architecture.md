
# Package Architecture

The package adopts a modular architecture, with each module responsible for a distinct component of the simulation workflow. 

## Core systems

Two fundamental data structures represent physical systems: the **Atoms** struct for atomic systems and the **Molecules** struct for molecular systems. 
Both structures inherit from the abstract **Particles** type, which itself is a subtype of **AriannaSystem** (see *Arianna* [types](https://thedisorderedorganization.github.io/Arianna.jl/stable/man/system/)). 
This design allows the use of the [*Arianna.jl*](https://github.com/TheDisorderedOrganization/Arianna.jl) framework. The **Atoms** struct encapsulates properties such as particle positions, 
atom types, and interaction potentials. The **Molecules** struct extends the **Atoms**  struct by including bond information (type and interaction potential).

Users can choose from a variety of potentials for calculating interactions between particles. The interactions depend solely on the type of the particles.
For molecular systems, users choose potentials for both bonded particles and non-bonded particles. The interactions depend on the type of particles for 
non-bonded interactions and on the type of bond for bonded interactions.

In both cases, the target distribution is the Boltzmann distribution and thus, 
simulations are performed in the canonical ensemble (constant number of particles, volume, and temperature). An extension to perform
simulations in the NPT ensemble (constant number of particles, pressure, and temperature) is in current development.

## Monte Carlo Moves

To use the [*Arianna.jl*](https://github.com/TheDisorderedOrganization/Arianna.jl) framework, the package implements moves defined through actions and policies, as described in the previous [section](https://thedisorderedorganization.github.io/ParticlesMC/stable/man/particlesmc). 
The package provides a collection of predefined actions and policies corresponding to standard and state-of-the-art Monte Carlo moves in atomic and molecular simulations:

- Standard displacement moves.
- Swap moves for particle exchange in atomic systems.
- Flip moves in molecular systems.
- Intramolecular swap moves in polymeric systems.

Each action supports multiple policy implementations. For instance, displacement moves can employ either normal 
or uniform probability distributions for trial move generation.

All moves are designed to be compatible with both the Metropolis-Hastings algorithm and the Policy-Guided Monte Carlo (PGMC) method of Arianna.
This means that all functions required by the *Arianna.jl* framework (see *Arianna* [moves](https://thedisorderedorganization.github.io/Arianna.jl/stable/man/system/)) are implemented for each action and policy.

## Neighbour Lists

The package incorporates efficient neighbor list algorithms[^1] to accelerate interaction computations between particles.
Users can choose two neighbor list implementations, **LinkedList** and **EmptyList**. Both structs are subtypes of the 
**NeighbourList** type. The **EmptyList** is a dummy neighbour list that does not store any information, and thus 
the interaction between particles is computed by looping over all the particles. The **LinkedList** struct is 
a cell list implementation that divides the simulation box into cells of size larger that the maximul cutoff of the interaction
potentials. This permits to calculate the potential energy of a particle by only calculating its interactions with particles 
inside the neighbouring cells. For large systems, using the **LinkedList** struct accelerates significantly simulations 
compared to using the **EmptyList** struct. 


## Input/Output

The package handle file operations, supporting configuration import from *.xyz*, *.exyz*, and *.lammps* formats.
The package can also write the output trajectories in these same formats, providing users flexibility in choosing their preferred file format 
for data storage and analysis.

Additionally, the package provides callback functions for on-the-fly calculation of thermodynamic and structural properties during simulations, 
including energy and pressure for both atomic and molecular systems. These callbacks are stored in the output directory for post-processing.

[^1]: *Computer simulation of liquids*, Allen, Michael P and Tildesley, Dominic J, 2017