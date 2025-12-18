# Molecular Monte Carlo

To see an introduction to Markov Chain Monte Carlo, please read *Arianna*'s [documentation](https://thedisorderedorganization.github.io/Arianna.jl/stable/man/montecarlo/).

## Birth of Molecular MCMC

The pioneering application of the Metropolis algorithm by Metropolis *et al.*[^1] focused on sampling the equilibrium distribution of a system consisting of $N = 224$ hard disks. A microstate of the system is fully characterized by the set of particle positions in space, $x = {\{ \mathbf{r}_j \}}_{j=1}^N$ 
where $\mathbf{r}_j$ denotes the position of particle $j$. Consequently, the configuration space is the $2N$-dimensional space of particle coordinates.

The target distribution $P$ corresponds to the Boltzmann distribution, whose probability density is given by:

```math
p_B(x) \propto \exp\left[-\frac{U(x)}{T}\right],
```

where $U(x)$ is the potential energy associated with configuration $x$, and $T$ denotes the temperature of the system.

The proposed Monte Carlo move in their study consists of sequentially displacing the disks by a random vector with uniformly distributed magnitude and direction.
Each particle displacement from $x$ to $x'$ is accepted or rejected according to the Metropolis acceptance probability,


```math
\alpha(x,x') = \min \left\{ 1, \exp\left[-\frac{\Delta U}{T}\right] \right\},
```

where $\Delta U = U(x') - U(x)$ represents the change in potential energy resulting from the proposed displacement.
In this hard disk case, the difference is either infinity or $0$ depending if the displaced disk overlaps with another disk or not.
We directly see that detailed balance is not respected as particles are not selected at random but sequentially, 
meaning that $q(x, x') \neq q(x', x)$.

That is why this displacement move has later been modified in two steps, such that it respects detailed balance: 
1. Select a particle uniformly at random from the system.
2. Displace the chosen particle by a random vector with uniformly distributed magnitude and direction.

This simple move type subsequently became the standard Monte Carlo move in liquid-state physics and remains widespread today. 
Notably, the original simulation determined numerically the equation of state for this system.

The simulation was performed on the Los Alamos MANIAC I computer (Mathematical Analyzer Numerical Integrator and Automatic Computer Model I). 
The researchers conducted a burn-in (equilibration) phase of 16 sweeps (where one sweep consists of $N$ attempted moves), followed by a production run of 48 to 64 sweeps. 
Each sweep required approximately 3 minutes of computation time, resulting in a total runtime of about 4 hours. By contrast, the same calculation would take approximately 
10 milliseconds on a modern laptop, representing a performance improvement of roughly six orders of magnitude.

```@raw html
<h1 align="center">
  <img src="https://raw.githubusercontent.com/TheDisorderedOrganization/ParticlesMC/main/docs/src/assets/maniac.jpg" width="500"/>
</h1>

<p align="center"> The MANIAC I computer at Los Alamos National Laboratory, where Arianna Rosenbluth implemented the first MCMC algorithm. 
This machine represented the cutting edge of computational capability in the 1950s. Image taken from Ref.~\cite{maniac}.</p>
```

## Monte Carlo Moves: Actions and Policies

The transition from state $x$ to state $x'$ in the MH algorithm is commonly referred to as a *Monte Carlo move*. Monte Carlo moves can be conceptually decomposed into two components[^2]:

- **Actions:** The type of modification applied to the system (e.g., particle displacement, rotation, swap, flip).
- **Policies:** The probabilistic rules governing how the action parameters are sampled (e.g., the distribution of displacement magnitudes and directions). Policies can be symmetric or non-symmetric. If a policy is symmetric, then the Metropolis algorithm is employed, else the Metropolis-Hastings algorithm is employed.


To illustrate this distinction, consider the move employed in the original Metropolis simulation. The *action* is the translation of a randomly selected particle. 
To satisfy the detailed balance condition with the symmetric Metropolis acceptance rule, the translation's direction and magnitude can
be sampled from a symmetric distribution, for example uniform, Gaussian, exponential, etc. The choice of this distribution constitutes the *policy*.

Transitions between successive states, $x_i \to x_{i+1}$, need not be governed by a single type of Monte Carlo move. In practice, it is often advantageous to define a pool of 
distinct move types and to select one randomly from this pool at each step of the Markov chain. To each move type $m$ in the pool, one assigns a selection probability 
$p_m$, with the normalization condition $\sum_{m} p_m = 1$.
This collection of moves, together with their associated probabilities, can be regarded as a single composite Monte Carlo move.

## Choosing the right parameters

```@raw html
<h1 align="center">
  <img src="https://raw.githubusercontent.com/TheDisorderedOrganization/ParticlesMC/main/docs/src/assets/goodpolicy.pdf" width="500"/>
</h1>

<p align="center"> Markov chain state space exploration for three different displacement policies. The system is a 1D particle in a potential well $U(x)$, with the target distribution being the Boltzmann distribution. The proposed move is the following: the action is the displacement of the particle and the policy is the distribution of the displacement magnitude and direction. Here, the policy is an uniform law $\sim \mathcal{U}(-u,u)$ with $u>0$. Each dot represents a state in the Markov chain, constructed from bottom to top. The horizontal position of a dot is the position of the particle. <b>Center:</b> When $u$ is too small, all moves are accepted but the state space is explored inefficiently due to small step sizes. <b>Right:</b> When $u$ is too large, nearly all moves are rejected, again leading to poor exploration. <b>Left:</b> An optimal choice of $u$ balances acceptance rate and step size, achieving efficient state space sampling. This illustrates the importance of tuning the proposal distribution.</p>
```

A MCMC simulation consists of two phases: a *burn-in* phase and an *equilibrium* phase[^3]. During burn-in, initial samples are discarded because the chain has not yet reached its stationary distribution. Starting from an arbitrary initial state, early samples may overweight 
low-probability regions of the state space. The burn-in period must be long enough to allow the chain to converge to the target distribution before collecting samples for analysis.

Once the chain reaches equilibrium, samples are collected to estimate properties of the system. The simulation must run long enough to 
adequately explore the state space $\mathcal{X}$. However, successive samples are typically correlated, reducing the effective sample size. To assess sampling quality, one should examine autocorrelation functions and apply convergence diagnostics to verify that the chain has sufficiently explored $\mathcal{X}$[^4].

The choice of actions and policies is critical for efficient sampling. As shown in the above Figure, a poorly chosen policy leads to 
inefficient exploration of the state space. Finding optimal actions and policies are essential for effective MCMC sampling.

[^1]: *Equation of State Calculations by Fast Computing Machines*, Metropolis, Nicholas and Rosenbluth, Arianna W. and Rosenbluth, Marshall N. and Teller, Augusta H. and Teller, Edward, 1953
[^2]: *Reinforcement learning: An introduction*, Sutton, Richard S and Barto, Andrew G and others, 1998
[^3]: *Monte Carlo methods*, Guiselin, Benjamin, 2024
[^4]: *Convergence diagnostics for markov chain monte carlo*, Roy, Vivekananda, 2020
