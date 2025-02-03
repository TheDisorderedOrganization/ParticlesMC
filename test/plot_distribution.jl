# TEST WITH GERHARD'S CONFIG
using Plots
using DelimitedFiles

temperature = 0.231
N = 1290
pswap = 0.0
M = 1
path_empty = "data/test/particles/KA2D_distribution/N$N/T$temperature/pswap$pswap/M$M"
path_linked = "data/test/particles/KA2D_distribution_LL/N$N/T$temperature/pswap$pswap/M$M"

energy_empty = readdlm(joinpath(path_empty, "energy.dat"))[:, 2]
energy_linked = readdlm(joinpath(path_linked, "energy.dat"))[:, 2]


@assert all(isapprox.(energy_empty, energy_linked, atol=1e-10))

plot()
plot!(energy_empty)
plot!(energy_linked)
plot(yaxis=(:log10, [0.1, :auto]))
stephist!(-energy_empty)
stephist!(-energy_linked)