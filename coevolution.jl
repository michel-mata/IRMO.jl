"set up!" |> println
num_workers = 9
# using Revise
include("./setup.jl")

# Game parameters
parameters = [1.0,0.2,0.02,0.02]

# Stable mixtures
num_ic = 100
for norm in social_norms, m in 1:3
    M = [4,6,8][m]
    ms = [M,M]
    qM = [[2,3],[4,5],[5,6]]
    # Paths for results
    evolving = "coevolution"
    path = "results/endpoints/$evolving/$norm/2-$M/"
    "\n\n"*path*"\n" |> println
    # Between thresholds and unconditional strategies
    Ms = [1,1,2,2,ms...]
    qs = [0,2,1,2,qM[m]...]
    # Obtain attractors
    steady_states_rndIC(path, parameters, norm, Ms, qs, num_ic, refine=true)
end