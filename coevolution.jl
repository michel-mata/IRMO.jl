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

begin
    using DifferentialEquations
    using SciMLNLSolve
    using NLsolve
    using LinearAlgebra
    using DelimitedFiles
    using Distributed
    using SharedArrays
    using ForwardDiff
end

# Game parameters
parameters = [1.0,0.2,0.02,0.02]

for norm in social_norms, m in 1:3
    M = [4,6,8][m]
    qM = [3,5,6]
    # Paths for results
    evolving = "coevolution"

    path = "results/trajectories/$evolving/$norm/2-$M/"
    mkpath(path)
    "\n\n"*path*"\n" |> println
    # Between thresholds and unconditional strategies
    Ms = [1,1,2,M]
    qs = [0,2,1,qM[m]]

    for i in 1:3
        freq = [repeat([1/S],S),[1/3,1/3,1/3,0],[1/3,1/3,0,1/3]][i]
        S = length(Ms)
        reps = reshape(repeat([1.0],S*S),(S,S))
        reputation, replicator = get_full_dynamics(norm,Ms,qs,parameters)
        replicator_dynamics = get_time_dynamics(norm,Ms,qs,parameters)
        T = 200
        @time trj = chop_pt.(solve(ODEProblem(replicator_dynamics,freq,(0.0,T),reps),abstol=1e-10,RadauIIA3(),saveat=1.0).u)
        trj = transpose(hcat(trj...))
        writedlm(path*"trj-$i-transient.csv", trj, ',')
    
        freq = [repeat([1/S],S),[1/3,1/3,1/3,0],[1/3,1/3,0,1/3]][i]
        S = length(Ms)
        reps = reshape(repeat([1.0],S*S),(S,S))
        reputation, replicator = get_full_dynamics(norm,Ms,qs,parameters)
        replicator_dynamics = get_time_dynamics(norm,Ms,qs,parameters)
        T = 200_000
        @time trj = chop_pt.(solve(ODEProblem(replicator_dynamics,freq,(0.0,T),reps),abstol=1e-10,RadauIIA5(),saveat=1e3).u)
        trj = transpose(hcat(trj...))
        writedlm(path*"trj-$i.csv", trj, ',')
    end
end