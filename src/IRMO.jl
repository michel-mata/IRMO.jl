module IRMO

    # Packages
    using DifferentialEquations
    using SciMLNLSolve
    using NLsolve
    using ForwardDiff
    using LinearAlgebra
    using DelimitedFiles
    using Distributed
    using SharedArrays

    # Objects:
    include("./module/functions.jl")
    export chop_pt
    export get_norm
    export get_rep
    export get_initial_frequencies
    export get_vertex, get_eigenvalues
    export get_full_dynamics
    export get_time_dynamics
    export steady_states_detIC
    export steady_states_rndIC
    export trajectories_detIC
    export refine_steady_states

    # Auxiliary constants:
    export social_norms
    const social_norms = ["SC","SJ","SS","SH"]

    export grid
	const grid = 0.05:0.05:0.95

end