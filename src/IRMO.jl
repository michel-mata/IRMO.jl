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

    include("./module/functions.jl")
    # Auxiliar functions
    export chop_pt
    export get_norm
    export get_rep
    export get_initial_frequencies
    export get_vertex, get_eigenvalues
    # Strategy competition
    export get_full_dynamics
    export get_time_dynamics
    export steady_states_detIC
    export steady_states_rndIC
    export trajectories_detIC
    export refine_steady_states
    # Norm competition
    export get_full_dynamics_n
    export get_time_dynamics_n
    export steady_states_detIC_n
    export steady_states_rndIC_n
    export refine_steady_states_n

    # Auxiliary constants:
    export social_norms
    const social_norms = ["SC","SJ","SS","SH"]

    export grid
	const grid = 0.05:0.05:0.95

end