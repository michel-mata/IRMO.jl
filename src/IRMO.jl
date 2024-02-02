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

    # Auxiliar functions
    include("./module/aux.jl")
    export chop_pt
    export get_norm
    export get_initial_frequencies
    # Strategy competition
    include("./module/fns-strats.jl")
    export get_full_dynamics
    export get_time_dynamics
    export get_rep
    export trajectories_detIC
    export steady_states_detIC
    export steady_states_rndIC
    export refine_steady_states
    export get_cooperation
    # Norm competition
    include("./module/fns-norms.jl")
    export get_full_dynamics_n
    export get_time_dynamics_n
    export get_rep_n
    export steady_states_detIC_n
    export steady_states_rndIC_n
    export refine_steady_states_n
    export get_cooperation_n
    # Strategy competition under fixed set of norms
    include("./module/fns-strats-fix-norms.jl")
    export get_full_dynamics_nf
    export get_time_dynamics_nf
    export get_rep_nf
    export steady_states_detIC_nf
    export steady_states_rndIC_nf
    export refine_steady_states_nf
    export get_cooperation_nf
    # Strategy competition with probabilistic cooperator
    include("./module/fns-probc.jl")
    export get_full_dynamics_s
    export get_time_dynamics_s
    export get_rep_s
    export steady_states_detIC_s
    export steady_states_rndIC_s
    export refine_steady_states_s
    export get_cooperation_s


    # Auxiliary constants:
    export social_norms
    const social_norms = ["SC","SJ","SS","SH"]

    export grid
	const grid = 0.05:0.05:0.95

end