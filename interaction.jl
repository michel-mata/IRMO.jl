"set up!" |> println
num_workers = 16
include("./setup.jl")

# Game parameters
parameters = [1.0,0.2,0.02,0.02]

# Baseline
for norm in social_norms
    # Paths for results
    evolving = "interaction"
    path = "results/endpoints/$evolving/$norm/baseline/"
    "\n\n"*path*"\n" |> println
    # Sets of strategies
    Ms = [1,1,1]
    qs = [0,2,1]
    # Obtain attractors
    steady_states_detIC(path, parameters, norm, Ms, qs)
    get_cooperation(path, parameters, norm, Ms, qs)
end

# Trios with unconditional strategies
for norm in social_norms, M in 2:2:10
    thresholds = [1,Int(M/2),M]
    ql = ["min","med","max"]

    for q in eachindex(thresholds)
        # Paths for results
        evolving = "interaction"
        path = "results/endpoints/$evolving/$norm/unconditional/M$M/$(ql[q])/"
        "\n\n"*path*"\n" |> println
        # M=2 INTERIOR
        Ms = [1,1,M]
        qs = [0,2,thresholds[q]]
        # Obtain attractors
        steady_states_detIC(path, parameters, norm, Ms, qs)
        get_cooperation(path, parameters, norm, Ms, qs)
    end
end

# Fixed M and 3 thresholds with unconditional
for norm in social_norms, M in 2:2:10
    # Paths for results
    evolving = "assignment"
    path = "results/endpoints/$evolving/$norm/$M/"
    "\n\n"*path*"\n" |> println
    # Between thresholds and unconditional strategies
    if M==2
        Ms = [1,M,M,1]
        qs = [0,1,M,2]
        # Obtain attractors
        steady_states_detIC(path, parameters, norm, Ms, qs)
    else
        num_ic = 1_000
        Ms = [1,M,M,M,1]
        qs = [0,1,Int(M/2),M,2]
        # Obtain attractors
        steady_states_rndIC(path, parameters, norm, Ms, qs, num_ic, refine=true)
    end
    get_cooperation(path, parameters, norm, Ms, qs)
end