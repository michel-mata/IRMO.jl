"set up!" |> println
num_workers = 16
include("./setup.jl")

# Game parameters
parameters = [1.0,0.2,0.02,0.02]

# Baseline
for norm in social_norms
    evolving = "interaction"
    # Sets of strategies
    Ms = [1,1,1]
    qs = [0,2,1]
    # Paths for results
    path = "results/endpoints/$evolving/$norm/baseline/"
    "\n\n"*path*"\n" |> println
    # Obtain attractors
    steady_states_detIC(path, parameters, norm, Ms, qs)
end

# Trios with unconditional strategies
for norm in social_norms, M in 2:2:10
    thresholds = [1,Int(M/2),M]
    ql = ["min","med","max"]

    for q in eachindex(thresholds)
        evolving = "interaction"
        # M=2 INTERIOR
        Ms = [1,1,M]
        qs = [0,2,thresholds[q]]
        # Paths for results
        path = "results/endpoints/$evolving/$norm/unconditional/M$M/$(ql[q])/"
        "\n\n"*path*"\n" |> println
        # Obtain attractors
        steady_states_detIC(path, parameters, norm, Ms, qs)
    end
end

# All thresholds - Trios with unconditional strategies
for norm in social_norms, M in 2:2:10
    for q in 1:M
        evolving = "interaction"
        Ms = [1,1,M]
        qs = [0,2,q]
        # Paths for results
        path = "results/endpoints/$evolving/$norm/unconditional-all/M$M/$q/"
        "\n\n"*path*"\n" |> println
        # Obtain attractors
        steady_states_detIC(path, parameters, norm, Ms, qs)
    end
end
