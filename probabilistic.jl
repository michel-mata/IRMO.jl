"set up!" |> println
num_workers = 16
include("./setup.jl")

# Game parameters
parameters = [1.0,0.2,0.02,0.02]
ps = 0.1:0.2:0.9


# Baseline
for norm in social_norms[1:2], p in [0.,ps...,1.0]
    evolving = "probabilistic-cooperator"
    # Sets of strategies
    Ms = [1,1]
    qs = [0,1]
    strategies = [[p,p],[0,1]]
    # Paths for results
    path = "results/endpoints/$evolving/$norm/baseline/PROB-DISC/$p/"
    "\n\n"*path*"\n" |> println
    # Obtain attractors
    steady_states_detIC_s(path, parameters, norm, Ms, qs, strategies)
    get_cooperation_s(path, parameters, norm, Ms, qs, strategies)
end

# LTFO
for norm in social_norms[1:2], p in [0.,ps...,1.0]
    evolving = "probabilistic-cooperator"
    # Sets of strategies
    Ms = [1,2]
    qs = [0,1]
    strategies = [[p,p],[0,1]]
    # Paths for results
    path = "results/endpoints/$evolving/$norm/2/PROB-LTFO/$p/"
    "\n\n"*path*"\n" |> println
    # Obtain attractors
    steady_states_detIC_s(path, parameters, norm, Ms, qs, strategies)
    get_cooperation_s(path, parameters, norm, Ms, qs, strategies)
end

# LTFN
for norm in social_norms[1:2], p in [0.,ps...,1.0]
    evolving = "probabilistic-cooperator"
    # Sets of strategies
    Ms = [1,2]
    qs = [0,2]
    strategies = [[p,p],[0,1]]
    # Paths for results
    path = "results/endpoints/$evolving/$norm/2/PROB-LTFN/$p/"
    "\n\n"*path*"\n" |> println
    # Obtain attractors
    steady_states_detIC_s(path, parameters, norm, Ms, qs, strategies)
    get_cooperation_s(path, parameters, norm, Ms, qs, strategies)
end

########### (PROBC, ALLD, DISC) ###########
for norm in social_norms[1:2], p in ps
    evolving = "probabilistic-cooperator"
    # Sets of strategies
    Ms = [1,1,1]
    qs = [0,2,1]
    strategies = [[p,p],[0,0],[0,1]]
    # Paths for results
    path = "results/endpoints/$evolving/$norm/baseline/PROB-ALLD-DISC/$p/"
    "\n\n"*path*"\n" |> println
    # Obtain attractors
    steady_states_detIC_s(path, parameters, norm, Ms, qs, strategies)
    get_cooperation_s(path, parameters, norm, Ms, qs, strategies)
end

########### (PROBC, ALLD, DISC_2,1, DISC_2,2) ###########
for norm in social_norms[1:1], p in ps
    evolving = "probabilistic-cooperator"
    # Sets of strategies
    Ms = [1,1,2,2]
    qs = [0,0,1,2]
    strategies = [[p,p],[0,0],[0,1],[0,1]]
    # Paths for results
    path = "results/endpoints/$evolving/$norm/2/PROB-ALLD-LTFO-LTFN/$p/"
    "\n\n"*path*"\n" |> println
    # Obtain attractors
    steady_states_detIC_s(path, parameters, norm, Ms, qs, strategies)
    get_cooperation_s(path, parameters, norm, Ms, qs, strategies)
end
