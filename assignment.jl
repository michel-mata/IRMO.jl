"set up!" |> println
num_workers = 9
# using Revise
include("./setup.jl")

# Game parameters
parameters = [1.0,0.2,0.02,0.02]


# Fixed M and 3 thresholds
num_ic = 1_000
for norm in social_norms, M in 2:2:10
    # Paths for results
    evolving = "assignment"
    path = "results/endpoints/$evolving/$norm/Qs-alone/$M/"
    "\n\n"*path*"\n" |> println
    # Between thresholds and unconditional strategies
    if M==2
        Ms = [M,M]
        qs = [1,M]
        # Obtain attractors
        steady_states_detIC(path, parameters, norm, Ms, qs)
    else
        Ms = [M,M,M]
        qs = [1,Int(M/2),M]
        # Obtain attractors
        steady_states_detIC(path, parameters, norm, Ms, qs)
    end
end

# Fixed M and 3 thresholds with unconditional
num_ic = 1_000
for norm in social_norms, M in 2:2:10
    # Paths for results
    evolving = "assignment"
    path = "results/endpoints/$evolving/$norm/Qs-uncond/$M/"
    "\n\n"*path*"\n" |> println
    # Between thresholds and unconditional strategies
    if M==2
        Ms = [1,M,M,1]
        qs = [0,1,M,2]
        # Obtain attractors
        steady_states_detIC(path, parameters, norm, Ms, qs)
    else
        Ms = [1,M,M,M,1]
        qs = [0,1,Int(M/2),M,2]
        # Obtain attractors
        steady_states_rndIC(path, parameters, norm, Ms, qs, num_ic, refine=true)
    end
end

# Fixed M and all thresholds
num_ic = 100
for norm in social_norms, M in 2:2:8
    # Paths for results
    evolving = "assignment"
    path = "results/endpoints/$evolving/$norm/Qs-opt/$M/"
    "\n\n"*path*"\n" |> println
    # Between thresholds and unconditional strategies
    Ms = [1,1,repeat([M],M)...]
    qs = [0,2,1:M...]
    # Obtain attractors
    steady_states_rndIC(path, parameters, norm, Ms, qs, num_ic, refine=true)
end