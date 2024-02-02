"set up!" |> println
num_workers = 16
include("./setup.jl")

# Game parameters
parameters = [1.0,0.2,0.02,0.02]

# Baseline-pairwise
for in1 in 1:length(social_norms)-1, in2 in (in1+1):length(social_norms)
    n1 = social_norms[in1]
    n2 = social_norms[in2]
    evolving = "norms"
    # Sets of strategies
    Ms = [1,1,1,1]
    qs = [0,2,1,1]
    # Social norms competing
    norms = ["ALLC","ALLC",n1,n2]
    # Paths for results
    path = "results/endpoints/$evolving/pairwise/baseline/$n1-$n2/"
    "\n\n"*path*"\n" |> println
    # Obtain attractors
    steady_states_detIC_n(path, parameters, norms, Ms, qs, refine=true)
    get_cooperation_n(path, parameters, norms, Ms, qs)
end

# Baseline-all
begin
    evolving = "norms"
    # Sets of strategies
    Ms = repeat([1],length(social_norms)+2)
    qs = [0,2,repeat([1],length(social_norms))...]
    # Social norms competing
    norms = ["ALLC","ALLD",social_norms...]
    # Paths for results
    path = "results/endpoints/$evolving/unconditional/baseline/"
    "\n\n"*path*"\n" |> println
    # Obtain attractors
    num_ic = 1_000
    steady_states_rndIC_n(path, parameters, norms, Ms, qs, num_ic, refine=true)
    get_cooperation_n(path, parameters, norms, Ms, qs)
end

# M2-pairwise
for in1 in 1:length(social_norms)-1, in2 in (in1+1):length(social_norms)
    n1 = social_norms[in1]
    n2 = social_norms[in2]
    evolving = "norms"
    # Sets of strategies
    Ms = [1,1,2,2]
    qs = [0,2,1,1]
    # Social norms competing
    norms = ["ALLC","ALLD",n1,n2]
    # Paths for results
    path = "results/endpoints/$evolving/pairwise/M2/$n1-$n2/"
    "\n\n"*path*"\n" |> println
    # Obtain attractors
    steady_states_detIC_n(path, parameters, norms, Ms, qs, refine=true)
    get_cooperation_n(path, parameters, norms, Ms, qs)
end

# M2-all
begin
    evolving = "norms"
    # Sets of strategies
    Ms = repeat([2],length(social_norms)+2)
    qs = [0,3,repeat([1],length(social_norms))...]
    # Social norms competing
    norms = ["ALLC","ALLD",social_norms...]
    # Paths for results
    path = "results/endpoints/$evolving/unconditional/M2/"
    "\n\n"*path*"\n" |> println
    # Obtain attractors
    num_ic = 1_000
    steady_states_rndIC_n(path, parameters, norms, Ms, qs, num_ic, refine=true)
    get_cooperation_n(path, parameters, norms, Ms, qs)
end


# Fixed norms
parameters = [1.0,0.2,0.02,0.02]

# Baseline
for in1 in 1:length(social_norms)-1, in2 in (in1+1):length(social_norms)
    n1 = social_norms[in1]
    n2 = social_norms[in2]
    for p in [0.0,0.1,0.5,0.9,1.0]
        evolving = "norms"
        # Sets of strategies
        Ms = [1,1,1]
        qs = [0,2,1]
        norms =norms = [n1,n2]
        n = [p,1-p]
        # Paths for results
        path = "results/endpoints/$evolving/fixed/baseline/$n1-$n2/$p/"
        "\n\n"*path*"\n" |> println
        # Obtain attractors
        steady_states_detIC_nf(path, parameters, norms, Ms, qs, n, refine=true)
        get_cooperation_nf(path, parameters, norms, Ms, qs, n)
    end
end

# M2
for in1 in 1:length(social_norms)-1, in2 in (in1+1):length(social_norms)
    n1 = social_norms[in1]
    n2 = social_norms[in2]
    for p in [0.0,0.1,0.5,0.9,1.0]
        evolving = "norms"
        # Sets of strategies
        Ms = [1,1,2,2]
        qs = [0,2,1,2]
        norms = [n1,n2]
        n = [p,1-p]
        # Paths for results
        path = "results/endpoints/$evolving/fixed/M2/$n1-$n2/$p/"
        "\n\n"*path*"\n" |> println
        # Obtain attractors
        steady_states_detIC_nf(path, parameters, norms, Ms, qs, n, refine=true)
        get_cooperation_nf(path, parameters, norms, Ms, qs, n)
    end
end