"set up!" |> println
num_workers = 16
include("./setup.jl")

# Pair of Scoring and Stern-Judging
begin
    # Parameters
    β = 0.1
    b, c, α, ϵ = 1.0,0.2,0.02,0.02
    # Sets of strategies
    n1 = social_norms[1]
    n2 = social_norms[2]
    Ns = [n1,n2]
    Ms = [1,1,2,2]
    Qs = [0,2,1,2]
    strategies = [[1,1],[0,0],[0,1],[0,1]]
    
    for η in [0.1,0.5,0.9], w in [0.1,0.5,0.9,1.0]
        parameters = [b, c, α, ϵ, β, η, w ]
        # Paths for results
        evolving = "norms"
        path = "results/endpoints/$evolving/independent/M2/$n1-$n2/$w-$η/"
        "\n\n"*path*"\n" |> println
        # Obtain attractors
        steady_states_normInd(path, parameters, Ns, Ms, Qs, strategies)
        get_cooperation_normInd(path, parameters, Ns, Ms, Qs, strategies)
    end
end

# All pairs of norms at equal initial proportions
for in1 in 1:length(social_norms)-1, in2 in (in1+1):length(social_norms)
    # Parameters
    β = 0.1
    η = 0.5
    b, c, α, ϵ = 1.0,0.2,0.02,0.02
    # Sets of strategies
    n1 = social_norms[in1]
    n2 = social_norms[in2]
    Ns = [n1,n2]
    Ms = [1,1,2,2]
    Qs = [0,2,1,2]
    strategies = [[1,1],[0,0],[0,1],[0,1]]
    
    for w in [0.1,0.5,0.9,1.0]
        parameters = [b, c, α, ϵ, β, η, w ]
        # Paths for results
        evolving = "norms"
        path = "results/endpoints/$evolving/independent/M2/$n1-$n2/$w-$η/"
        "\n\n"*path*"\n" |> println
        # Obtain attractors
        steady_states_normInd(path, parameters, Ns, Ms, Qs, strategies)
        get_cooperation_normInd(path, parameters, Ns, Ms, Qs, strategies)
    end
end