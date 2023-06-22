"set up!" |> println
num_workers = 16
# using Revise
include("./setup.jl")

# Game parameters
parameters = [1.0,0.2,0.02,0.02]

# Random samples
num_ic = 1_000

# Between numbers of observations >1, classic DISC, and unconditional strategies
for norm in social_norms
    evolving = "observations"
    Ms = [1,1,1,2:2:10...]
    thresholds = [[0,2,1,1,1,1,1,1],[0,2,1,1,2,3,4,5],[0,2,1,2:2:10...]]
    ql = ["min","med","max"]
    for q in eachindex(thresholds)
        # Paths for results
        path = "results/endpoints/$evolving/$norm/Ms-all/Q-$(ql[q])/"
        "\n\n"*path*"\n" |> println
        # Obtain attractors
        steady_states_rndIC(path, parameters, norm, Ms, thresholds[q], num_ic, refine=true)
    end
end