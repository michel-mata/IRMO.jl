# Load packages
include("./setup.jl")
# Game parameters
parameters = [1.0,0.2,0.02,0.02]
# Social norm
norm = "SJ"
# Sets of strategies
Ms = [1,1,2]    # Number of observations of each strategy
qs = [0,2,1]    # Strictness thresholds: q=0 corresponds to ALLC, q>M corresponds to ALLD, 0<q<=M corresponds to DISC
# Paths for results
path = "demo/$norm/"
# Obtain attractors
@time steady_states_detIC(path, parameters, norm, Ms, qs)
