"set up!" |> println
num_workers = 16
include("./setup.jl")

# Parameters to sweep
using StatsBase
_Ms = 1:10
_qm = [ [q/m for q in 1:m] for m in _Ms ]
_qs = [ k  for (k,v) in countmap(vcat(_qm...)) if v>1 ]|>sort
_bs = [1.1,1.5,2,5,10]
_es = [0.01,0.02,0.05,0.1,0.2]
_ms = 2:6

# Evolution of q
for norm in social_norms, M in _ms, b in _bs, ϵ in _es
    # Paths for results
    evolving = "robustness"
    path = "results/endpoints/$evolving/q/$norm/$b-$ϵ/$M/"
    "\n\n"*path*"\n" |> println
    # Select parameters
    Ms = [1,1,repeat([M],M)...]
    qs = [0,2,1:M...]
    # Obtain attractors
    num_ic = 100
    steady_states_rndIC(path, [b,1,ϵ,ϵ], norm, Ms, qs, num_ic, refine=true)
    get_cooperation(path, [b,1,ϵ,ϵ], norm, Ms, qs)
end

# Evolution of M
for norm in social_norms, q in eachindex(_qs), b in _bs, ϵ in _es
    # Paths for results
    evolving = "robustness"
    path = "results/endpoints/$evolving/M/$norm/$b-$ϵ/$q/"
    "\n\n"*path*"\n" |> println
    # Select parameters
    ms = _Ms[findall(_qs[q] .∈ _qm)]
    Ms = [1,1,ms...]
    qs = [0,2,Int.(_qs[q].*ms)...]
    # Obtain attractors
    num_ic = 100
    steady_states_rndIC(path, [b,1,ϵ,ϵ], norm, Ms, qs, num_ic, refine=true)
    get_cooperation(path, [b,1,ϵ,ϵ], norm, Ms, qs)
end