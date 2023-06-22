"set up!" |> println
include("./setup.jl")
begin
    using DifferentialEquations
    using SciMLNLSolve
    using NLsolve
    using LinearAlgebra
    using DelimitedFiles
    using Distributed
    using SharedArrays
    using ForwardDiff
end

for norm in social_norms, M in [2,4,6,8,10]
    for q in 1:3, uncond in [0,2]
        path = "results/stability/$norm/uncond/$M/"
        path *= q == 1 ? "min-" : q==3 ? "max-" : "med-"
        path *= uncond == 0 ? "ALLC/" : "ALLD/"
        "\n\n"*path*"\n" |> println
        mkpath(path)
        Ms = [1,M]
        qs = [uncond,[1,Int(M/2),M][q]]
        evals = []
        for c in 0.0:0.1:1.0, α in 0.0:0.02:0.3, ϵ in 0.0:0.02:0.3
            parameters = [1.0,c,α,ϵ]
            S = length(Ms)
            reps = reshape(repeat([1.0],S*S),(S,S))
            dr, df = get_full_dynamics(norm,Ms,qs,parameters)
            f0 = [0,1]
            r0 = dr(f0,reps)
            _df(frq) = df(frq,r0)
            ev = eigvals(ForwardDiff.jacobian(_df,f0))
            stable = all(sign.(ev) .< 0)
            "($M,$q,$uncond,$stable)\t$(round.(ev,digits=10))"|>println
            push!(evals,[c,α,ϵ,ev...])
        end
        writedlm(path*"evals.csv",evals,',')
    end
end

for norm in social_norms, M in [2,4,6,8,10]
    for q in 1:3
        path = "results/stability/$norm/uncond-all/$M/"
        path *= q == 1 ? "min/" : q==3 ? "max/" : "med/"
        "\n\n"*path*"\n" |> println
        mkpath(path)
        Ms = [1,1,M]
        qs = [0,2,[1,Int(M/2),M][q]]
        evals = []
        for b in 1:20, α in 0.0:0.02:0.5, ϵ in 0.0:0.02:0.5
            parameters = [b,1.0,α,ϵ]
            S = length(Ms)
            reps = reshape(repeat([1.0],S*S),(S,S))
            dr, df = get_full_dynamics(norm,Ms,qs,parameters)
            f0 = [0,0,1]
            r0 = dr(f0,reps)
            _df(frq) = df(frq,r0)
            ev = eigvals(ForwardDiff.jacobian(_df,f0))
            stable = all(sign.(ev) .< 0)
            "($M,$q,$stable)\t$(round.(ev,digits=10))"|>println
            push!(evals,[b,α,ϵ,ev...])
        end
        writedlm(path*"evals.csv",evals,',')
    end
end