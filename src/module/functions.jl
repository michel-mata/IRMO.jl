"""
"""
function get_norm(
    norm_ID::String
    )
    # Select norm
    (norm_ID == "SJ") && (norm = [1 0; 0 1])
    (norm_ID == "SH") && (norm = [0 0; 0 1])
    (norm_ID == "SC") && (norm = [0 0; 1 1])
    (norm_ID == "SS") && (norm = [1 0; 1 1])

    return norm
end

"""
"""
_sample(n) = [0,sort(rand(n-1))...,1]
_diff(x) = circshift(x,-1)[1:end-1]-x[1:end-1]
get_sample_simplex(n) = _diff(_sample(n))

function get_initial_frequencies(S;random=false,grid=0.05:0.05:0.95) 
    if random
        freqs = round.(get_sample_simplex(S);digits=3)#[1:end-1]
        while (sum(freqs) !== 1.0) || (0.0 in freqs)
            freqs = round.(get_sample_simplex(S);digits=3)#[1:end-1]
        end
    else
        freqs = [ round.([fq...,1-sum(fq)];digits=3) for fq in Base.product(repeat([grid],S-1)...)  if sum(fq)<1 ]
    end

    return freqs
end

function chop_pt(ss)
    # ss = [ss...,1-sum(ss)]
    ss[abs.(ss) .< 1e-4] .= 0
    ss ./= sum(ss)
    return ss
end

function get_vertex(S)
    fqs=[]
    for i in 1:S
        zs = zeros(S)
        zs[i]=1 
        push!(fqs,zs)
    end
    return fqs
end

function get_eigenvalues(path,norm,Ms,qs,parameters)
    S = length(Ms)
    reps = reshape(repeat([1.0],S*S),(S,S))
    dr, df = get_full_dynamics(norm,Ms,qs,parameters)
    ss = readdlm(path*"ss-f.csv",',')
    ev = []
    for f0 in eachrow(ss)
        r0 = dr(f0,reps)
        _df(f0) = df(f0,r0)
        push!(ev,eigvals(ForwardDiff.jacobian(_df,f0)))
    end
    ev = transpose(hcat(ev...))
    writedlm(path*"ev.csv",ev,',')
    return ev
end

"""
"""
H(p,q,M) = sum( binomial(big(M),m) * p^m * (1-p)^(M-m) for m in q:M ; init=0)

"""
"""
function get_full_dynamics(norm,Ms,qs,parameters)
    # Game and errors
    b, c, α, ϵ = parameters
    # Norm
    norm_d = get_norm(norm)[1,:]
    norm_c = get_norm(norm)[2,:]
    # Reputation dynamics
    H(p,q,M) = sum( binomial(big(M),m) * p^m * (1-p)^(M-m) for m in q:M ; init=0)
    obs(r,i,j,k) = (1-ϵ)*r[j,k]*(norm_c⋅[1-r[i,k],r[i,k]]) .+ (1-((1-ϵ)*r[j,k]))*(norm_d⋅[1-r[i,k],r[i,k]])
    rep(f,r,q,M,i,j) = H( α+(1-2α)*sum( f[k]*obs(r,i,j,k) for k in eachindex(f) ), q, M )
    reputation_dynamics(dr,r,f,t) = dr[1:length(f),1:length(f)] = [ rep(f,r,qs[i],Ms[i],i,j)-r[i,j] for i in eachindex(f), j in eachindex(f) ]
    reputation(f,r) = solve(ODEProblem(reputation_dynamics,r,(0.,1e+3),f), reltol=1e-8, abstol=1e-8, saveat=1e+3).u |> last
    
    payoff(f,r) = [ sum( (1-ϵ)*f[j]*( b*r[j,i] - c*r[i,j] ) for j in eachindex(f) ) for i in eachindex(f) ]
    replicator(f,r) = f .* ( payoff(f,r) .- sum(f.*payoff(f,r)) )
    return reputation, replicator
end

"""
"""
function get_rep(norm,Ms,qs,parameters)
    reputation,_= get_full_dynamics(norm,Ms,qs,parameters)
    return reputation
end

"""
"""
function get_time_dynamics(norm,Ms,qs,parameters)

    reputation,replicator = get_full_dynamics(norm,Ms,qs,parameters)
    replicator_dynamics(df,f,r,t) = df[1:length(f)] = replicator(f,reputation(f,r))#replicator([f...,1-sum(f)],reputation([f...,1-sum(f)],r))[1:end-1]

    return replicator_dynamics
end

"""
"""
function steady_states_detIC(path, parameters, norm, Ms, qs; r0=1.0, refine=true)
    # Paths for results
    !ispath(path) && mkpath(path)
    path_res = path*"ss.csv"
    path_f = path*"f0.csv"
    path_c = path*"coop.csv"

    # Settings    
    S = length(Ms)
    reps = reshape(repeat([r0],S*S),(S,S))
    # Start or load
    if !isfile(path_res)
        # Sample initial frequencies
        freqs = get_initial_frequencies(S)
    elseif refine
        # Refine
        @warn("Results exist, refining!")
        refine_steady_states(path, parameters, norm, Ms, qs)
        return
    else
        # Skip
        @error("Results exist, not overwriting!")
        return 
    end
    # Dynamics
    reputation,replicator = get_full_dynamics(norm,Ms,qs,parameters)
    # Storage
    steady_states = SharedArray{Float64,2}(length(freqs),S)
    cooperation = SharedArray{Float64,2}(length(freqs),S)
    # Progress
    done = SharedArray{Float64,1}(length(freqs))
    # Find attractors
    @sync @distributed for i in eachindex(freqs)
        "starting freq $i: $(round.(freqs[i];digits=3))" |> println
        # Get dynamics
        replicator_dynamics = get_time_dynamics(norm,Ms,qs,parameters)
        # Transient
        freq = copy(freqs[i])
        df = replicator(freq,reputation(freq,reps))
        T = 1e0
        while !all( abs.(df) .< 1e-5)
            f0 = chop_pt(last(solve(ODEProblem(replicator_dynamics,freq,(0.0,T),reps),saveat=T,RadauIIA3(),abstol=1e-15).u))
            if all( 0 .<= f0 .<= 1)
                freq = f0
                df = replicator(freq,reputation(freq,reps))
                all(abs.(df) .< 5e-2) && (T < 1e1) && (T *= 10)
                all(abs.(df) .< 1e-3) && (T < 1e2) && (T *= 10)
                all(abs.(df) .< 1e-4) && (T < 1e3) && (T *= 10)
                all(abs.(df) .< 5e-5) && (T < 1e4) && (T *= 10)
            else
                T /= 10
            end
            df = replicator(freq,reputation(freq,reps))
            "$T - $(round.(df;digits=8))\t$(round.(freq;digits=4))" |> println
        end
        rep = get_rep(norm,Ms,qs,parameters)(freq,reps)
        freq = chop_pt(freq)
        coop = freq.*(rep*freq)
        # Save
        steady_states[i,:] .= freq
        cooperation[i,:] .= coop
        # Record progress and report
        done[i] += 1
        "steady state: \t\t $(round.(steady_states[i,:];digits=3))\t\t $(round(100sum(done)/length(freqs);digits=3))%" |> println
    end
    
    # Save
    writedlm(path_f,freqs,',')
    writedlm(path_res,steady_states,',')
    writedlm(path_c, cooperation, ',')
end

"""
"""
delete_zeros(fs) = hcat([row for row in eachrow(fs) if !all(row .== 0)]...)'|>Matrix

"""
"""
function steady_states_rndIC(path, parameters, norm, Ms, qs, num_ic; r0=1.0, refine=false)
    # Paths for results
    !ispath(path) && mkpath(path)
    path_res = path*"ss.csv"
    path_f = path*"f0.csv"
    path_c = path*"coop.csv"

    # Settings    
    S = length(Ms)
    reps = reshape(repeat([r0],S*S),(S,S))

    # Storage
    frequencies = SharedArray{Float64,2}(num_ic,S)
    steady_states = SharedArray{Float64,2}(num_ic,S)
    cooperation = SharedArray{Float64,2}(num_ic,S)
    # Initialize or expand
    if !isfile(path_res)
        # Sarting 
        n0 = 0
    else
        f0 = readdlm(path_f,',')|>delete_zeros
        ss = readdlm(path_res,',')|>delete_zeros
        coop = readdlm(path_c,',')|>delete_zeros
        @show size(f0)
        # starting
        n0 = size(f0,1)
        # Load previous
        if n0 < num_ic
            frequencies[1:n0,:] = f0
            steady_states[1:n0,:] = ss
            cooperation[1:n0,:] = coop
        elseif refine
            # Refine
            @warn("Results exist, refining!")
            refine_steady_states(path, parameters, norm, Ms, qs)
            return
        else
            # Skip
            @error("Larger sample size exists: (current - $n0, arg - $num_ic )")
            return 
        end
    end
    # Dynamics
    reputation,replicator = get_full_dynamics(norm,Ms,qs,parameters)
    # Progress
    done = SharedArray{Float64,1}(num_ic)
    done[1:n0,:] .+= 1
    n0 += 1
    # Find attractors
    @sync @distributed for i in n0:num_ic
        # Get initial frequency
        f0 = get_initial_frequencies(S, random=true)
        while vec(f0) ∈ eachrow(frequencies)
            @warn("Repeated initial frequency: $f0")
            f0 = get_initial_frequencies(S, random=true)
        end
        freq = f0
        # Integrate
        "starting freq $(round.(freq;digits=3))" |> println
        # Get dynamics
        replicator_dynamics = get_time_dynamics(norm,Ms,qs,parameters)
        # Solve fine-grain
        df = replicator(freq,reputation(freq,reps))
        # Integrate
        T = 1e0
        while !all( abs.(df) .< 1e-5)
            tr = chop_pt(last(solve(ODEProblem(replicator_dynamics,freq,(0.0,T),reps),saveat=T,RadauIIA5(),abstol=1e-10).u))
            if all( 0 .<= tr .<= 1)
                freq = tr
                df = replicator(freq,reputation(freq,reps))
                all(abs.(df) .< 5e-2) && (T < 1e1) && (T *= 10)
                all(abs.(df) .< 5e-3) && (T < 1e2) && (T *= 10)
                all(abs.(df) .< 5e-4) && (T < 1e3) && (T *= 10)
                all(abs.(df) .< 5e-5) && (T < 1e4) && (T *= 10)
            else
                T /= 10
            end
            df = replicator(freq,reputation(freq,reps))
            "$T - $(round.(df;digits=8))\t$(round.(freq;digits=4))" |> println
        end
        rep = get_rep(norm,Ms,qs,parameters)(freq,reps)
        freq = chop_pt(freq)
        coop = freq.*(rep*freq)
        # Save
        frequencies[i,:] .= f0
        steady_states[i,:] .= freq
        cooperation[i,:] .= coop
        writedlm(path_f, frequencies, ',')
        writedlm(path_res, steady_states, ',')
        writedlm(path_c, cooperation, ',')
        # Record progress and report
        done[i] += 1
        "steady state: \t\t\t $(round.(steady_states[i,:];digits=3))\t $(round(100sum(done)/num_ic;digits=3))%" |> println
    end

    # Save
    writedlm(path_f, frequencies, ',')
    writedlm(path_res, steady_states, ',')
    writedlm(path_c, cooperation, ',')
end

"""
"""
function trajectories_detIC(path, parameters, norm, Ms, qs; T=1e1, r0=1.0, expand=false)
    # Paths for results
    !ispath(path) && mkpath(path)
    path_f = path*"f0.csv"
    # Settings    
    S = length(Ms)
    reps = reshape(repeat([r0],S*S),(S,S))
    # Sample initial frequencies
    freqs = get_initial_frequencies(S)
    # Progress
    done = SharedArray{Float64,1}(length(freqs))
    # Find attractors
    @sync @distributed for i in eachindex(freqs)
        "starting freq $i" |> println
        # Path
        path_res = path*"$i.csv"
        # Start or load
        if !isfile(path_res)
            freq = freqs[i]
        elseif expand
            freq = readdlm(path_res,',')[end,:]
        else
            @error("Results exist, not overwriting!")
            continue 
        end
        # Dynamics
        reputation,replicator = get_full_dynamics(norm,Ms,qs,parameters)
        # Time dynamics
        replicator_dynamics = get_time_dynamics(norm,Ms,qs,parameters)
        trajectory = []
        # Transient
        @time trj = chop_pt.(solve(ODEProblem(replicator_dynamics,freq,(0.0,1.0),reps),abstol=1e-15,RadauIIA5()).u)
        freq = last(trj)
        trj = transpose(hcat(trj...))
        push!(trajectory,trj)
        df = replicator(freq,reputation(freq,reps))
        # Integrate
        @time while !all(df .< 1e-8)
            @time trj = chop_pt.(solve(ODEProblem(replicator_dynamics,freq,(0.0,T),reps),abstol=1e-15,RadauIIA5()).u)
            if all([ all( 0 .<= t .<= 1) for t in trj ])
                freq = last(trj)
                trj = transpose(hcat(trj...))
                push!(trajectory,trj)
                all(abs.(df) .< 5e-2) && (T < 1e1) && (T *= 10)
                all(abs.(df) .< 5e-3) && (T < 1e2) && (T *= 10)
                all(abs.(df) .< 5e-4) && (T < 1e3) && (T *= 10)
                all(abs.(df) .< 5e-5) && (T < 1e4) && (T *= 10)
                all(abs.(df) .< 5e-6) && (T < 1e5) && (T *= 10)
                all(abs.(df) .< 5e-7) && (T < 1e6) && (T *= 10)
            else
                T /= 10
            end
            df = replicator(freq,reputation(freq,reps))
            "$T - $(round.(df;digits=12))\t$(round.(freq;digits=4))" |> println
        end
        if length(trajectory)>1
            trajectory = transpose(hcat(unique(eachrow(vcat(trajectory...)))...))
        else
            trajectory = first(trajectory)
        end
        
        # Save steps
        if !isfile(path_res)
            # Save
            writedlm(path_res,trajectory,',')
        else
            # Add new steps
            f0 = readdlm(path_res,',')
            writedlm(path_res, transpose(hcat(unique(eachrow(vcat(f0,trajectory)))...)), ',')
        end
        # Record progress and report
        done[i] += 1
        "final frequency: $(round.(last(eachrow(trj));digits=3))\t\t $(round(100sum(done)/length(freqs);digits=3))%" |> println
    end
    # Save
    writedlm(path_f,freqs,',')
end

"""
"""
function refine_steady_states(path, parameters, norm, Ms, qs; r0=1.0)
    # Paths for results
    path_res = path*"ss.csv"
    path_ref = path*"ss-f.csv"
    # Settings    
    S = length(Ms)
    reps = reshape(repeat([r0],S*S),(S,S))
    # Start or load
    if !isfile(path_res)
        # Sample initial frequencies
        @error("Results do not exist, nothing to refine!")
        return 
    elseif isfile(path_ref)
        # Read results
        freqs = readdlm(path_res,',')|>eachrow|>collect
        ref_freqs = readdlm(path_ref,',')
        # Count refined versus results
        n0 = length(freqs)
        n1 = size(ref_freqs|>delete_zeros,1)
        # If not all refined
        if n0 == n1
            # Sample initial frequencies
            @error("Results already refined, not overwriting!")
            return
        elseif n1 < n0
            # Storage
            steady_states = SharedArray{Float64,2}(n0,S)
            steady_states .= ref_freqs
            inx = findall([ fq == zeros(S) for fq in eachrow(ref_freqs) ])
        else
            @error("Something wrong with files: results $n0 and refined $n1 !")
            return
        end
    else
        # Load results
        freqs = readdlm(path_res,',')|>eachrow|>collect
        n0 = length(freqs)
        # Storage
        steady_states = SharedArray{Float64,2}(n0,S)
        inx = collect(1:n0)
    end
    # Progress
    done = SharedArray{Float64,1}(n0)
    done .+= 1
    done[inx] .= 0
    # Find attractors
    @sync @distributed for i in inx
        "starting freq $i: $(round.(freqs[i];digits=3))" |> println
        # Get dynamics
        replicator_dynamics = get_time_dynamics(norm,Ms,qs,parameters)
        # Solve
        T = 1e0
        freq = copy(freqs[i])
        reputation, replicator = get_full_dynamics(norm,Ms,qs,parameters)
        df = replicator(freq,reputation(freq,reps))
        while !all( abs.(df) .< 1e-10)
            tr = chop_pt(last(solve(ODEProblem(replicator_dynamics,freq,(0.0,T),reps),saveat=T,RadauIIA5(),abstol=1e-15).u))
            if all( 0 .<= tr .<= 1)
                freq = tr
                df = replicator(freq,reputation(freq,reps))
                # all(abs.(df) .< 5e-5) && (T < 1e3) && (T *= 10)
                all(abs.(df) .< 5e-5) && (T < 1e6) && (T *= 10)
                all(abs.(df) .< 5e-7) && (T < 1e8) && (T *= 10)
            else
                T /= 10
            end
            df = replicator(freq,reputation(freq,reps))
            "$T - $(round.(df;digits=12))\t$(round.(freq;digits=4))" |> println
        end
        prob = NonlinearProblem(ODEProblem(replicator_dynamics,freq,(0.0,1e1),reps))
        freq = chop_pt(solve(prob,DynamicSS(RadauIIA5()),abstol=1e-10))
        # Save
        steady_states[i,:] .= chop_pt(freq)
        writedlm(path_ref,steady_states,',')
        # Record progress and report
        done[i] += 1
        "\nsteady state: \t\t $(round.(steady_states[i,:];digits=3))\t\t $(round(100sum(done)/length(freqs);digits=3))%\n" |> println
    end    
    # Save
    writedlm(path_ref,steady_states,',')
end
