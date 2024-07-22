
"""
Reputation and replicator dynamics for strategy competition under fixed mixture of norms
"""
function get_full_dynamics_nf(norms,Ms,qs,parameters,n)
    # Game and errors
    b, c, α, ϵ = parameters
    # Expand parameters by norms of size n
    nrms = repeat(norms,length(Ms))
    Ms = vcat(fill.(Ms, length(n))...)
    qs = vcat(fill.(qs, length(n))...)
    # Norm
    norm_d = [get_norm(nrm)[1,:] for nrm in nrms]
    norm_c = [get_norm(nrm)[2,:] for nrm in nrms]
    # Reputation dynamics
    H(p,q,M) = sum( binomial(big(M),m) * p^m * (1-p)^(M-m) for m in q:M ; init=0)
    obs(r,i,j,k) = (1-ϵ)*r[j,k]*(norm_c[i]⋅[1-r[i,k],r[i,k]]) .+ (1-((1-ϵ)*r[j,k]))*(norm_d[i]⋅[1-r[i,k],r[i,k]])
    rep(f,r,q,M,i,j) = H( α+(1-2α)*sum( f[k]*obs(r,i,j,k) for k in eachindex(f) ), q, M )
    reputation_dynamics(dr,r,f,t) = dr[1:length(f),1:length(f)] = [ rep(f,r,qs[i],Ms[i],i,j)-r[i,j] for i in eachindex(f), j in eachindex(f) ]
    solve_reputation(f,r) = solve(ODEProblem(reputation_dynamics,r,(0.,1e+3),f), reltol=1e-8, abstol=1e-8, saveat=1e+3).u |> last 
    function reputation(f,r) 
        freq = (n.*repeat(f,1,length(n))')[:]
        reps = solve_reputation(freq,r)
        ns = length(n)
        return [ n⋅(reps[i:i+ns-1,j:j+ns-1]*n) for i in 1:ns:length(freq), j in 1:ns:length(freq) ]
    end
    # Replicator dynamics
    payoff(f,r) = [ sum( (1-ϵ)*f[j]*( b*r[j,i] - c*r[i,j] ) for j in eachindex(f) ) for i in eachindex(f) ]
    replicator(f,r) = f .* ( payoff(f,r) .- sum(f.*payoff(f,r)) )

   return reputation, replicator
end

"""
Get reputation dynamics for strategy competition under fixed mixture of norms
"""
function get_rep_nf(norms,Ms,qs,parameters,n)
   reputation,_= get_full_dynamics_nf(norms,Ms,qs,parameters,n)
   return reputation
end

"""
Get replicator dynamics for time trajectories for strategy competition under fixed mixture of norms
"""
function get_time_dynamics_nf(norms,Ms,qs,parameters,n)

   reputation,replicator = get_full_dynamics_nf(norms,Ms,qs,parameters,n)
   replicator_dynamics(df,f,r,t) = df[1:length(f)] = replicator(f,reputation(f,r))

   return replicator_dynamics
end

"""
Integrate numerically until a steady state is reached.
The initial frequencies are deterministic, evenly spaced, and determined by the global variable `grid`.
Function fitted for strategy competition under fixed mixture of norms.
"""
function steady_states_detIC_nf(path, parameters, norms, Ms, qs, n; r0=1.0, refine=false)
   # Paths for results
   !ispath(path) && mkpath(path)
   path_res = path*"ss.csv"
   path_f = path*"f0.csv"
   path_c = path*"coop.csv"

   # Settings    
   S = length(Ms)
   reps = reshape(repeat([r0],S^2*length(n)^2),(S*length(n),S*length(n)))
   # Storage
   freqs = get_initial_frequencies(S)
   frequencies = SharedArray{Float64,2}(length(freqs),S)
   steady_states = SharedArray{Float64,2}(length(freqs),S)
   cooperation = SharedArray{Float64,2}(length(freqs),S)
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
       if n0 < length(freqs)
           frequencies[1:n0,:] = f0
           steady_states[1:n0,:] = ss
           cooperation[1:n0,:] = coop
       elseif refine
           # Refine
           @warn("Results exist, refining!")
           refine_steady_states_nf(path, parameters, norms, Ms, qs,n)
           return
       else
           # Skip
           @error("Larger sample size exists: (current - $n0, arg - $(length(freqs)) )")
           return 
       end
   end
   # Dynamics
   reputation,replicator = get_full_dynamics_nf(norms,Ms,qs,parameters,n)
   # Progress
   done = SharedArray{Float64,1}(length(freqs))
   done[1:n0,:] .+= 1
   n0 += 1
   # Find attractors
   @sync @distributed for i in n0:length(freqs)
       "starting freq $i: $(round.(freqs[i];digits=3))" |> println
       # Get dynamics
       replicator_dynamics = get_time_dynamics_nf(norms,Ms,qs,parameters,n)
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
               all(abs.(df) .< 5e-3) && (T < 1e2) && (T *= 10)
               all(abs.(df) .< 5e-4) && (T < 1e3) && (T *= 10)
               all(abs.(df) .< 5e-5) && (T < 1e4) && (T *= 10)
           else
               T /= 10
           end
           df = replicator(freq,reputation(freq,reps))
           "$T - $(round.(df;digits=8))\t$(round.(freq;digits=4))" |> println
       end
       rep = get_rep_nf(norms,Ms,qs,parameters,n)(freq,reps)
       freq = chop_pt(freq)
       coop = freq.*(rep*freq)
       # Save
       frequencies[i,:] .= freqs[i]
       steady_states[i,:] .= freq
       cooperation[i,:] .= coop
       writedlm(path_f, frequencies, ',')
       writedlm(path_res, steady_states, ',')
       writedlm(path_c, cooperation, ',')
       # Record progress and report
       done[i] += 1
       "steady state: \t\t $(round.(steady_states[i,:];digits=3))\t\t $(round(100sum(done)/length(freqs);digits=3))%" |> println
   end
   
   # Save
   writedlm(path_f,frequencies,',')
   writedlm(path_res,steady_states,',')
   writedlm(path_c, cooperation, ',')
end

"""
Expand integration to a lower tolerance for strategy competition under fixed mixture of norms.
"""
function refine_steady_states_nf(path, parameters, norms, Ms, qs, n; r0=1.0)
   # Paths for results
   path_res = path*"ss.csv"
   path_ref = path*"ss-f.csv"
   # Settings    
   S = length(Ms)
   reps = reshape(repeat([r0],S^2*length(n)^2),(S*length(n),S*length(n)))
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
       replicator_dynamics = get_time_dynamics_nf(norms,Ms,qs,parameters,n)
       # Solve
       T = 1e0
       freq = copy(freqs[i])
       reputation, replicator = get_full_dynamics_nf(norms,Ms,qs,parameters,n)
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

"""
Calculate the rate of cooperation for strategy competition under fixed mixture of norms
"""
function get_cooperation_nf(path, parameters, norms, Ms, qs, n; r0=1.0)
    # Paths for results
    path_res = path*"ss-f.csv"
    path_c = path*"coop-f.csv"

    # Read steady-states
    if !isfile(path_res)
        # Sarting 
        @error("no results yet!")
    else
        ss = readdlm(path_res,',')
    end

    S = length(Ms)
    reps = reshape(repeat([r0],S^2*length(n)^2),(S*length(n),S*length(n)))
    num_ic = size(ss,1)

    reputation_dynamics = get_rep_nf(norms,Ms,qs,parameters,n)
    cooperation = SharedArray{Float64,2}(num_ic,S)
    
    @sync @distributed for i in 1:num_ic
        freq = ss[i,:]
        rep = reputation_dynamics(freq,reps)
        coop = freq.*(rep*freq)
        cooperation[i,:] .= coop
    end

    writedlm(path_c, cooperation, ',')
end