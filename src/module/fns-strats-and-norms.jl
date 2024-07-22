"""
Integrate numerically until a steady state is reached.
The initial frequencies are deterministic, evenly spaced, and determined by the global variable `grid`.
Function fitted for strategy and norm independent evolution.
"""
function steady_states_normInd(path, parameters, Ns, Ms, Qs, strategies; r0=1.0, refine=false)
    # Paths for results
    !ispath(path) && mkpath(path)
    path_res = path*"ss.csv"
    path_f = path*"f0.csv"
    path_c = path*"coop.csv"

    # Settings
    begin
        b, c, α, ϵ, β, η, w = parameters
        S = length(Ms)
        N = length(Ns)
        ns = 1:N
        ss = 1:S
        # Expand parameters by norms
        ms = hcat(repeat([Ms], N)...)
        qs = hcat(repeat([Qs], N)...)
        # Norm
        norm_d = [get_norm(n)[1,:] for n in Ns]
        norm_c = [get_norm(n)[2,:] for n in Ns]
        norm_d = reshape(repeat(norm_d,S),(N,S))' |> Matrix
        norm_c = reshape(repeat(norm_c,S),(N,S))' |> Matrix
        # Reputation dynamics
        H(p,q,M) = sum( binomial(big(M),m) * p^m * (1-p)^(M-m) for m in q:M ; init=0)
        R(r,i,j,ni,nj) = [1-r[i,j,ni,nj],r[i,j,ni,nj]]
        a(r,i,j,ni,nj) = (1-ϵ)*(strategies[i]⋅R(r,i,j,ni,nj))
        obs(r,i,j,k,ni,nj,nk) = [a(r,j,k,nj,nk),1-a(r,j,k,nj,nk)]⋅[norm_c[i,ni]⋅R(r,i,k,ni,nk), norm_d[i,ni]⋅R(r,i,k,ni,nk)]
        rep(f,r,q,M,i,j,ni,nj) = H( α+(1-2α)*sum( f[k,nk]*obs(r,i,j,k,ni,nj,nk) for k in ss, nk in ns ), q, M )
        reputation_dynamics(dr,r,f,t) = dr[ss,ss,ns,ns] = [ rep(f,r,qs[i,ni],ms[i,ni],i,j,ni,nj)-r[i,j,ni,nj] for i in ss, j in ss, ni in ns, nj in ns ]
        reputation(f,r) = solve(ODEProblem(reputation_dynamics,r,(0.,1e3),f), reltol=1e-10, abstol=1e-10, saveat=1e3).u |> last        
        # Replication dynamics
        payoff(f,r,i,ni) = sum( (1-ϵ)*f[j,nj]*( b*r[j,i,nj,ni] - c*r[i,j,ni,nj] ) for j in ss, nj in ns )
        payoff(f,r) = [payoff(f,r,i,ni) for i in ss, ni in ns]
        ϕ(f,r,i,j,s,n) = 1/(1+exp( β*(payoff(f,r,i,j)-payoff(f,r,s,n)) ))
        function dynamics(f,r,s,n,w)
            ds = sum( f[i,n]*f[s,j]*ϕ(f,r,i,n,s,j) - f[s,n]*f[i,j]*ϕ(f,r,s,n,i,j) for j in ns, i in ss if i!=s; init=0.0)
            dn = sum( f[i,n]*f[s,j]*ϕ(f,r,s,j,i,n) - f[s,n]*f[i,j]*ϕ(f,r,s,n,i,j) for j in ns, i in ss if j!=n; init=0.0)
            return w*ds + (1-w)*dn
        end
        dynamics(f,r) = [dynamics(f,r,i,ni,w) for i in ss, ni in ns]
        replicator_dynamics(df,f,r,t) = df[ss,ns] = [dynamics(f,r,s,n,w) for s in ss, n in ns]
    end

    # Parameters
    norm_prop = [η...,1-sum(η)]
    freqs = [ hcat([f for _ in ns].*norm_prop...) for f in get_initial_frequencies(S) ]
    # Storage
    frequencies = SharedArray{Float64,2}(length(freqs),S*N)
    steady_states = SharedArray{Float64,2}(length(freqs),S*N)
    cooperation = SharedArray{Float64,2}(length(freqs),S*N)
    if !isfile(path_res)
        # Sarting 
        n0 = 0
        inx = eachindex(freqs)
    else
        f0 = readdlm(path_f,',')
        sts = readdlm(path_res,',')
        coop = readdlm(path_c,',')
        # select the missing ICs
        inx = findall(sum(sign.(sts),dims=2)[:] .== 0)
        @show length(inx)
        n0 = size(f0,1)-length(inx)
        # Load previous
        if n0 < length(freqs)
            frequencies .= f0
            steady_states .= sts
            cooperation .= coop
        elseif refine
            # Refine
            @warn("Results exist, refining!")
            refine_normInd(path, parameters, Ns, Ms, Qs, strategies)
            return
        else
            # Skip
            @error("Larger sample size exists: (current - $n0, arg - $(length(freqs)) )")
            return 
        end
    end
    # Progress
    done = SharedArray{Float64,1}(length(freqs))
    done[1:n0,:] .+= 1
    n0 += 1
    # Find attractors
    @sync @distributed for i in inx
        # "starting freq $i" |> println
        # Initial values
        reps = reshape(repeat([r0],S^2*N^2),(S,S,N,N))
        freq = copy(freqs[i])
        df = dynamics(freq,reputation(freq,reps))
        # Integrate
        T = 1e0
        @time while !all( abs.(df) .< 1e-12)
            f0 = chop_pt(last(solve(ODEProblem(replicator_dynamics,freq,(0.0,T),reps),saveat=T,RadauIIA3(),abstol=1e-15).u))
            if all( 0 .<= f0 .<= 1)
                freq = f0
                reps = reputation(freq,reps)
                df = dynamics(freq,reps)
                all(abs.(df) .< 5e-2) && (T < 1e1) && (T *= 10)
                all(abs.(df) .< 5e-3) && (T < 1e2) && (T *= 10)
                all(abs.(df) .< 5e-4) && (T < 1e3) && (T *= 10)
                all(abs.(df) .< 5e-5) && (T < 1e4) && (T *= 10)
                all(abs.(df) .< 5e-6) && (T < 1e5) && (T *= 10)
                all(abs.(df) .< 5e-7) && (T < 1e6) && (T *= 10)
                all(abs.(df) .< 5e-8) && (T < 1e7) && (T *= 10)
                all(abs.(df) .< 5e-9) && (T < 1e8) && (T *= 10)
                all(abs.(df) .< 5e-10) && (T < 1e9) && (T *= 10)
                all(abs.(df) .< 5e-11) && (T < 1e10) && (T *= 10)
            else
                T /= 10
            end
            df = dynamics(freq,reputation(freq,reps))
            # "$T - $(round.(df;digits=13))\t$(round.(freq;digits=4))\n" |> println
        end
        r = reputation(freq,reps)
        coop = [ freq[i,ni]*sum(freq[j,nj]*a(r,i,j,ni,nj) for j in ss, nj in ns) for i in ss, ni in ns]

        # Save
        frequencies[i,:] .= freqs[i][:]
        steady_states[i,:] .= freq[:]
        cooperation[i,:] .= coop[:]
        writedlm(path_f, frequencies, ',')
        writedlm(path_res, steady_states, ',')
        writedlm(path_c, cooperation, ',')
        # Record progress and report
        done[i] += 1
        "$i)\t\t$(round.(freq;digits=3))\t\t $(round(100sum(done)/length(freqs);digits=3))%" |> println
    end
    
    # Save
    writedlm(path_f,frequencies,',')
    writedlm(path_res,steady_states,',')
    writedlm(path_c, cooperation, ',')
end

"""
Expand integration to a lower tolerance for strategy and norm independent evolution.
"""
function refine_normInd(path, parameters, Ns, Ms, Qs, strategies; r0=1.0)
    path_res = path*"ss.csv"
    path_f = path*"f0.csv"
    path_c = path*"coop.csv"

    # Settings
    begin
        b, c, α, ϵ, β, η, w = parameters
        S = length(Ms)
        N = length(Ns)
        ns = 1:N
        ss = 1:S
        # Expand parameters by norms
        ms = hcat(repeat([Ms], N)...)
        qs = hcat(repeat([Qs], N)...)
        # Norm
        norm_d = [get_norm(n)[1,:] for n in Ns]
        norm_c = [get_norm(n)[2,:] for n in Ns]
        norm_d = reshape(repeat(norm_d,S),(N,S))' |> Matrix
        norm_c = reshape(repeat(norm_c,S),(N,S))' |> Matrix
        # Reputation dynamics
        H(p,q,M) = sum( binomial(big(M),m) * p^m * (1-p)^(M-m) for m in q:M ; init=0)
        R(r,i,j,ni,nj) = [1-r[i,j,ni,nj],r[i,j,ni,nj]]
        a(r,i,j,ni,nj) = (1-ϵ)*(strategies[i]⋅R(r,i,j,ni,nj))
        obs(r,i,j,k,ni,nj,nk) = [a(r,j,k,nj,nk),1-a(r,j,k,nj,nk)]⋅[norm_c[i,ni]⋅R(r,i,k,ni,nk), norm_d[i,ni]⋅R(r,i,k,ni,nk)]
        rep(f,r,q,M,i,j,ni,nj) = H( α+(1-2α)*sum( f[k,nk]*obs(r,i,j,k,ni,nj,nk) for k in ss, nk in ns ), q, M )
        reputation_dynamics(dr,r,f,t) = dr[ss,ss,ns,ns] = [ rep(f,r,qs[i,ni],ms[i,ni],i,j,ni,nj)-r[i,j,ni,nj] for i in ss, j in ss, ni in ns, nj in ns ]
        reputation(f,r) = solve(ODEProblem(reputation_dynamics,r,(0.,1e3),f), reltol=1e-10, abstol=1e-10, saveat=1e3).u |> last        
        # Replication dynamics
        payoff(f,r,i,ni) = sum( (1-ϵ)*f[j,nj]*( b*r[j,i,nj,ni] - c*r[i,j,ni,nj] ) for j in ss, nj in ns )
        payoff(f,r) = [payoff(f,r,i,ni) for i in ss, ni in ns]
        ϕ(f,r,i,j,s,n) = 1/(1+exp( β*(payoff(f,r,i,j)-payoff(f,r,s,n)) ))
        function dynamics(f,r,s,n,w)
            ds = sum( f[i,n]*f[s,j]*ϕ(f,r,i,n,s,j) - f[s,n]*f[i,j]*ϕ(f,r,s,n,i,j) for j in ns, i in ss if i!=s; init=0.0)
            dn = sum( f[i,n]*f[s,j]*ϕ(f,r,s,j,i,n) - f[s,n]*f[i,j]*ϕ(f,r,s,n,i,j) for j in ns, i in ss if j!=n; init=0.0)
            return w*ds + (1-w)*dn
        end
        dynamics(f,r) = [dynamics(f,r,i,ni,w) for i in ss, ni in ns]
        replicator_dynamics(df,f,r,t) = df[ss,ns] = [dynamics(f,r,s,n,w) for s in ss, n in ns]
    end

    # Parameters
    norm_prop = [η...,1-sum(η)]
    freqs = [ hcat([f for _ in ns].*norm_prop...) for f in get_initial_frequencies(S) ]
    # Storage
    frequencies = SharedArray{Float64,2}(length(freqs),S*N)
    steady_states = SharedArray{Float64,2}(length(freqs),S*N)
    cooperation = SharedArray{Float64,2}(length(freqs),S*N)
    if !isfile(path_res)
        @error("Results do not exist!")
        return 
    else
        f0 = readdlm(path_f,',')
        sts = readdlm(path_res,',')
        coop = readdlm(path_c,',')
        # Refine when more than two strategies coexist
        inx1 = findall(sum(sign.(sts),dims=2)[:] .== 0)
        inx2 = findall(sum(sign.(sts),dims=2)[:] .> 2)
        inx = [inx1...,inx2...]


        frequencies .= f0
        steady_states .= sts
        cooperation .= coop

        if length(inx) == 0 
            @error("Already refined!")
            return
        end
    end
    # Progress
    done = SharedArray{Float64,1}(length(inx))
    forloop = enumerate(inx) |> collect
    @show forloop
    # Find attractors
    @sync @distributed for x in forloop
        (n0,i) = x
        # Initial values
        reps = reshape(repeat([r0],S^2*N^2),(S,S,N,N))
        freq = copy(freqs[i])
        df = dynamics(freq,reputation(freq,reps))
        # Integrate
        T = 1e0
        @time while !all( abs.(df) .< 1e-15)
            f0 = chop_pt(last(solve(ODEProblem(replicator_dynamics,freq,(0.0,T),reps),saveat=T,RadauIIA3(),abstol=1e-15).u))
            if all( 0 .<= f0 .<= 1)
                freq = f0
                reps = reputation(freq,reps)
                df = dynamics(freq,reps)
                all(abs.(df) .< 5e-2) && (T < 1e1) && (T *= 10)
                all(abs.(df) .< 5e-3) && (T < 1e2) && (T *= 10)
                all(abs.(df) .< 5e-4) && (T < 1e3) && (T *= 10)
                all(abs.(df) .< 5e-6) && (T < 1e4) && (T *= 10)

                all(abs.(df) .< 1e-7) && (T < 1e5) && (T *= 10)
                all(abs.(df) .< 1e-8) && (T < 1e6) && (T *= 10)
                all(abs.(df) .< 1e-10) && (T < 1e8) && (T *= 10)
                all(abs.(df) .< 1e-12) && (T < 1e10) && (T *= 10)
            else
                T /= 10
            end
            df = dynamics(freq,reputation(freq,reps))
            # "$T - $(round.(df;digits=16))\t$(round.(freq;digits=4))\n" |> println
        end
        r = reputation(freq,reps)
        coop = [ freq[i,ni]*sum(freq[j,nj]*a(r,i,j,ni,nj) for j in ss, nj in ns) for i in ss, ni in ns]

        # Save
        frequencies[i,:] .= freqs[i][:]
        steady_states[i,:] .= freq[:]
        cooperation[i,:] .= coop[:]
        writedlm(path_f, frequencies, ',')
        writedlm(path_res, steady_states, ',')
        writedlm(path_c, cooperation, ',')
        # Record progress and report
        done[n0] += 1
        "$i)\t\t$(round.(freq;digits=3))\t\t $(round(100sum(done)/length(inx);digits=3))%" |> println
    end
    
    # Save
    writedlm(path_f,frequencies,',')
    writedlm(path_res,steady_states,',')
    writedlm(path_c, cooperation, ',')
end

"""
Calculate the rate of cooperation for strategy and norm independent evolution.
"""
function get_cooperation_normInd(path, parameters, Ns, Ms, Qs, strategies; r0=1.0)
    path_res = path*"ss.csv"
    path_c = path*"coop-f.csv"

    # Settings
    begin
        b, c, α, ϵ, β, η, w = parameters
        S = length(Ms)
        N = length(Ns)
        ns = 1:N
        ss = 1:S
        # Expand parameters by norms
        ms = hcat(repeat([Ms], N)...)
        qs = hcat(repeat([Qs], N)...)
        # Norm
        norm_d = [get_norm(n)[1,:] for n in Ns]
        norm_c = [get_norm(n)[2,:] for n in Ns]
        norm_d = reshape(repeat(norm_d,S),(N,S))' |> Matrix
        norm_c = reshape(repeat(norm_c,S),(N,S))' |> Matrix
        # Reputation dynamics
        H(p,q,M) = sum( binomial(big(M),m) * p^m * (1-p)^(M-m) for m in q:M ; init=0)
        R(r,i,j,ni,nj) = [1-r[i,j,ni,nj],r[i,j,ni,nj]]
        a(r,i,j,ni,nj) = (1-ϵ)*(strategies[i]⋅R(r,i,j,ni,nj))
        obs(r,i,j,k,ni,nj,nk) = [a(r,j,k,nj,nk),1-a(r,j,k,nj,nk)]⋅[norm_c[i,ni]⋅R(r,i,k,ni,nk), norm_d[i,ni]⋅R(r,i,k,ni,nk)]
        rep(f,r,q,M,i,j,ni,nj) = H( α+(1-2α)*sum( f[k,nk]*obs(r,i,j,k,ni,nj,nk) for k in ss, nk in ns ), q, M )
        reputation_dynamics(dr,r,f,t) = dr[ss,ss,ns,ns] = [ rep(f,r,qs[i,ni],ms[i,ni],i,j,ni,nj)-r[i,j,ni,nj] for i in ss, j in ss, ni in ns, nj in ns ]
        reputation(f,r) = solve(ODEProblem(reputation_dynamics,r,(0.,1e3),f), reltol=1e-10, abstol=1e-10, saveat=1e3).u |> last        
        # Replication dynamics
        payoff(f,r,i,ni) = sum( (1-ϵ)*f[j,nj]*( b*r[j,i,nj,ni] - c*r[i,j,ni,nj] ) for j in ss, nj in ns )
        payoff(f,r) = [payoff(f,r,i,ni) for i in ss, ni in ns]
        ϕ(f,r,i,j,s,n) = 1/(1+exp( β*(payoff(f,r,i,j)-payoff(f,r,s,n)) ))
        function dynamics(f,r,s,n,w)
            ds = sum( f[i,n]*f[s,j]*ϕ(f,r,i,n,s,j) - f[s,n]*f[i,j]*ϕ(f,r,s,n,i,j) for j in ns, i in ss if i!=s; init=0.0)
            dn = sum( f[i,n]*f[s,j]*ϕ(f,r,s,j,i,n) - f[s,n]*f[i,j]*ϕ(f,r,s,n,i,j) for j in ns, i in ss if j!=n; init=0.0)
            return w*ds + (1-w)*dn
        end
        dynamics(f,r) = [dynamics(f,r,i,ni,w) for i in ss, ni in ns]
        replicator_dynamics(df,f,r,t) = df[ss,ns] = [dynamics(f,r,s,n,w) for s in ss, n in ns]
    end

    # Read steady-states
    if !isfile(path_res)
        # Sarting 
        @error("no results yet!")
    else
        fqs = readdlm(path_res,',')
        num_ic = size(fqs,1)
    end
    # Storage
    cooperation = SharedArray{Float64,1}(num_ic)
    reps = reshape(repeat([r0],S^2*N^2),(S,S,N,N))
    # Get cooperation
    @sync @distributed for i in 1:num_ic
        freq = reshape(fqs[i,:],(S,N))
        r = reputation(freq,reps)
        coop = [ freq[i,ni]*sum(freq[j,nj]*a(r,i,j,ni,nj) for j in ss, nj in ns) for i in ss, ni in ns]
        cooperation[i] = sum(coop)
    end
    # Save
    writedlm(path_c, cooperation, ',')
end