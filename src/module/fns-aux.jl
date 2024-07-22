"""
Get the social norm matrix
"""
function get_norm(
    norm_ID::String
    )
    # Select norm
    (norm_ID == "SJ") && (norm = [1 0; 0 1])
    (norm_ID == "SH") && (norm = [0 0; 0 1])
    (norm_ID == "SC") && (norm = [0 0; 1 1])
    (norm_ID == "SS") && (norm = [1 0; 1 1])
    # For norm competition
    (norm_ID == "ALLC") && (norm = [1 1; 1 1])
    (norm_ID == "ALLD") && (norm = [0 0; 0 0])

    return norm
end

"""
Functions for sampling from the simplex
"""
_sample(n) = [0,sort(rand(n-1))...,1]
_diff(x) = circshift(x,-1)[1:end-1]-x[1:end-1]
get_sample_simplex(n) = _diff(_sample(n))

function get_initial_frequencies(S;random=false,grid=0.05:0.05:0.95) 
    if random
        freqs = round.(get_sample_simplex(S);digits=3)
        while (sum(freqs) !== 1.0) || (0.0 in freqs)
            freqs = round.(get_sample_simplex(S);digits=3)
        end
    else
        freqs = [ round.([fq...,1-sum(fq)];digits=3) for fq in Base.product(repeat([grid],S-1)...)  if sum(fq)<1 ]
    end

    return freqs
end

"""
Remove numerical zeros
"""
function chop_pt(ss)
    ss[abs.(ss) .< 1e-4] .= 0
    ss ./= sum(ss)
    return ss
end

"""
Remove rows of zeros
"""
delete_zeros(fs) = hcat([row for row in eachrow(fs) if !all(row .== 0)]...)'|>Matrix

"""
CDF of the binomial distribution
"""
H(p,q,M) = sum( binomial(big(M),m) * p^m * (1-p)^(M-m) for m in q:M ; init=0)