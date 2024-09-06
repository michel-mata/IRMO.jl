"set up!" |> println
include("./setup.jl")
using DelimitedFiles

parameters = [1.0,0.2,0.05,0.05]
ms = 2:2:10
freq = [1,0]
reps = reshape(repeat([1.0],2*2),(2,2))

for norm in social_norms[1:2], r in ms
    path = "results/reputations/$norm/"
    !ispath(path) && mkpath(path)
    rs = hcat([ get_rep(norm,[r,m],[r==1 ? r : r÷2,m==1 ? m : m÷2],parameters)(freq,reps)[:] for m in ms ]...)
    writedlm(path*"$r.csv",rs,',')
end

errors = [0.01:0.01:0.1...]
for norm in social_norms[1:2], m in ms
    rps_min = zeros(length(errors),4)
    rps_max = zeros(length(errors),4)
    for e in eachindex(errors)
        ϵ = errors[e]
        # Min
        Ms = [m,m+2]
        qs = [m÷2,m÷2+1]
        rps_min[e,:] .= get_rep(norm,Ms,qs,[1.,.2,ϵ,ϵ])([1,0],reps)[:]
        # Max
        Ms = [m+2,m]
        qs = [m÷2+1,m÷2]
        rps_max[e,:] .= get_rep(norm,Ms,qs,[1.,.2,ϵ,ϵ])([1,0],reps)[:]
    end
    lns_min = (x-> isnan(x) ? 0.0 : x).([ (rps_min[i,1].-rps_min[i,2])./(rps_min[i,1].-rps_min[i,3]) for i in eachindex(errors) ])
    lns_max = (x-> isnan(x) ? 0.0 : x).([ (rps_max[i,1].-rps_max[i,2])./(rps_max[i,1].-rps_max[i,3]) for i in eachindex(errors) ])
    
    path = "results/reputations/$norm/$m-"
    !ispath(path) && mkpath(path)
    writedlm(path*"min.csv",lns_min,',')
    writedlm(path*"max.csv",lns_max,',')
end