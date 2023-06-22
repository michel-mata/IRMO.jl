"set up!" |> println
num_workers = 9
include("./setup.jl")

# Game parameters
parameters = [1.0,0.2,0.02,0.02]

# For all social norms
for norm in social_norms
    
    # Baseline
    Ms = [1,1,1]
    qs = [0,2,1]
    # Path for trajectories
    path = "results/trajectories/baseline/$norm/"
    "\n\n"*path*"\n" |> println
    # Obtain trajectories
    trajectories_detIC(path, parameters, norm, Ms, qs)

    # Trios with unconditional strategies
    for M in 2:2:10
        thresholds = [1,Int(M/2),M]
        ql = ["min","med","max"]

        for q in eachindex(thresholds)
            Ms = [1,1,M]
            qs = [0,2,thresholds[q]]
            # Paths for results
            path = "results/trajectories/unconditional/$norm/M$M/$(ql[q])/"
            "\n\n"*path*"\n" |> println
            # Obtain trajectories
            trajectories_detIC(path, parameters, norm, Ms, qs)
        end
    end

    # Look twice (competition of tolerance)
    begin
        Ms = [2,2,2]
        thresholds = [[2,1,0],[1,2,3],[3,0,1],[0,3,2]]
        # Paths for results
        path = "results/trajectories/assignment/$norm/Qs/"
        "\n\n"*path*"\n" |> println
        # Obtain trajectories
        for q in eachindex(thresholds)
            trajectories_detIC(path*"$q/", parameters, norm, Ms, thresholds[q])
        end
    end

    # Look twice forgive once (competition of observations)
    begin
        # Sets of strategies
        Ms = [[1,2,2],[2,1,2],[2,2,2],[2,2,1]]
        thresholds = [[1,1,0],[1,1,3],[3,0,1],[0,3,1]]
        # Paths for results
        path = "results/trajectories/observations/$norm/Ms/"
        "\n\n"*path*"\n" |> println
        # Obtain trajectories
        for q in eachindex(thresholds)
            trajectories_detIC(path*"$q/", parameters, norm, Ms[q], thresholds[q])
        end
    end
end