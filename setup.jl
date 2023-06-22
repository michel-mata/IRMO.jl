begin
    "setting directory..." |> print
    (@__DIR__) == pwd() || cd(@__DIR__)
    any(LOAD_PATH .== pwd()) || push!(LOAD_PATH, pwd())
    "done, set in \n\t'$(pwd())'!" |> println
    " "|>println

    "activating environment..." |> println
    using Pkg
    Pkg.activate(".")
    Pkg.instantiate()
    Pkg.resolve()
    "done!" |> println
    " "|>println
end

begin
    using Distributed
    "initial workers: $(nprocs())" |> println
    num_workers = @isdefined(num_workers) ? num_workers : 0
end

begin
    "loading module in Main... " |> print
    using IRMO
    "done!" |> println
end

begin
    "loading workers:" |> println
    "\t"|> print
    addprocs( num_workers-nprocs() , exeflags="--project=$(Base.active_project())" )
    @everywhere using IRMO
    "done, $(nprocs()) workers loaded!" |> println
end
