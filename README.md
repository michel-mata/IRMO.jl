# IRMO.jl

[![Build status (Github Actions)](https://github.com/michel-mata/IRMO.jl/workflows/CI/badge.svg)](https://github.com/michel-mata/IRMO.jl/actions)
[![Documentation](https://github.com/michel-mata/IRMO.jl/actions/workflows/Documentation.yml/badge.svg)](https://github.com/michel-mata/IRMO.jl/actions/workflows/Documentation.yml)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://michel-mata.github.io/IRMO.jl)

Indirect Reciprocity model under private monitoring that aggregates multiple observations of a recipient to assign a reputation.


## System requirements
The code was produced on Julia v1.8.5 and uses the following packages:
```
DifferentialEquations
ForwardDiff
NLsolve
SciMLNLSolve
StaticArrays
StatsBase
DelimitedFiles
Distributed
LinearAlgebra
SharedArrays
```
For more information, check `Project.toml` and `Manifest.toml`.

## Installation guide
For installation of all packages and precompilation of the `IRMO.jl` module, run: `setup.jl`.

Expected installation and precompilation time is 5 minutes.

## Demo
For a demo, run the file `demo.jl` or the following code:
```julia
# Load packages
include("./setup.jl")
# Game parameters
parameters = [1.0,0.2,0.02,0.02]
# Social norm
norm = "SJ"
# Sets of strategies
Ms = [1,1,1]
qs = [0,2,1]
# Paths for results
path = "demo/$norm/"
# Obtain attractors
steady_states_detIC(path, parameters, norm, Ms, qs)
```

The expected output is a folder `demo/SJ` containing the files:
- `f0.csv`: containing the 171 initial frequencies of the simulation.
- `ss.csv`: containing the reached steady states after integration of such initial frequencies.
- `coop.csv`: containing the cooperation level of a population in such steady state.

Expected running time of demo with a single worker is 15 minutes.

## Instructions for use and reproduction of data
- For the evolutionary dynamics between $ALLC$, $ALLD$, and $DISC_{q,M}$ run: `interaction.jl`.
- For the time trajectories of the evolutionary dynamics between $ALLC$, $ALLD$, and $DISC_{q,M}$ run: `trajectories.jl`.
- For the evolution of elements of the aggregation rule ($q$ or $M$) run: `robustness.jl`.
- For the coevolution of the aggregation rule ($q$,$M$) run: `coevolution.jl`.
- For the evolution of social norms or strategy competition under heterogeneous norms run: `norms.jl`.
- For the evolutionary dynamics between $DISC_{q,M}$ and $probC$ run: `probabilistic.jl`.

All data used to produce figures of the paper is available in the `results/` folder.


---
Copyright (c) 2023 Sebastian Michel-Mata