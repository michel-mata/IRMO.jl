# IRMO.jl

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://michel-mata.github.io/IRMO.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://michel-mata.github.io/IRMO.jl/dev)

Indirect Reciprocity model under private monitoring that aggregates multiple observations of a recipient to assign a reputation.

## System requirements
The code was produced on Julia v1.8.5 and uses the following packages:
```
DifferentialEquations v7.7.0
ForwardDiff v0.10.35
NLsolve v4.5.1
SciMLNLSolve v0.1.4
StaticArrays v1.5.19
StatsBase v0.33.21
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
For a demo, run the following code:
```
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
path = "demo/$norm"
# Obtain attractors
steady_states_detIC(path, parameters, norm, Ms, qs)
```
Or the file `demo.jl`.

The expected output is a folder `demo/SJ` containing the files:
`f0.csv`: containing the 171 initial frequencies of the simulation.
`ss.csv`: containing the reached steady states after integration of such initial frequencies.
`coop.csv`: containing the cooperation rate 

## Instructions for use and reproduction of data
For the evolutionary dynamics between $ALLC$, $ALLD$, and $DISC_{q,M}$ run: `interaction.jl`.

For the evolution of tolerance in the aggregation process run: `assignment.jl`.

For the evolution of the number of observations used to assign reputations run: `obsevations.jl`.

All data used to produce figures of the paper is available in the `results/` folder.

Copyright (c) 2023 Sebastian Michel-Mata