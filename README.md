# IRMO.jl

[![Build status (Github Actions)](https://github.com/michel-mata/IRMO.jl/workflows/CI/badge.svg)](https://github.com/michel-mata/IRMO.jl/actions)
[![Documentation](https://github.com/michel-mata/IRMO.jl/actions/workflows/Documentation.yml/badge.svg)](https://github.com/michel-mata/IRMO.jl/actions/workflows/Documentation.yml)
[![DOCS](https://img.shields.io/badge/docs-stable-blue.svg)](https://michel-mata.github.io/IRMO.jl)
[![DOI](https://zenodo.org/badge/657317785.svg)](https://zenodo.org/doi/10.5281/zenodo.12795781)

Indirect Reciprocity model under private monitoring that aggregates multiple observations of a recipient to assign a reputation.

Zenodo DOI: 10.5281/zenodo.12795782

## System requirements

The code was produced on Julia v1.8.5 and the last release is compatible up to Julia v1.10.4.

Module uses the following packages:

```console
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

## Julia installation guide

Download and install Julia by following these [instructions](https://julialang.org/downloads/), or by running:

> Linux and MacOS:
>
> ```console
> $ curl -fsSL https://install.julialang.org | sh
> ```
>
> Windows:
>
> ```console
> > winget install julia -s msstore
> ```

This will install the `juliaup` installation manager.
Make sure to be up to date by running:

> ```console
> $ juliaup update
> $ juliaup default release
> ```

To install different versions or explore more options run `juliaup --help`.

### Run scripts

Once installed, Julia will be available via the command line interface. Then, a script like `my_script.jl` can be run as:

```console
$ julia my_script.jl
```

## Module installation guide

### Download the repository

Clone this repository either by the `Download ZIP` option under the `Code` dropdown menu above, or by running:

```console
$ git clone https://github.com/michel-mata/IRMO.jl.git
```

For more information on how to clone a repository visit the [documentation](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository).

### Install packages

For installation of all packages and precompilation of the `IRMO.jl` module, run `setup.jl` as:

```console
$ julia setup.jl
```

Assuming the command line is in the directory of the repository. Alternatively, run the script using the relative path to the repository as:

```console
$ julia path/to/repo/setup.jl
```

Expected installation and precompilation time is 5 minutes.

## Demo

For a demo, run the file `demo.jl` as:

```console
$ julia demo.jl
```

Assuming the command line is in the directory of the repository. Alternatively, run the script using the relative path to the repository as:

```console
$ julia path/to/repo/demo.jl
```

Or the following code in Julia:

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
- For the independent evolution of social norms and strategies run: `independent.jl`.
- For the evolutionary dynamics between $DISC_{q,M}$ and $probC$ run: `probabilistic.jl`.

All data used to produce figures of the paper is available in the `results/` folder.

The Mathematica notebook to plot the panels is available in the `reproduce_figures/` folder.
The panels are in the `reproduce_figures/figures/` folder.

---
Copyright (c) 2023 Sebastian Michel-Mata
