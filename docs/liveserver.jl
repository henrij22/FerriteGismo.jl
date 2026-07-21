#!/usr/bin/env julia

# Note, Revise.jl and LiveServer.jl do not come with this environment.
# It is advised to install them in the global environment.

using Revise
# Root of the repository
const repo_root = dirname(@__DIR__)

# Make sure docs environment is active and instantiated
import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

# Communicate with docs/make.jl that we are running in live mode
push!(ARGS, "liveserver")

# Run LiveServer.servedocs(...)
import LiveServer
LiveServer.servedocs(;
    host = "0.0.0.0",
    # Documentation root where make.jl and src/ are located
    foldername = joinpath(repo_root, "docs"),
    # Extra source folder to watch for changes so docstrings can be Revise'd
    include_dirs = [
        joinpath(repo_root, "src"),
    ]
)
