# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:hydrogen
#     text_representation:
#       extension: .jl
#       format_name: hydrogen
#       format_version: '1.3'
#       jupytext_version: 1.10.3
#   kernelspec:
#     display_name: Julia 1.9.0
#     language: julia
#     name: julia-1.9
# ---

# %%
using Random
using BenchmarkTools
versioninfo()

# %%
rng = @btime Random.Xoshiro(4649373)
@btime rand($rng, 2^20);

# %%
rng = @btime Random.MersenneTwister(4649373)
@btime rand($rng, 2^20);

# %%
rng = @btime Random.Xoshiro(4649373)
@btime rand($rng, UInt64, 2^20);

# %%
rng = @btime Random.MersenneTwister(4649373)
@btime rand($rng, UInt64, 2^20);

# %%
