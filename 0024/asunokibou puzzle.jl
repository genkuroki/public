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
#     display_name: Julia 1.6.3
#     language: julia
#     name: julia-1.6
# ---

# %% [markdown]
# * https://twitter.com/asunokibou/status/1457694534016507922
# * https://twitter.com/mathsorcerer/status/1457937816235634691
# * https://twitter.com/sophia8wendel/status/1457980892966576128

# %%
ENV["LINES"] = 100
using Combinatorics
sol = [p for p in permutations(1:9, 7) if sum(i//p[i] for i in 1:7) == 7]

# %%
