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
#     display_name: Julia 1.10.0
#     language: julia
#     name: julia-1.10
# ---

# %%
@time sum(randn(10^8)); Base.GC.gc()
@time sum(randn(10^8)); Base.GC.gc()
@time sum(randn(10^8))

# %%
