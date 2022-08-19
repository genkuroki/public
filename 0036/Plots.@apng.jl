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
#     display_name: Julia 1.9.0-DEV
#     language: julia
#     name: julia-1.9
# ---

# %% [markdown]
# https://github.com/JuliaPlots/Plots.jl/pull/4295

# %%
using Plots
using Plots: apng, @apng
using Random

# %%
@apng for a in 1:100
    plot(x -> sin(a*x); label="")
end

# %%
anim = @animate for i in 1:50
    Random.seed!(123)
    scatter(cumsum(randn(i)), ms=i, lab="", alpha = 1 - i/50,
        xlim=(0,50), ylim=(-5, 7))
end

apng(anim)

# %%
@apng for i in 1:50
    Random.seed!(123)
    scatter(cumsum(randn(i)), ms=i, lab="", alpha = 1 - i/50,
        xlim=(0,50), ylim=(-5, 7))
end

# %%
