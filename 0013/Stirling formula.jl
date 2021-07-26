# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:hydrogen
#     text_representation:
#       extension: .jl
#       format_name: hydrogen
#       format_version: '1.3'
#       jupytext_version: 1.11.2
#   kernelspec:
#     display_name: Julia 1.8.0-DEV
#     language: julia
#     name: julia-1.8
# ---

# %%
using SpecialFunctions, Plots

logfact(n) = loggamma(n + 1)
stirling(n, ::Val{1}) = n*log(n)
stirling(n, ::Val{2}) = n*log(n) - n
stirling(n, ::Val{3}) = n*log(n) - n + (1/2)*log(n)
stirling(n, ::Val{4}) = n*log(n) - n + (1/2)*log(n) + log(√(2π))

PP = []
for (m, legend) in ((0, :bottomright), (1, :topleft), (2, :topleft), (23, :topleft))
    xmax = 10.0^m
    P = plot(; legend, xlim=(0, 1.05xmax), tickfontsize=6)
    plot!(logfact, 0.01, 1.05xmax; label="log n!")
    for k in 1:4
        plot!(x -> stirling(x, Val(k)), 0.01, 1.05xmax; label="stirling $k", ls=:auto)
    end
    push!(PP, P)
end
plot(PP...; size=(800, 600))

# %% [markdown]
# ```julia
# logfact(n) = loggamma(n + 1)
# stirling(n, ::Val{1}) = n*log(n)
# stirling(n, ::Val{2}) = n*log(n) - n
# stirling(n, ::Val{3}) = n*log(n) - n + (1/2)*log(n)
# stirling(n, ::Val{4}) = n*log(n) - n + (1/2)*log(n) + log(√(2π))
# ```

# %%
