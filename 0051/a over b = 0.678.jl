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
#     display_name: Julia 1.11.1
#     language: julia
#     name: julia-1.11
# ---

# %%
ENV["LINES"] = 1000
using Plots

# %%
[(a, b, round(a/b; digits=3)) for a in 1:80 for b in a:80 if 0.67 < a/b < 0.68]

# %%
[(a, b, round(a/b; digits=3)) for a in 1:80 for b in a:80 if round(a/b; digits=3) == 0.678]

# %%
x = 1:80
y = 1:80
z = log.(y ./ x')
heatmap(x, y, z)
plot!(x, 0.678x; label="")
scatter!([(59, 40)]; label="")

# %%
z .|> exp

# %%
