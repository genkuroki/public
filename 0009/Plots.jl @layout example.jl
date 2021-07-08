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
using Plots
default(legend = false, tickfontsize=6)

A = plot(cumsum(randn(100)))
B = heatmap(randn(30, 10))
C = heatmap(randn(50, 50))
D = plot([cumsum(randn(100)) cumsum(randn(100))])
E = heatmap(randn(10, 30))
F = plot(cumsum(randn(50)))

layout = @layout [
    [a; b] d{0.4w, 0.4h}
    c      [e f]
]
plot(A, B, D, C, E, F; layout)

# %%
