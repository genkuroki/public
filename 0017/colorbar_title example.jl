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
#     display_name: Julia 1.6.2
#     language: julia
#     name: julia-1.6
# ---

# %%
using Plots
p = heatmap(rand(100, 100);
    colorbar_title = "colorbar title",
    thickness_scaling = 1.6,
)

# %%
using Plots
p = heatmap(rand(100, 100);
    colorbar_title = " \ncolorbar title",
    thickness_scaling = 1.6,
    right_margin = 2Plots.mm,
)

# %%
