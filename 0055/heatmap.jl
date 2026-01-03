# -*- coding: utf-8 -*-
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
#     display_name: Julia
#     language: julia
#     name: julia
# ---

# %% [markdown]
# https://x.com/dannchu/status/2007050063747666414
#
# https://colab.research.google.com/github/genkuroki/public/blob/main/0055/heatmap.ipynb

# %%
using Plots

f(x, y) = 1 ≤ abs(abs(x) - 2) + abs(abs(y) - 2) ≤ 3
xs = -6:0.01:6
ys = -6:0.01:6
heatmap(xs, ys, f; aspectratio=:equal, c=:binary, colorbar=false)
plot!(xlim=extrema(xs), ylim=extrema(ys))
plot!(xguide="x", yguide="y")
title!("Region: 1 ≤ |│x│-2| + |│y│-2| ≤ 3")
plot!(size=(400, 400))

# %%
