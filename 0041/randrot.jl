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
#     display_name: Julia 1.9.0-beta3
#     language: julia
#     name: julia-1.9
# ---

# %%
using StaticArrays
using Plots

function randrot(σ = 0.01)
    x, y, z = σ*randn(), σ*randn(), σ*randn()
    A = SMatrix{3, 3}(0, -z, y, z, 0, -x, -y, x, 0)
    exp(A)
end

function anim_randrot(L = 200, σ = 0.2)
    v = SVector(1.0, 0.0, 0.0)
    anim = @animate for _ in 1:L
        plot3d([0, v[1]], [0, v[2]], [0, v[3]]; label="", lw=3)
        plot!(xlim=(-1.1, 1.1), ylim=(-1.1, 1.1), zlim=(-1.1, 1.1))
        plot!(size=(500, 450))
        v = randrot(σ)*v
    end
    gif(anim, "randrot.gif")
end

# %%
U = randrot()

# %%
U'U

# %%
V = randrot(1)

# %%
V'V

# %%
anim_randrot()

# %%
