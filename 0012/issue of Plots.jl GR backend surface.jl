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

# %% [markdown]
# https://discourse.julialang.org/t/wrong-image-when-trying-to-create-a-simple-3d-surface/64767

# %%
using Plots

x = collect(0:0.1:1);
y = collect(0:0.1:1);
z = collect(0:0.1:1);
function make_plane(x,y,z)
    a = []
    b = []
    c = []
    for i in x, j in y, k in z
        if (i + j + k) == 1
            push!(a,i)
            push!(b,j)
            push!(c,k)
        end
    end
    return a,b,c
end

plane = make_plane(x,y,z)

plot(plane, camera = (60,60)) #I know this is wrong, but adding this 2nd for reference!-

# %%
using Plots
gr()

function make_plane(x, y, z)
    a = eltype(x)[]
    b = eltype(y)[]
    c = eltype(z)[]
    for i in x, j in y, k in z
        if i + j + k ≈ 1
            push!(a, i)
            push!(b, j)
            push!(c, k)
        end
    end
    return a,b,c
end

x = y = z = 0:0.1:1
surface(make_plane(x, y, z), camera = (60, 60), c = :blues)

# %%
gr()
f(x, y) = x ≥ 0 && y ≥ 0 && x + y ≤ 1 ? 1 - x - y : NaN
x = y = 0:0.1:1
surface(x, y, f; camera = (60, 60), c = :blues)

# %%
pyplot(fmt = :svg)
x = y = z = 0:0.1:1
surface(make_plane(x, y, z), camera = (70, 20), c = :blues)

# %%
