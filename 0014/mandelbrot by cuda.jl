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
# https://nbviewer.jupyter.org/gist/antimon2/dd23b1e62d5fcd08d91fd420d7ba996d/Mandelbrot.mt.cu.jl.ipynb

# %%
using BenchmarkTools
using Plots
gr(fmt=:png)

plotmandelbrot(m) = heatmap(m; c=reverse(cgrad(:jet1)), 
    size=(300, 300), colorbar=false, ticks=false, frame=false)

function mandelbrot(c; maxiters=2^10, infabs2=2^10)
    z = zero(c)
    for i in 1:maxiters
        z = z * z + c
        abs2(z) ≥ infabs2 && return i
    end
    maxiters + 1
end

n = 2^8
x = range(-0.714689, -0.714679; length=n)
y = range( 0.299872,  0.299882; length=n)
z = complex.(x', y)

@time m = mandelbrot.(z)
@time m = mandelbrot.(z)
plotmandelbrot(m)

# %%
x32 = range(-0.714689f0, -0.714679f0; length=n)
y32 = range( 0.299872f0,  0.299882f0; length=n)
z32 = complex.(x32', y32)

@time m32 = mandelbrot.(z32)
@time m32 = mandelbrot.(z32)
plotmandelbrot(m32)

# %%
using CUDA
z_cuda = cu(z32)

@time m_cuda = collect(mandelbrot.(z_cuda))
@time m_cuda = collect(mandelbrot.(z_cuda))
plotmandelbrot(m_cuda)

# %%
@benchmark mandelbrot.($z)

# %%
@benchmark mandelbrot.($z32)

# %%
@benchmark collect(mandelbrot.($z_cuda))

# %%
@benchmark collect(mandelbrot.(cu($z)))

# %%
