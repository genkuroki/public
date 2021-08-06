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

function mandelbrot(c; maxiters=2^10, threshold_abs2=2^10)
    z = zero(c)
    for i in 1:maxiters
        z = z * z + c
        abs2(z) ≥ threshold_abs2 && return i
    end
    maxiters + 1
end

n = 2^8
x = range(-0.714689, -0.714679; length=n)
y = range( 0.299872,  0.299882; length=n)
c = complex.(x', y)

@show typeof(c)
@time m = mandelbrot.(c)
@time m = mandelbrot.(c)
@show typeof(m)
plotmandelbrot(m)

# %%
x32 = range(-0.714689f0, -0.714679f0; length=n)
y32 = range( 0.299872f0,  0.299882f0; length=n)
c32 = complex.(x32', y32)

@show typeof(c32)
@time m32 = mandelbrot.(c32)
@time m32 = mandelbrot.(c32)
@show typeof(m32)
plotmandelbrot(m32)

# %%
using CUDA
c_cuda = cu(c32)

@show typeof(c_cuda)
@time m_cuda = collect(mandelbrot.(c_cuda))
@time m_cuda = collect(mandelbrot.(c_cuda))
@show typeof(mandelbrot.(c_cuda))
@show typeof(collect(mandelbrot.(c_cuda)))
plotmandelbrot(m_cuda)

# %%
function mandelbrot_threads(c)
    m = similar(c, typeof(mandelbrot(c[end])))
    Threads.@threads for i in keys(c)
        m[i] = mandelbrot(c[i])
    end
    m
end

@show Threads.nthreads()
@time m_th = mandelbrot_threads(c)
@time m_th = mandelbrot_threads(c)
plotmandelbrot(m_th)

# %%
@benchmark mandelbrot.($c) # CPU Float64

# %%
@benchmark mandelbrot.($c32) # CPU Float32

# %%
@benchmark collect(mandelbrot.($c_cuda)) # GPU Float32 (GPU → GPU)

# %%
@benchmark collect(mandelbrot.(cu($c))) # GPU Float32 (CPU → GPU → GPU → CPU)

# %%
@show Threads.nthreads()
@benchmark mandelbrot_threads($c) # CPU Float64 multi-threaded version with nthreads = 12

# %%
