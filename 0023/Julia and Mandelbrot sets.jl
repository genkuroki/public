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
#     display_name: Julia 1.6.3
#     language: julia
#     name: julia-1.6
# ---

# %%
using BenchmarkTools
using Plots
gr(fmt=:png)
using CUDA

# %%
plotjulia(j) = heatmap(j; c=:gist_earth,
    size=(540, 540), colorbar=false, ticks=false, frame=false)

plotmandelbrot(m) = heatmap(m; c=reverse(cgrad(:jet1)), 
    size=(540, 540), colorbar=false, ticks=false, frame=false)

function julia(z, c, maxiters=2^10, threshold_abs2=Inf)
    for i in Base.OneTo(maxiters)
        z = z * z + c
        abs2(z) ≥ threshold_abs2 && return i
    end
    maxiters + oneunit(maxiters)
end

mandelbrot(c, maxiters=2^10, threshold_abs2=Inf) =
    julia(zero(c), c, maxiters, threshold_abs2)

# %%
n = 2^9
x = range(-1.5, 1.5; length=n)
y = range(-1.5, 1.5; length=n)
z = complex.(x', y)
c = -0.380 + 0.605im

@show typeof(z)
@time j = julia.(z, c, 2^8)
@time j = julia.(z, c, 2^8)
@show typeof(j)
plotjulia(j)

# %%
n = 2^9
x = range(-1.5, 1.5; length=n)
y = range(-1.5, 1.5; length=n)
z = complex.(x', y)
c = -0.380 + 0.605im

c32 = ComplexF32(c)
z_cuda32 = cu(z)

@show typeof(z_cuda32)
@time j_cuda32 = julia.(z_cuda32, c32, Int32(2^8), Inf32)
@time j_cuda32 = julia.(z_cuda32, c32, Int32(2^8), Inf32)
@show typeof(j_cuda32)
@time cj_cuda32 = collect(j_cuda32)
plotjulia(cj_cuda32)

# %%
n = 2^9
x = range(-1.5, 1.5; length=n)
y = range(-1.5, 1.5; length=n)
z = complex.(x', y)
c = -0.380 + 0.605im

z_cuda = CuMatrix{ComplexF64}(z)

@show typeof(z_cuda)
@time j_cuda = julia.(z_cuda, c, 2^8)
@time j_cuda = julia.(z_cuda, c, 2^8)
@show typeof(j_cuda)
@time cj_cuda = collect(j_cuda)
plotjulia(cj_cuda)

# %%
n = 2^9
x = range(-1.6, 0.5; length=n)
y = range(-1.1, 1.1; length=n)
c = complex.(x', y)

@show typeof(c)
@time m = mandelbrot.(c, 2^6)
@time m = mandelbrot.(c, 2^6)
@show typeof(m)
plotmandelbrot(m)

# %%
n = 2^9
x = range(-1.6, 0.5; length=n)
y = range(-1.1, 1.1; length=n)
c = complex.(x', y)
c_cuda = CuMatrix{ComplexF64}(c)

@show typeof(c_cuda)
@time m_cuda = mandelbrot.(c_cuda, 2^6)
@time m_cuda = mandelbrot.(c_cuda, 2^6)
@show typeof(m_cuda)
@time cm_cuda = collect(m_cuda)
@show typeof(cm_cuda)
plotmandelbrot(cm_cuda)

# %%
n = 2^9
x = range(-0.714689, -0.714679; length=n)
y = range( 0.299872,  0.299882; length=n)
c = complex.(x', y)

@show typeof(c)
@time m = mandelbrot.(c, 2^12)
@time m = mandelbrot.(c, 2^12)
@show typeof(m)
plotmandelbrot(m)

# %%
n = 2^9
x = range(-0.714689, -0.714679; length=n)
y = range( 0.299872,  0.299882; length=n)
c = complex.(x', y)
c_cuda = CuMatrix{ComplexF64}(c)

@show typeof(c_cuda)
@time m_cuda = mandelbrot.(c_cuda, 2^12)
@time m_cuda = mandelbrot.(c_cuda, 2^12)
@show typeof(m_cuda)
@time cm_cuda = collect(m_cuda)
@show typeof(cm_cuda)
plotmandelbrot(cm_cuda)

# %%
