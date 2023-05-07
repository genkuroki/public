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
#     display_name: Julia 1.8.5
#     language: julia
#     name: julia-1.8
# ---

# %%
versioninfo()

# %%
using Plots
default(fmt = :png)

# %%
Base.@kwdef struct Param{I, F}
    xmin::F = -1.75
    xmax::F = 0.75
    width::I = 4096
    ymin::F = -1.25
    ymax::F = 1.25
    height::I = 4096
    max_iter::I = 500
end

function makeaxes(param::Param)
    (; xmin, xmax, width, ymin, ymax, height) = param
    x = range(xmin, xmax, width + 1)[begin:end-1]
    y = range(ymin, ymax, height + 1)[begin:end-1]
    x, y
end

function makeinput(param::Param)
    x, y = makeaxes(param)
    complex.(x, y')
end

function makeoutput(param::Param)
    (; width, height) = param
    zeros(Int32, width, height)
end

function mandelbrot_kernel(param::Param, c)::Int32
    (; max_iter) = param
    z = c
    for i in 0:max_iter-1
        z = z * z + c
        # zが閾値を超えたら終了します
        abs(z) > 2 && return i 
    end
    max_iter
end

function compute_mandelbrot!(param::Param, input, output)
    output .= mandelbrot_kernel.((param,), input)
end

@time param = Param()
@time z = makeinput(param)
@time image = makeoutput(param)
@time compute_mandelbrot!(param, z, image)
@show typeof(z) typeof(image);

# %%
@time compute_mandelbrot!(param, z, image);

# %%
@time compute_mandelbrot!(param, z, image);

# %%
x, y = makeaxes(param)
m = image'
heatmap(x, y, m; size=(800, 720), clim=(0, 30), xguide="x", yguide="y")

# %%
using CUDA

# %%
@time param = Param()
@time z_cuda64 = CuArray{ComplexF64, 2}(makeinput(param))
@time image_cuda = cu(makeoutput(param))
@time compute_mandelbrot!(param, z_cuda64, image_cuda)
@show typeof(z_cuda64) typeof(image_cuda);

# %%
@time compute_mandelbrot!(param, z_cuda64, image_cuda);

# %%
@time compute_mandelbrot!(param, z_cuda64, image_cuda);

# %%
x, y = makeaxes(param)
m = Matrix(image_cuda)'
heatmap(x, y, m; size=(800, 720), clim=(0, 30), xguide="x", yguide="y")

# %%
@time z_cuda32 = cu(makeinput(param))
@time compute_mandelbrot!(param, z_cuda32, image_cuda)
@show typeof(z_cuda32) typeof(image_cuda);

# %%
@time compute_mandelbrot!(param, z_cuda32, image_cuda);

# %%
@time compute_mandelbrot!(param, z_cuda32, image_cuda);

# %%
