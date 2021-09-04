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
#     display_name: Julia 1.6.2
#     language: julia
#     name: julia-1.6
# ---

# %% [markdown]
# https://cometscome.github.io/JuliaFromFortran/build/chapter2/01/

# %% [markdown]
# https://docs.julialang.org/en/v1/manual/variables-and-scoping/#Loops-and-Comprehensions

# %%
function f()
    i = -99
    for i in 1:3 end
    i
end

@show f();

# %%
function g()
    local i
    for outer i in 1:3 end
    i
end

@show g();

# %%
function mandel()
    nx = 61; ny = 31; maxiter = 90
    x0 = -2.0; x1 = 2.0; y0 = -2.0; y1 = 2.0
    mandelbrot = zeros(Int64,nx, ny)
    mandelbrot .= maxiter
    for iy=1:ny
        y = y0 + (y1 - y0) * (iy - 1) / real(ny - 1) 
        for ix=1:nx
            x = x0 + (x1 - x0) * (ix - 1) / real(nx - 1)
            c = x + im*y
            z = 0im
            for iter=0:maxiter
                z = z * z + c
                if abs(z) > 2
                    mandelbrot[ix,iy] = iter
                    break
                end
            end
            
        end
    end
    for iy=1:ny
        for ix=1:nx
            if mandelbrot[ix,iy] == maxiter
                print("*")
            else
                print((mandelbrot[ix,iy]+9) ÷ 10)
            end
        end
        println("\t")
    end
end
mandel()

# %%
function mandel_outer()
    nx = 61; ny = 31; maxiter = 90
    x0 = -2.0; x1 = 2.0; y0 = -2.0; y1 = 2.0
    mandelbrot = zeros(Int64,nx, ny)
    for iy=1:ny
        y = y0 + (y1 - y0) * (iy - 1) / real(ny - 1) 
        for ix=1:nx
            x = x0 + (x1 - x0) * (ix - 1) / real(nx - 1)
            c = x + im*y
            z = 0im
            local iter; for outer iter=0:maxiter
                z = z * z + c
                abs(z) > 2 && break
            end
            mandelbrot[ix,iy] = iter
        end
    end
    for iy=1:ny
        for ix=1:nx
            if mandelbrot[ix,iy] == maxiter
                print("*")
            else
                print((mandelbrot[ix,iy]+9) ÷ 10)
            end
        end
        println("\t")
    end
end
mandel_outer()

# %%
