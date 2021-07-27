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
# https://discourse.julialang.org/t/fast-performance-of-array-comprehension-without-allocations/65352
#
# https://github.com/JuliaSIMD/LoopVectorization.jl/issues/307

# %%
println("Julia v$VERSION")
using Pkg
Pkg.status("SLEEFPirates")
Pkg.status("LoopVectorization")

# %%
using BenchmarkTools
using LoopVectorization

function assign!(u, x, y, z)
    @turbo for k in 1:length(z)
        for j in 1:length(y)
            for i in 1:length(x)
                u[i, j, k] = sin(x[i]) + sin(y[j]) + sin(z[k])
            end 
        end 
    end 
end

function assign_broadcast!(uv, X, Y, Z, xx, yy, zz)
    @. X = sin(xx)
    @. Y = sin(yy)
    @. Z = sin(zz)
    @. uv = X + Y + Z
end

itot = 384 
dx = 1. / itot
x = dx*collect(0:itot-1); y = dx*collect(0:itot-1); z = dx*collect(0:itot-1)
u = zeros(itot+8, itot+8, itot+8)

uv = @view u[5:5+itot-1, 5:5+itot-1, 5:5+itot-1]
xx, yy, zz = reshape(x, (:, 1, 1)), reshape(y, (1, :, 1)), reshape(z, (1, 1, :))
X, Y, Z = similar(xx), similar(yy), similar(zz)

assign!(uv, x, y, z)
a = deepcopy(uv)
uv[:, :, :] = [ sin(x) + sin(y) + sin(z) for x=x, y=y, z=z ]
b = deepcopy(uv)
@. uv = sin(xx) + sin(yy) + sin(zz)
c = deepcopy(uv)
assign_broadcast!(uv, X, Y, Z, xx, yy, zz)
d = deepcopy(uv)
@show a ≈ b ≈ c ≈ d

print("LoopVectorization.@turbo:")
@btime assign!($uv, $x, $y, $z)
print("comprehension:           ")
@btime $uv[:, :, :] = [ sin(x) + sin(y) + sin(z) for x=$x, y=$y, z=$z ]
print("simple broadcast:        ")
@btime @. $uv = sin($xx) + sin($yy) + sin($zz)
print("assign_broadcast!:       ")
@btime assign_broadcast!($uv, $X, $Y, $Z, $xx, $yy, $zz);

# %%
