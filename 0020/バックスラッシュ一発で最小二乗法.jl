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

# %%
using Plots

f(x) = 1 - 2x
x = range(0, 2; length=21)
y = @. f(x) + randn()
X = x .^ (0:1)'

β̂ = X \ y # \beta TAB \hat TAB
g(x) = β̂[1] + β̂[2]*x

xs = range(extrema(x)...; length=300)
scatter(x, y; label="data", color=1)
plot!(xs, f; label="true line", color=1, ls=:dash)
plot!(xs, g; label="regression line", color=2, lw=2)

# %%
x = range(0, 2; length=21)

# %%
collect(x)

# %%
X = x .^ (0:1)'

# %%
x .^ (0:3)'

# %%
X \ y

# %%
using LinearAlgebra
qr(X, ColumnNorm()) \ y

# %%
(X'X)\X'y

# %% tags=[]
@which X \ y

# %% [markdown]
# https://github.com/JuliaLang/julia/blob/master/stdlib/LinearAlgebra/src/generic.jl#L1128

# %%
pinv(X)y

# %%
pinv(X)

# %%
U, S, Vt = svd(X)
Vt' * Diagonal(inv.(S)) * U'

# %%
@which pinv(X)

# %% [markdown]
# https://github.com/JuliaLang/julia/blob/master/stdlib/LinearAlgebra/src/dense.jl#L1430

# %% [markdown]
# https://discourse.julialang.org/t/differences-in-a-b-for-sparse-and-nonsparse-rank-deficient-a/66917

# %%
using LinearAlgebra, SparseArrays

Is = [9, 8, 9, 12, 2, 5, 6, 11, 8, 12]
Js = [1, 3, 4, 4, 5, 5, 7, 7, 10, 10]
Vs = [0.046668312772398135, 0.0732172552606527, 0.07228019735547542, 0.27500506619596665, 0.8821260519923395, 0.188438924730004, 0.22883463110234215, 0.17705262933291666, 0.297278962545636, 0.07770904672896628]
A = sparse(Is, Js, Vs)
b = [0.0
 0.26743848063303416
 0.0
 0.0
 0.057129952809003536
 0.20824640373847988
 0.0
 0.1580262484804872
 0.09667497046894227
 0.0
 0.16112322314768995
 0.27723584012581337]

@show qr(A) \ b ≈ A \ b # True
@show A \ b ≈ pinv(Matrix(A))*b # False
@show qr(Matrix(A), Val(true)) \ b ≈ pinv(Matrix(A))*b # True;

# %%
typeof(A), typeof(b)

# %%
@show A \ b;

# %%
@show Matrix(A) \ b;

# %%
@show pinv(Matrix(A)) * b;

# %%
