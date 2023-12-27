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
#     display_name: Julia 1.10.0
#     language: julia
#     name: julia-1.10
# ---

# %%
function f(x_prev, v, dt)
    x_prev + v*dt
end

@code_native debuginfo=:none dump_module=false f(1.0, 2.0, 0.1)

# %%
@fastmath function f_fastmath(x_prev, v, dt)
    x_prev + v*dt
end

@code_native debuginfo=:none dump_module=false f_fastmath(1.0, 2.0, 0.1)

# %%
"""
Solve x(0) = 0, x'(t) = x(t) + 1, 0 ≤ t ≤ 1.

Exact solution: x(t) = exp(t) - 1, x(1) = exp(1) - 1 = 1.718281828⋯.
"""
function solveode(f, n=10^6, x0=0.0, dt=1/n)
    x = x0
    for _ in 1:n
        x = f(x, x+1, dt)
    end
    x
end

@show solveode(f)
@show solveode(f_fastmath)
@show exp(1) - 1;

# %%
using BenchmarkTools

@btime solveode(f)
@btime solveode(f_fastmath);

# %%
@code_native debuginfo=:none dump_module=false solveode(f)

# %%
@code_native debuginfo=:none dump_module=false solveode(f_fastmath)

# %%
