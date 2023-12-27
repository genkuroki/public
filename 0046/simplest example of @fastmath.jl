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
function repeatf(f, n=10^6, x0=0.0, v=1.0, dt=1/n)
    x = x0
    for _ in 1:n
        x = f(x, v, dt)
    end
    x
end

@show repeatf(f)
@show repeatf(f_fastmath);

# %%
using BenchmarkTools

@btime repeatf(f)
@btime repeatf(f_fastmath);

# %%
