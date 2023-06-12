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
#     display_name: Julia 1.9.1
#     language: julia
#     name: julia-1.9
# ---

# %% [markdown]
# https://twitter.com/yosaytouchgot/status/1668240209354592257

# %%
m = Threads.nthreads()

# %%
using Plots
default(fmt=:png)
using SharedArrays
using Distributed

using SymPy
@show x = symbols("x")
@show F = sin(x)
@show G = F.integrate(x)
@show Gexpr = Meta.parse(string(G))

addprocs(m)
@everywhere begin
    @eval g(x) = $Gexpr
end

function makedata(m = m)
    data = SharedArray{Float64}(100, m)
    t = range(0, 2pi, 101)[1:end-1]
    @sync @distributed for i in 1:m
        data[:,i] .= g.(i .* t)
    end
    data
end

data = @time makedata(m)
data = @time makedata(m)
data = @time makedata(m)
rmprocs(workers())

plot()
for i in 1:4 plot!(data[:,i]) end
plot!()

# %%
using Plots
default(fmt=:png)

using SymPy
@show x = symbols("x")
@show F = sin(x)
@show G = F.integrate(x)
@show Gexpr = Meta.parse(string(G))
@eval g(x) = $Gexpr

function makedata(m = m)
    data = Matrix{Float64}(undef, 100, m)
    t = range(0, 2pi, 101)[1:end-1]
    for i in 1:m
        data[:,i] .= g.(i .* t)
    end
    data
end

data = @time makedata(m)
data = @time makedata(m)
data = @time makedata(m)

plot()
for i in 1:4 plot!(data[:,i]) end
plot!()

# %%
using Plots
default(fmt=:png)

using SymPy
@show x = symbols("x")
@show F = sin(x)
@show G = F.integrate(x)
@show Gexpr = Meta.parse(string(G))
@eval g(x) = $Gexpr

function makedata(m = m)
    data = Matrix{Float64}(undef, 100, m)
    t = range(0, 2π, 101)[1:end-1]
    Threads.@threads for i in 1:m
        data[:,i] .= g.(i .* t)
    end
    data
end

data = @time makedata(m)
data = @time makedata(m)
data = @time makedata(m)
plot()
for i in 1:4 plot!(data[:,i]) end
plot!()

# %%
using Plots
default(fmt=:png)

using SymPy
@show x = symbols("x")
@show F = sin(x)
@show G = F.integrate(x)
@show Gexpr = Meta.parse(string(G))
@eval g(x) = $Gexpr

function makedata(m = m)
    t = range(0, 2π, 101)[1:end-1]
    g.((1:m)'.* t)
end

data = @time makedata(m)
data = @time makedata(m)
data = @time makedata(m)
plot()
for i in 1:4 plot!(data[:,i]) end
plot!()

# %%
