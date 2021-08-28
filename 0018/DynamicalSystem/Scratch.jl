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

# %%
using Base.Threads
using ProgressMeter
using Plots
nthreads()
xoxp1(x) = x/(x+1) # {x > 0} -> {0 < y < 1}

# %%
using DynamicalSystems
using StaticArrays

# %%
function f(x, p)
   r, k = p
   a, mu, d = 5.0, 0.5, 0.2
   x1new = x[1]*exp(r*(1-x[1]/k)-x[2]/(a+x[1]^2))
   x2new = x[2]*exp(mu*x[1]/(a+x[1]^2)-d)
   x1new, x2new
end

f(x, p, t) = SVector(f(x, p))

# %%
function traj2d!(tr1, tr2, f, p, x0, niters, nwarmup)
    x = x0
    for _ in 1:nwarmup
        x = f(x, p)
    end
    @inbounds tr1[1], tr2[1] = x
    for i in 2:niters+1
        x = f(x, p)
        @inbounds tr1[i], tr2[i] = x
    end
end

# %%
function traj2d_ds!(tr1, tr2, dds, niters, nwarmup)
    integ = integrator(dds, dds.u0)
    nwarmup != 0 && step!(integ, nwarmup)
    @inbounds tr1[1], tr2[1] = DynamicalSystems.DynamicalSystemsBase.obtain_access(integ.u, nothing)
    for i in 2:niters+1
        step!(integ)
        @inbounds tr1[i], tr2[i] = DynamicalSystems.DynamicalSystemsBase.obtain_access(integ.u, nothing)
    end
end

# %%
niters = 10
nwarmup = 10
tr1 = zeros(Int, niters+1)
tr2 = similar(tr1)
traj2d!(tr1, tr2, ((a, b), p) -> (b, a+b), nothing, (0, 1), niters, nwarmup)
@show tr1 tr2;

# %%
function findperiod(x; tol=1e-3)
    n = length(x)
    @inbounds for k in 2:(n ÷ 2 + 1)
        if abs(x[k] - x[1]) ≤ tol
            all(j -> abs(x[j] - x[j-k+1]) ≤ tol, k:n) && return k - 1
        end
    end 
    return n
end

# %%
x = [1, 2, 3, 4, 1, 2, 3]
findperiod(x)

# %%
x = [1, 2, 3, 4, 1, 2, 3, 4]
findperiod(x)

# %%
x = [1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3]
findperiod(x)

# %%
x = @. sinpi(0:0.1:10) + 1e-4randn()
findperiod(x)

# %%
niters = 2000
nwarmup = 50000
tr1_ds = zeros(niters+1)
tr2_ds = similar(tr1)
x0 = (0.4, 0.5)
p = (2.5, 4.0)
dds = DiscreteDynamicalSystem(f, SVector(x0), collect(p))
traj2d_ds!(tr1_ds, tr2_ds, dds, niters, nwarmup)
@show tr1_ds[1:10];

# %%
niters = 2000
nwarmup = 50000
tr1 = zeros(niters+1)
tr2 = similar(tr1)
x0 = (0.4, 0.5)
p = (2.5, 4.0)
traj2d!(tr1, tr2, f, p, x0, niters, nwarmup)
@show tr1[1:10];

# %%
tr1_ds == tr1

# %%
@btime traj2d!($tr1, $tr2, f, p, x0, niters, nwarmup)
@btime traj2d_ds!($tr1_ds, $tr2_ds, dds, niters, nwarmup);

# %%
function calcperiod2d!(tr1, tr2, f, p, x0, niters, nwarmup; tol=1e-3)
    traj2d!(tr1, tr2, f, p, x0, niters, nwarmup)
    findperiod(tr1; tol)
end

# %%
using BenchmarkTools
@btime calcperiod2d!($(zeros(2001)), $(zeros(2001)), f, (2.6, 4.0), (0.4, 0.5), 2000, 50000; tol=1e-3)

# %%
function calcperiod2d(f = f, niters = 2000, nwarmup = 50000, x0 = (0.4, 0.5),
        r = (1:0.05:5)[1:end-1], k = (2:0.05:5)[1:end-1], tol=1e-3
    )
    tr1 = Vector{Float64}(undef, niters+1)
    tr2 = similar(tr1)
    period = Matrix{Int}(undef, length(r), length(k))
    prog = Progress(length(r)*length(k))
    for j in eachindex(k)
        for i in eachindex(r)
            @inbounds period[i, j] = calcperiod2d!(tr1, tr2, f, (r[i], k[j]), x0, niters, nwarmup; tol)
            next!(prog)
        end
    end
    period, r, k
end

period, r, k = @time calcperiod2d()
heatmap(r, k, xoxp1.(period'); xlabel="r", ylabel="k", title="period/(period + 1)")

# %%
function calcperiod2d_threads(f = f, niters = 2000, nwarmup = 50000, x0 = (0.4, 0.5),
        r = (1:0.05:5)[1:end-1], k = (2:0.05:5)[1:end-1], tol=1e-3
    )
    tr1 = [Vector{Float64}(undef, niters+1) for _ in 1:nthreads()]
    tr2 = [Vector{Float64}(undef, niters+1) for _ in 1:nthreads()]
    period = Matrix{Int}(undef, length(r), length(k))
    #prog = Progress(length(r)*length(k))
    @threads for j in eachindex(k)
        tid = threadid()
        for i in eachindex(r)
            @inbounds period[i, j] = calcperiod2d!(tr1[tid], tr2[tid], f, (r[i], k[j]), x0, niters, nwarmup; tol)
            #next!(prog)
        end
    end
    period, r, k
end

# %%
period, r, k_th, r, k = @time calcperiod2d_threads(); flush(stdout)
period, r, k_th, r, k = @time calcperiod2d_threads(); flush(stdout)
period, r, k_th, r, k = @time calcperiod2d_threads(); flush(stdout)
heatmap(r, k, xoxp1.(period_th'); xlabel="r", ylabel="k", title="period/(period + 1)")

# %%
period, r, k_th = @time calcperiod2d_threads(); flush(stdout)
period, r, k_th = @time calcperiod2d_threads(); flush(stdout)
period, r, k_th = @time calcperiod2d_threads(); flush(stdout)
heatmap(r, k, xoxp1.(period_th'); xlabel="r", ylabel="k", title="period/(period + 1)")

# %%
period == period_th

# %%
period'

# %%
@btime calcperiod2d_threads(f; niters = 2000, nwarmup = 50000, x0 = (0.4, 0.5),
    r = (1:0.05:5)[1:end-1], k = (2:0.05:5)[1:end-1], tol=1e-3)

# %%
@code_warntype calcperiod2d_threads(f; niters = 2000, nwarmup = 50000, x0 = (0.4, 0.5),
    r = (1:0.05:5)[1:end-1], k = (2:0.05:5)[1:end-1], tol=1e-3)

# %%
