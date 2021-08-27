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

# %%
function f(x, p)
   r, k = p
   a, mu, d = 5.0, 0.5, 0.2
   x1new = x[1]*exp(r*(1-x[1]/k)-x[2]/(a+x[1]^2))
   x2new = x[2]*exp(mu*x[1]/(a+x[1]^2)-d)
   x1new, x2new
end

# %%
function traj2d!(tr1, tr2, f, p, x0, lentr, nskips)
    x = x0
    for _ in 1:nskips
        x = f(x, p)
    end
    @inbounds tr1[1], tr2[1] = x
    for i in 1:lentr
        @inbounds tr1[i+1], tr2[i+1] = f((tr1[i], tr2[i]), p)
    end
    tr1, tr2
end

# %%
tr1, tr2 = traj2d!(zeros(Int, 11), zeros(Int, 11), ((a, b), p) -> (b, a+b), nothing, (0, 1), 10, 10)
@show tr1 tr2;

# %%
function findperiod(x; tol=1e-3)
    @inbounds for k in 2:length(x)
        if abs(x[k] - x[1]) ≤ tol
            all(j -> abs(x[j] - x[j-k+1]) ≤ tol, k:length(x)) && return k - 1
        end
    end 
    return length(x)
end

# %%
x = @. sinpi(0:0.1:10) + 1e-4randn()
findperiod(x)

# %%
function calcperiod2d!(tr1, tr2, f, p, x0, lentr, nskips; tol=1e-3)
    traj2d!(tr1, tr2, f, p, x0, lentr, nskips)
    findperiod(tr1; tol)
end

# %%
using BenchmarkTools
@btime calcperiod2d!($(zeros(2001)), $(zeros(2001)), f, (2.6, 4.0), (0.4, 0.5), 2000, 50000; tol=1e-3)

# %%
function calcperiod2d(f; lentr = 2000, nskips = 50000, x0 = (0.4, 0.5),
        r = (1:0.05:5)[1:end-1], k = (2:0.05:5)[1:end-1], tol=1e-3
    )
    tr1 = Vector{Float64}(undef, lentr+1)
    tr2 = Vector{Float64}(undef, lentr+1)
    period = Matrix{Int}(undef, length(r), length(k))
    prog = Progress(length(r)*length(k))
    for j in eachindex(k)
        for i in eachindex(r)
            @inbounds period[i, j] = calcperiod2d!(tr1, tr2, f, (r[i], k[j]), x0, lentr, nskips; tol)
            next!(prog)
        end
    end
    period, r, k
end

period, r, k = @time calcperiod2d(f)
heatmap(r, k, log10.(period'); xlabel="r", ylabel="k", title="log10(period)")

# %%
function calcperiod2d_threads(f; lentr = 2000, nskips = 50000, x0 = (0.4, 0.5),
        r = (1:0.05:5)[1:end-1], k = (2:0.05:5)[1:end-1], tol=1e-3
    )
    tr1 = [Vector{Float64}(undef, lentr+1) for _ in 1:nthreads()]
    tr2 = [Vector{Float64}(undef, lentr+1) for _ in 1:nthreads()]
    period = Matrix{Int}(undef, length(r), length(k))
    # prog = Progress(length(r)*length(k))
    @threads for j in eachindex(k)
        tid = threadid()
        for i in eachindex(r)
            @inbounds period[i, j] = calcperiod2d!(tr1[tid], tr2[tid], f, (r[i], k[j]), x0, lentr, nskips; tol)
            # next!(prog)
        end
    end
    period, r, k
end

period_th, r, k = @time calcperiod2d_threads(f); flush(stdout)
period_th, r, k = @time calcperiod2d_threads(f); flush(stdout)
period_th, r, k = @time calcperiod2d_threads(f); flush(stdout)
heatmap(r, k, log10.(period_th'); xlabel="r", ylabel="k", title="log10(period)")

# %%
period == period_th

# %%
period

# %%
@btime calcperiod2d_threads(f; lentr = 2000, nskips = 50000, x0 = (0.4, 0.5),
    r = (1:0.05:5)[1:end-1], k = (2:0.05:5)[1:end-1], tol=1e-3)

# %%
@code_warntype calcperiod2d_threads(f; lentr = 2000, nskips = 50000, x0 = (0.4, 0.5),
    r = (1:0.05:5)[1:end-1], k = (2:0.05:5)[1:end-1], tol=1e-3)

# %%
