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
# There are a few bugs in the code of https://discourse.julialang.org/t/solving-difference-equation-part-2/67057.
#
# * The code `for k in 2:length(x)` in the function `seqper_new(x; tol=1e-3)` must be corrected as `for k in 2:(length(x) ÷ 2 + 1)`.
# * `nxblock * (ii - 1) + jj` must be `nyblock * (ii - 1) + jj`.
# * `(i - 1) * xpts + j` must be `(i - 1) * ypts + j`
# * The `Threads.@threads for` loop of the function `isoperiodic_test` is not thread-safe, so that the result changes with each execution.  Threads are conflicting at the line `set_parameter!(ds, [par1range[i], par2range[j]])`, where the same `ds` will be modified by multiple threads.
#
# Memory allocations are also unnecessarily too much.
#
# * The size of the array `sol_last` in memory is (NIter + 1) * xpts * ypts * 8 bytes = 2001 * 100 * 100 * 8 bytes = 153 MB, which is quite large, but it is discarded at the end of the function `isoperiodic_test`.
# * Each execution of the `DynamicalSystems.trajectory` function results in a wasted memory allocation, which is repeated tens of thousands of times.
#
# Other points that need improvement
#
# * The calculation results are saved to a file, and the results cannot be reused without going through the file.
# * The code is unnecessarily complicated.
# * You have created the complicated function `isoperiodic_test` to do most of all the work you want to do.
#
# Furthermore, it is not necessary to use `DynamicalSystems.jl` in this case.
#
# __My test implementation:__

# %%
using Base.Threads
using Plots
pyplot() # https://github.com/JuliaPlots/Plots.jl/issues/3560
using CSV
using DataFrames

function f(x, p)
   r, k = p
   a, mu, d = 5.0, 0.5, 0.2
   xnew1 = x[1]*exp(r*(1-x[1]/k)-x[2]/(a+x[1]^2))
   xnew2 = x[2]*exp(mu*x[1]/(a+x[1]^2)-d)
   xnew1, xnew2
end

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

function findperiod(x; tol=1e-3)
    n = length(x)
    @inbounds for k in 2:(n ÷ 2 + 1)
        if abs(x[k] - x[1]) ≤ tol
            all(j -> abs(x[j] - x[j-k+1]) ≤ tol, k:n) && return k - 1
        end
    end 
    return n
end

function calcperiod2d(f = f, niters = 2000, nwarmup = 50000, x0 = (0.4, 0.5),
        r = (1:0.05:5)[1:end-1], k = (2:0.05:5)[1:end-1], tol=1e-3
    )
    tmp1 = [Vector{Float64}(undef, niters+1) for _ in 1:nthreads()]
    tmp2 = [Vector{Float64}(undef, niters+1) for _ in 1:nthreads()]
    period = Matrix{Int}(undef, length(r), length(k))
    @threads for j in eachindex(k)
        tid = threadid()
        tr1, tr2 = tmp1[tid], tmp2[tid]
        for i in eachindex(r)
            p = (r[i], k[j])
            traj2d!(tr1, tr2, f, p, x0, niters, nwarmup)
            @inbounds period[i, j] = findperiod(tr1; tol)
        end
    end
    period, r, k
end

xm1ox(x) = (x - 1)/x # x minus 1 over x : {x ≥　1} -> {0 ≤ y < 1}
inv1my(y) = 1/(1 - y) # inverse of 1 minus 1 : {0 ≤ y < 1} -> {x ≥ 1}

function plot_period(r, k, period; kwargs...)
    cbtick = [1, 2, 3, 4, 5, 6, 10, 30]
    colorbar_ticks = (xm1ox.(cbtick), string.(cbtick))
    heatmap(r, k, xm1ox.(period)'; xlabel="r", ylabel="k", title="period",
        colorbar_ticks, kwargs...)
end

function save_period(r, k, period; fn = "result.csv")
    ii, jj = reim(complex.(axes(r, 1), axes(k, 1)'))
    rr, kk = reim(complex.(r, k'))
    df = DataFrame(i = vec(ii), j = vec(jj), r = vec(rr), k = vec(kk), period = vec(period))
    CSV.write(fn, df)
end

function load_period(fn = "result.csv")
    df = CSV.read(fn, DataFrame)
    imax = maximum(df.i)
    jmax = maximum(df.j)
    r = reshape(df.r, imax, jmax)[1:imax]
    k = reshape(df.k, imax, jmax)[1:imax:end]
    period = reshape(df.period, imax, jmax)
    period, r, k
end

# %%
using Logging
disable_logging(Logging.Warn)

# %%
period, r, k = @time calcperiod2d()
period, r, k = @time calcperiod2d()
period, r, k = @time calcperiod2d()
plot_period(r, k, period)

# %%
period, r, k = calcperiod2d()
fn = save_period(r, k, period; fn = "result.csv")
df = CSV.read(fn, DataFrame)

# %%
period, r, k = load_period(fn)
plot_period(r, k, period)

# %%
niters = 2000
nwarmup = 50000
x0 = (0.4, 0.5)
r = range(1, 5, length=100+1)[1:end-1]
k = range(2, 5, length=100+1)[1:end-1]
period, = @time calcperiod2d(f, niters, nwarmup, x0, r, k)
period, = @time calcperiod2d(f, niters, nwarmup, x0, r, k)
period, = @time calcperiod2d(f, niters, nwarmup, x0, r, k)
plot_period(r, k, period)

# %%
niters = 2000
nwarmup = 50000
x0 = (0.4, 0.5)
r = range(1, 5, length=200+1)[1:end-1]
k = range(2, 5, length=200+1)[1:end-1]
period, = @time calcperiod2d(f, niters, nwarmup, x0, r, k)
period, = @time calcperiod2d(f, niters, nwarmup, x0, r, k)
period, = @time calcperiod2d(f, niters, nwarmup, x0, r, k)
plot_period(r, k, period)

# %%
