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
#     display_name: Julia 1.6.1
#     language: julia
#     name: julia-1.6
# ---

# %% [markdown]
# # Getting even faster by reducing allocations and by multi-threading
#
# * Gen Kuroki
# * 2021-06-28

# %%
VERSION

# %% [markdown]
# ## Original version
#
# * https://twitter.com/elbersb/status/1409454534666207232
# * https://elbersb.com/public/posts/interaction_simulation/

# %%
using Distributions 
using DataFrames     

const std_normal = Normal(0, 1)

function gen_data()
    ## total time periods in the the panel = 500
    tt = 500

    # x1 and x2 covariates
    x1_A = 1 .+ rand(std_normal, tt)
    x1_B = 1/4 .+ rand(std_normal, tt)
    x2_A = 1 .+ x1_A .+ rand(std_normal, tt)
    x2_B = 1 .+ x1_B .+ rand(std_normal, tt)

    # outcomes (notice different slope coefs for x2_A and x2_B)
    y_A = x1_A .+ 1*x2_A + rand(std_normal, tt)
    y_B = x1_B .+ 2*x2_B + rand(std_normal, tt)

    # combine
    DataFrame(
        id = vcat(fill(0, length(x1_A)), fill(1, length(x1_B))),
        x1 = vcat(x1_A, x1_B),
        x2 = vcat(x2_A, x2_B),
        x1_dmean = vcat(x1_A .- mean(x1_A), x1_B .- mean(x1_B)),
        x2_dmean = vcat(x2_A .- mean(x2_A), x2_B .- mean(x2_B)),
        y = vcat(y_A, y_B))
end

# %%
using GLM

function coefs_lm_formula(data)
    mod_level = lm(@formula(y ~ id + x1 * x2), data)
    mod_dmean = lm(@formula(y ~ id + x1_dmean * x2_dmean), data)
    (coef(mod_level)[end], coef(mod_dmean)[end])
end

# example
data = gen_data()
coefs_lm_formula(data)

# %%
function run_simulations(nsim)
    sims = zeros(nsim, 2);
    for i in 1:nsim
        data = gen_data()
        sims[i, :] .= coefs_lm_formula(data)
    end
    sims
end

using BenchmarkTools
n = 20000
@btime run_simulations($n);

# %%
using Plots
sims = run_simulations(n)
histogram(sims, label = ["level" "dmean"])

# %%
function coefs_lm_formula(data)
    constant = fill(1, nrow(data))
    X = Float64[constant data.id data.x1 data.x2 data.x1 .* data.x2]
    mod_level = fit(LinearModel, X, data.y)
    X[:, 5] .= data.x1_dmean .* data.x2_dmean
    mod_dmean = fit(LinearModel, X, data.y)
    (coef(mod_level)[end], coef(mod_dmean)[end])
end

@btime run_simulations($n);

# %%
using LinearAlgebra 
fastfit(X, y) = cholesky!(Symmetric(X' * X)) \ (X' * y)

function coefs_lm_formula(data)
    constant = fill(1, nrow(data))
    X = Float64[constant data.id data.x1 data.x2 data.x1 .* data.x2]
    mod_level = fastfit(X, data.y)
    X[:, 5] .= data.x1_dmean .* data.x2_dmean
    mod_dmean = fastfit(X, data.y)
    (mod_level[end], mod_dmean[end])
end

@btime run_simulations($n);

# %% [markdown]
# ## Optimized version: Reducing allocations
#
# Ref. [Performance Tips: Pre-allocating outputs](https://docs.julialang.org/en/v1/manual/performance-tips/#Pre-allocating-outputs)

# %%
using DataFrames
using Distributions
using LinearAlgebra
using Random: rand!
using BenchmarkTools
using Plots

const std_normal = Normal(0, 1)

function undef_data(tt)
    id = vcat(fill(0, tt), fill(1, tt))
    x1 = similar(id, Float64)
    x2 = similar(x1)
    x1_dmean = similar(x1)
    x2_dmean = similar(x2)
    y = similar(x1)
    DataFrame(id = id, x1 = x1, x2 = x2, x1_dmean = x1_dmean, x2_dmean = x2_dmean, y = y)
end

function gen_data!(data, tmp)
    tt = length(tmp) # Assume 2tt = nrow(data)
    A, B = 1:tt, tt+1:2tt
    
    # Warning: dataframe is type-unstable
    x1::Vector{Float64} = data.x1
    x2::Vector{Float64} = data.x2
    x1_dmean::Vector{Float64} = data.x1_dmean
    x2_dmean::Vector{Float64} = data.x2_dmean
    y::Vector{Float64} = data.y
    
    @views x1[A] .= 1   .+ rand!(std_normal, tmp)
    @views x1[B] .= 1/4 .+ rand!(std_normal, tmp)
    @views x2[A] .= 1   .+ x1[A] .+ rand!(std_normal, tmp)
    @views x2[B] .= 1   .+ x1[B] .+ rand!(std_normal, tmp)

    # outcomes (notice different slope coefs for x2[A] and x2[B])
    @views y[A] .= x1[A] .+ 1 .* x2[A] .+ rand!(std_normal, tmp)
    @views y[B] .= x1[B] .+ 2 .* x2[B] .+ rand!(std_normal, tmp)
    
    @views x1_dmean[A] .= x1[A] .- mean(x1[A])
    @views x1_dmean[B] .= x1[B] .- mean(x1[B])
    @views x2_dmean[A] .= x2[A] .- mean(x2[A])
    @views x2_dmean[B] .= x2[B] .- mean(x2[B])
    
    data
end

function fastfit!(mod, X, y, XtX, Xty)
    mul!(XtX, X', X)
    mul!(Xty, X', y)
    ldiv!(mod, cholesky!(Symmetric(XtX)), Xty)
end

function coefs_lm_formula!(result, data, X, mod_level, mod_dmean, XtX, Xty)
    # Warning: dataframe is type-unstable.
    id::Vector{Int} = data.id
    x1::Vector{Float64} = data.x1
    x2::Vector{Float64} = data.x2
    x1_dmean::Vector{Float64} = data.x1_dmean
    x2_dmean::Vector{Float64} = data.x2_dmean
    y::Vector{Float64} = data.y
    
    X[:, 1] .= 1
    X[:, 2] .= id
    X[:, 3] .= x1
    X[:, 4] .= x2
    X[:, 5] .= x1 .* x2
    fastfit!(mod_level, X, y, XtX, Xty)

    X[:, 5] .= x1_dmean .* x2_dmean
    fastfit!(mod_dmean, X, y, XtX, Xty)
    
    result .= (mod_level[end], mod_dmean[end])
end

function run_simulations_optimized(nsim, tt=500)
    sims = Matrix{Float64}(undef, nsim, 2)
    data = undef_data(tt)
    tmp = Vector{Float64}(undef, tt)
    X = Matrix{Float64}(undef, 2tt, 5)
    mod_level = Vector{Float64}(undef, 5)
    mod_dmean = similar(mod_level)
    XtX = Matrix{Float64}(undef, 5, 5)
    Xty = Vector{Float64}(undef, 5)
    for i in 1:nsim
        gen_data!(data, tmp)
        @views coefs_lm_formula!(sims[i, :], data, X, mod_level, mod_dmean, XtX, Xty)
    end
    sims
end

@show n = 20000
@time sims = run_simulations_optimized(n)
histogram(sims, label = ["level" "dmean"])

# %%
tt = 500
data = undef_data(tt)
tmp = Vector{Float64}(undef, tt)
X = Matrix{Float64}(undef, 2tt, 5)
mod_level = Vector{Float64}(undef, 5)
mod_dmean = similar(mod_level)
XtX = Matrix{Float64}(undef, 5, 5)
Xty = Vector{Float64}(undef, 5)
sim = zeros(2)
n = 20000

println("Optimized version:")
@btime gen_data!($data, $tmp)
@btime coefs_lm_formula!($sim, $data, $X, $mod_level, $mod_dmean, $XtX, $Xty)
@btime run_simulations_optimized($n);

# %%
println("Original version:")
@btime gen_data()
@btime coefs_lm_formula($data)
@btime run_simulations($n);

# %% [markdown]
# ## Multi-thread original version

# %%
function run_simulations_multithread(nsim)
    sims = zeros(nsim, 2);
    Threads.@threads for i in 1:nsim
        data = gen_data()
        sims[i, :] .= coefs_lm_formula(data)
    end
    sims
end

@show Threads.nthreads()
@time sims = run_simulations_multithread(n)
histogram(sims, label = ["level" "dmean"])

# %%
println("Multi-thread original version: ")
@show Threads.nthreads()
@btime run_simulations_multithread($n);

# %% [markdown]
# ## Multi-thread optimized version

# %%
# pkg> dev https://github.com/genkuroki/MyUtils.jl
using MyUtils: @my_threads

function run_simulations_multithread_optimized(nsim, tt=500)
    sims = Matrix{Float64}(undef, nsim, 2)
    @my_threads begin
        data = undef_data(tt)
        tmp = Vector{Float64}(undef, tt)
        X = Matrix{Float64}(undef, 2tt, 5)
        mod_level = Vector{Float64}(undef, 5)
        mod_dmean = similar(mod_level)
        XtX = Matrix{Float64}(undef, 5, 5)
        Xty = Vector{Float64}(undef, 5)
    end for i in 1:nsim
        gen_data!(data, tmp)
        @views coefs_lm_formula!(sims[i, :], data, X, mod_level, mod_dmean, XtX, Xty)
    end begin
    end
    sims
end

@show Threads.nthreads()
@time sims = run_simulations_multithread_optimized(n)
histogram(sims, label = ["level" "dmean"])

# %%
println("Multi-thread optimized version: ")
@show Threads.nthreads()
@btime run_simulations_multithread_optimized($n);

# %% [markdown]
# ## non-DataFrame version
#
# https://discourse.julialang.org/t/simulation-of-regression-coefficients/63683/2

# %%
using Distributions
using LinearAlgebra
using BenchmarkTools
using Plots
using Base.Threads

const std_normal = Normal(0, 1)
fastfit(X, y) = cholesky!(Symmetric(X' * X)) \ (X' * y)

function gen_data2!(X, y, x1_dmean, x2_dmean)
    ## total time periods in the the panel = 500
    tt = 500

    # x1 and x2 covariates
    x1_A = 1 .+ rand(std_normal, tt)
    x1_B = 1/4 .+ rand(std_normal, tt)
    x2_A = 1 .+ x1_A .+ rand(std_normal, tt)
    x2_B = 1 .+ x1_B .+ rand(std_normal, tt)

    # outcomes (notice different slope coefs for x2_A and x2_B)
    y_A = x1_A .+ 1*x2_A + rand(std_normal, tt)
    y_B = x1_B .+ 2*x2_B + rand(std_normal, tt)
    x1 = vcat(x1_A, x1_B)
    x2 = vcat(x2_A, x2_B)

    X[:, 2] .= vcat(fill(0, length(x1_A)), fill(1, length(x1_B)))
    X[:, 3] .= x1
    X[:, 4] .= x2
    X[:, 5] .= (x1 .* x2)

    y .= vcat(y_A, y_B)

    x1_dmean .= vcat(x1_A .- mean(x1_A), x1_B .- mean(x1_B))
    x2_dmean .= vcat(x2_A .- mean(x2_A), x2_B .- mean(x2_B))
end

function coefs_lm_formula4(X, y, x1_dmean, x2_dmean)
    mod_level = fastfit(X, y)
    @inbounds X[:, 5] .= x1_dmean .* x2_dmean
    mod_dmean = fastfit(X, y)
    (mod_level[end], mod_dmean[end])
end

function run_simulations5(nsim)
    sims = zeros(nsim, 2);
    X = [zeros(Float64, 1000, 5) for _ in 1:nthreads()]
    y = [zeros(Float64, 1000) for _ in 1:nthreads()]
    x1_dmean = [zeros(Float64, 1000) for _ in 1:nthreads()]
    x2_dmean = [zeros(Float64, 1000) for _ in 1:nthreads()]

    # set constant
    for i in 1:nthreads()
        X[i][:, 1] .= 1
    end

    @threads for i in 1:nsim
        gen_data2!(X[threadid()], y[threadid()], x1_dmean[threadid()], x2_dmean[threadid()])
        @inbounds sims[i, 1], sims[i, 2] = coefs_lm_formula4(X[threadid()], y[threadid()], x1_dmean[threadid()], x2_dmean[threadid()])
    end
    sims
end

# %%
n = 20000
@show Threads.nthreads()
@time sims = run_simulations5(n)
histogram(sims, label = ["level" "dmean"])

# %%
println("Multi-thread non-dataframe version (run_simulations5): ")
@show Threads.nthreads()
@btime run_simulations5($n);

# %% [markdown]
# ## Optimized non-DataFrame version

# %%
using Distributions
using LinearAlgebra
using Random: rand!
using BenchmarkTools
using Plots

const std_normal = Normal(0, 1)

# use NamedTuple indtead of DataFrame
function undef_data_namedtuple(tt)
    id = vcat(fill(0, tt), fill(1, tt))
    x1 = similar(id, Float64)
    x2 = similar(x1)
    x1_dmean = similar(x1)
    x2_dmean = similar(x2)
    y = similar(x1)
    (; id, x1, x2, x1_dmean, x2_dmean, y)
end

function gen_data!(data, tmp)
    tt = length(tmp) # Assume 2tt = nrow(data)
    A, B = 1:tt, tt+1:2tt
    
    # Warning: dataframe is type-unstable
    x1::Vector{Float64} = data.x1
    x2::Vector{Float64} = data.x2
    x1_dmean::Vector{Float64} = data.x1_dmean
    x2_dmean::Vector{Float64} = data.x2_dmean
    y::Vector{Float64} = data.y
    
    @views x1[A] .= 1   .+ rand!(std_normal, tmp)
    @views x1[B] .= 1/4 .+ rand!(std_normal, tmp)
    @views x2[A] .= 1   .+ x1[A] .+ rand!(std_normal, tmp)
    @views x2[B] .= 1   .+ x1[B] .+ rand!(std_normal, tmp)

    # outcomes (notice different slope coefs for x2[A] and x2[B])
    @views y[A] .= x1[A] .+ 1 .* x2[A] .+ rand!(std_normal, tmp)
    @views y[B] .= x1[B] .+ 2 .* x2[B] .+ rand!(std_normal, tmp)
    
    @views x1_dmean[A] .= x1[A] .- mean(x1[A])
    @views x1_dmean[B] .= x1[B] .- mean(x1[B])
    @views x2_dmean[A] .= x2[A] .- mean(x2[A])
    @views x2_dmean[B] .= x2[B] .- mean(x2[B])
    
    data
end

function fastfit!(mod, X, y, XtX, Xty)
    mul!(XtX, X', X)
    mul!(Xty, X', y)
    ldiv!(mod, cholesky!(Symmetric(XtX)), Xty)
end

function coefs_lm_formula!(result, data, X, mod_level, mod_dmean, XtX, Xty)
    # Warning: dataframe is type-unstable.
    id::Vector{Int} = data.id
    x1::Vector{Float64} = data.x1
    x2::Vector{Float64} = data.x2
    x1_dmean::Vector{Float64} = data.x1_dmean
    x2_dmean::Vector{Float64} = data.x2_dmean
    y::Vector{Float64} = data.y
    
    X[:, 1] .= 1
    X[:, 2] .= id
    X[:, 3] .= x1
    X[:, 4] .= x2
    X[:, 5] .= x1 .* x2
    fastfit!(mod_level, X, y, XtX, Xty)

    X[:, 5] .= x1_dmean .* x2_dmean
    fastfit!(mod_dmean, X, y, XtX, Xty)
    
    result .= (mod_level[end], mod_dmean[end])
end

function run_simulations_optimized_nondataframe(nsim, tt=500)
    sims = Matrix{Float64}(undef, nsim, 2)
    data = undef_data_namedtuple(tt) # use NamedTuple indtead of DataFrame
    tmp = Vector{Float64}(undef, tt)
    X = Matrix{Float64}(undef, 2tt, 5)
    mod_level = Vector{Float64}(undef, 5)
    mod_dmean = similar(mod_level)
    XtX = Matrix{Float64}(undef, 5, 5)
    Xty = Vector{Float64}(undef, 5)
    for i in 1:nsim
        gen_data!(data, tmp)
        @views coefs_lm_formula!(sims[i, :], data, X, mod_level, mod_dmean, XtX, Xty)
    end
    sims
end

@show n = 20000
@time sims = run_simulations_optimized_nondataframe(n)
histogram(sims, label = ["level" "dmean"])

# %%
println("Optimized non-dataframe version (single thread): ")
@btime run_simulations_optimized_nondataframe($n);

# %% [markdown]
# ## Multi-thread optimized non-DataFrame version

# %%
# pkg> dev https://github.com/genkuroki/MyUtils.jl
using MyUtils: @my_threads

function run_simulations_multithread_optimized_nondataframe(nsim, tt=500)
    sims = Matrix{Float64}(undef, nsim, 2)
    @my_threads begin
        data = undef_data_namedtuple(tt)
        tmp = Vector{Float64}(undef, tt)
        X = Matrix{Float64}(undef, 2tt, 5)
        mod_level = Vector{Float64}(undef, 5)
        mod_dmean = similar(mod_level)
        XtX = Matrix{Float64}(undef, 5, 5)
        Xty = Vector{Float64}(undef, 5)
    end for i in 1:nsim
        gen_data!(data, tmp)
        @views coefs_lm_formula!(sims[i, :], data, X, mod_level, mod_dmean, XtX, Xty)
    end begin
    end
    sims
end

@show Threads.nthreads()
@time sims = run_simulations_multithread_optimized_nondataframe(n)
histogram(sims, label = ["level" "dmean"])

# %%
println("Multi-thread optimized non-dataframe version: ")
@show Threads.nthreads()
@btime run_simulations_multithread_optimized_nondataframe($n);

# %% [markdown]
# ## Summary

# %%
n = 20000

@show VERSION
@show n
println()

println("Original version (with DataFrame, single thread):")
@btime run_simulations($n)
println()

println("Optimized version (with DataFrame, single thread):")
@btime run_simulations_optimized($n)
println()

println("Multi-thread optimized version (with DataFrame): ")
@show Threads.nthreads()
@btime run_simulations_multithread_optimized($n)
println()

println("Multi-thread non-dataframe version (run_simulations5): ")
@show Threads.nthreads()
@btime run_simulations5($n)
println()

println("Optimized non-dataframe version (single thread): ")
@btime run_simulations_optimized_nondataframe($n)
println()

println("Multi-thread optimized non-dataframe version: ")
@show Threads.nthreads()
@btime run_simulations_multithread_optimized_nondataframe($n)
println()

# %%
versioninfo()

# %%
