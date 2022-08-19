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
#     display_name: Julia 1.8.0
#     language: julia
#     name: julia-1.8
# ---

# %%
using Profile
using ProfileCanvas
using BenchmarkTools

# %%
using Distributions

function sim_t_test(; dist = Normal(), n = 10, μ = mean(dist), L = 10^5)
    pval = Float64[]
    for i in 1:L
        X = rand(dist, n)
        X̄ = mean(X)
        S = std(X)
        t = √n * (X̄ - μ)/S
        p = 2ccdf(TDist(n-1), abs(t))
        push!(pval, p)
    end
    pval
end

pval = @time sim_t_test(L=100)
pval = @time sim_t_test(L=100)
pval = @time sim_t_test(L=100)
pval = @btime sim_t_test(L=100);

# %%
@profview_allocs sim_t_test(L=100) sample_rate=1

# %%
using Distributions
using Random

function sim_t_test_rev(; dist = Normal(), n = 10, μ = mean(dist), L = 10^5,
        pval = Vector{Float64}(undef, L),
        tmp = Vector{Float64}(undef, n)
    )
    for i in 1:L
        X = rand!(dist, tmp)
        X̄ = mean(X)
        S = std(X)
        t = √n * (X̄ - μ)/S
        pval[i] = 2ccdf(TDist(n-1), abs(t))
    end
    pval
end

n = 10
L = 100
pval =  Vector{Float64}(undef, L)
tmp =  Vector{Float64}(undef, n)
pval = @time sim_t_test_rev(; n, L, pval, tmp)
pval = @time sim_t_test_rev(; n, L, pval, tmp)
pval = @time sim_t_test_rev(; n, L, pval, tmp)
pval = @btime sim_t_test_rev(; n=$n, L=$L, pval=$pval, tmp=$tmp);

# %%
@profview_allocs sim_t_test_rev(; n, L, pval, tmp) sample_rate=1

# %%
