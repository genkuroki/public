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
#     display_name: Julia 1.11.3
#     language: julia
#     name: julia-1.11
# ---

# %%
module O

using Distributions
using LinearAlgebra
using Random
using Printf

function pvalue_N(X)
    n = length(X)
    z = sum(X) / √n
    ccdf(Chisq(1), z^2)
end

function pvalue_F(X)
    n = length(X)
    χ² = dot(X, X)
    ccdf(Chisq(n), χ²)
end

struct SimPval{M,P,T}
    μ::M
    pval_N::P
    pval_F::P
    α::T
end

function Base.show(io::IO, ::MIME"text/plain", sp::SimPval)
    (; μ, pval_N, pval_F, α) = sp
    Nrej = pval_N .< α
    Nacc = .!Nrej
    Frej = pval_F .< α
    Facc = .!Frej
    if iszero(μ)
        print(io, @sprintf("alpha error rate for test of \"norm(μ)=0\": %5.1f%%\n", 100mean(Frej)))
    else
        print(io, @sprintf("power for test of \"norm(μ)=0\": %5.1f%%    (norm(μ) = %.3f)\n", 100mean(Frej), norm(μ)))
    end
    if iszero(mean(μ))
        print(io, @sprintf("alpha error rate for test of \"mean(μ)=0\": %5.1f%%\n", 100mean(Nrej)))
    else
        print(io, @sprintf("power for test of \"mean(μ)=0\": %5.1f%%    (mean(μ) = %.3f)\n", 100mean(Nrej), mean(μ)))
    end
    print(io, "\n")
    print(io, @sprintf("                       %18s  %22s\n", "reject \"mean(μ)=0\"", "not reject \"mean(μ)=0\""))
    print(io, @sprintf("reject     \"norm(μ)=0\" %11.1f%%  %18.1f%%\n", 100mean(Frej .& Nrej), 100mean(Frej .& Nacc)))
    print(io, @sprintf("not reject \"norm(μ)=0\" %11.1f%%  %18.1f%%\n", 100mean(Facc .& Nrej), 100mean(Facc .& Nacc)))
end

function sim_pval(; n=20, μ=zeros(n), niters=10^6, α=0.05)
    pval_N = zeros(niters)
    pval_F = zeros(niters)
    dist = product_distribution(Normal.(μ))
    Xtmp = [rand(dist) for _ in 1:Threads.nthreads()]
    Threads.@threads :static for i in 1:niters
        tid = Threads.threadid()
        X = rand!(dist, Xtmp[tid])
        pval_N[i] = pvalue_N(X)
        pval_F[i] = pvalue_F(X)
    end
    O.SimPval(μ, pval_N, pval_F, α)
end

end

# %% [markdown]
# モデル:
# $
# X_i \sim \mathrm{Normal}(\mu_i, 1) \quad \text{($i=1,2,\ldots,n$) かつ $X_i$達は独立}
# $
#
# 帰無仮説N:
# $
# \dfrac{\mu_1+\cdots+\mu_n}{n} = 0
# $
#
# 帰無仮説F:
# $
# (\mu_1,\ldots,\mu_n)=(0,\ldots,0)
# $

# %%
O.sim_pval(; n=5, niters=10^7, α=0.05)

# %%
O.sim_pval(; n=20, niters=10^7, α=0.05)

# %%
O.sim_pval(; n=100, niters=10^7, α=0.05)

# %%
O.sim_pval(; n=1000, niters=10^7, α=0.05)

# %% tags=[]
O.sim_pval(; μ=fill(0.5, 20), niters=10^7, α=0.05)

# %%
O.sim_pval(; μ=[fill(0.5, 10); fill(-0.5, 10)], niters=10^7, α=0.05)

# %%
O.sim_pval(; μ=[5; zeros(19)], niters=10^7, α=0.05)

# %%
