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

struct SimPval{P,T}
    pval_N::P
    pval_F::P
    α::T
end

function Base.show(io::IO, ::MIME"text/plain", x::SimPval)
    (; pval_N, pval_F, α) = x
    Nrej = pval_N .< α
    Nacc = .!Nrej
    Frej = pval_F .< α
    Facc = .!Frej
    print(io, @sprintf("alpha error rate for N: %6.3f\n", mean(Nrej)))
    print(io, @sprintf("alpha error rate for F: %6.3f\n", mean(Frej)))
    print(io, "\n")
    print(io, @sprintf("          %8s  %8s\n", "reject F", "accept F"))
    print(io, @sprintf("reject N %8.3f  %8.3f\n", mean(Nrej .& Frej), mean(Nrej .& Facc)))
    print(io, @sprintf("accept N %8.3f  %8.3f\n", mean(Nacc .& Frej), mean(Nacc .& Facc)))
end

function sim_pval(; n=20, niters=10^6, α=0.05)
    pval_N = zeros(niters)
    pval_F = zeros(niters)
    Xtmp = [randn(n) for _ in 1:Threads.nthreads()]
    Threads.@threads :static for i in 1:niters
        tid = Threads.threadid()
        X = rand!(Normal(), Xtmp[tid])
        pval_N[i] = pvalue_N(X)
        pval_F[i] = pvalue_F(X)
    end
    O.SimPval(pval_N, pval_F, α)
end

end

# %% [markdown]
# モデル
# $$
# \mathrm{Normal}(\mu_1, 1)\times\cdots\times\mathrm{Normal}(\mu_n, 1)
# $$
#
# 帰無仮説N:
# $$
# \frac{\mu_1+\cdots+\mu_n}{n} = 0
# $$
#
# 帰無仮説F:
# $$
# (\mu_1,\ldots,\mu_n)=(0,\ldots,0)
# $$

# %%
O.sim_pval(; n=5, niters=10^7, α=0.05)

# %%
O.sim_pval(; n=20, niters=10^7, α=0.05)

# %%
O.sim_pval(; n=100, niters=10^7, α=0.05)

# %%
O.sim_pval(; n=1000, niters=10^7, α=0.05)

# %%
