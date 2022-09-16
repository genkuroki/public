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
#     display_name: Julia 1.8.1
#     language: julia
#     name: julia-1.8
# ---

# %% [markdown] toc=true
# <h1>目次<span class="tocSkip"></span></h1>
# <div class="toc"><ul class="toc-item"><li><span><a href="#データ生成" data-toc-modified-id="データ生成-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>データ生成</a></span></li><li><span><a href="#平坦事前分布の場合" data-toc-modified-id="平坦事前分布の場合-2"><span class="toc-item-num">2&nbsp;&nbsp;</span>平坦事前分布の場合</a></span></li><li><span><a href="#正規事前分布の場合" data-toc-modified-id="正規事前分布の場合-3"><span class="toc-item-num">3&nbsp;&nbsp;</span>正規事前分布の場合</a></span></li></ul></div>

# %% [markdown]
# # Bayesian Poisson回帰
#
# * 黒木玄
# * 2020-09-16
#
# Bayesian Poisson回帰と最尤法(+ Wald検定で使う正規分布近似)によるPoisson回帰の結果が近似的に一致することを数値的に確認する.

# %%
using Distributions
using Optim
using Random
Random.seed!(4649373)
using StaticArrays
using StatsPlots
default(fmt=:png, titlefontsize=10, plot_titlefontsize=10)
using Turing

# %%
# 最尤推定 + Wald法で使う正規分布近似のβとβ̂の役割を交換したもの

function mvnormal_approx_posterior(x, y)
    λ(β₀, β₁, xᵢ) = exp(β₀ + β₁*xᵢ)
    f(w) = -sum(logpdf(Poisson(λ(w[1], w[2], x)), y) for (x, y) in zip(x, y))
    o = optimize(f, zeros(2), LBFGS())
    β̂₀, β̂₁ = β̂ = o.minimizer # maximum likelihood estimate
    λ̂(i) = λ(β̂₀, β̂₁, x[i])
    â = sum(λ̂(i) for i in eachindex(x))
    b̂ = sum(λ̂(i)*x[i] for i in eachindex(x))
    ĉ = sum(λ̂(i)*x[i]^2 for i in eachindex(x))
    Î = SMatrix{2,2}(â, b̂, b̂, ĉ) # Fisher information at β̂
    MvNormal(β̂, inv(Î)) # mvnormal approximation of the posterior of β
end

# %% [markdown]
# ## データ生成

# %%
n = 20
x = sort(rand(Uniform(0, 3), n))
β₀, β₁ = -0.2, 1.0
pois = @. Poisson(exp(β₀ + β₁*x))
y = rand.(pois)
scatter(x, y; label="data (x, y)", msw=0, legend=:topleft)

# %% [markdown]
# ## 平坦事前分布の場合

# %%
@model function poissonreg(x, y)
    β₀ ~ Turing.Flat()
    β₁ ~ Turing.Flat()
    for i in eachindex(x, y)
        y[i] ~ Poisson(exp(β₀ + β₁*x[i]))
    end
end

# %%
chn = sample(poissonreg(x, y), NUTS(), MCMCThreads(), 10^5, 3);

# %%
chn

# %%
plot(chn[1001:5000]; lw=0.5, leftmargin=4Plots.mm, bottommargin=4Plots.mm)

# %%
B0, B1 = vec(chn[:β₀]), vec(chn[:β₁])

# %%
mvnormal_approx = mvnormal_approx_posterior(x, y)
m = 5000
MV = rand(mvnormal_approx, m)
P1 = scatter(B0[1:m], B1[1:m]; label="MCMC", ms=2, msw=0, alpha=0.3)
plot!(xguide="β₀", yguide="β₁")
P2 = scatter(MV[1,:], MV[2,:]; label="mvnormal approx", ms=2, msw=0, alpha=0.3)
plot!(xguide="β₀", yguide="β₁")
plot(P1, P2; size=(800, 400), plot_title="posterior of (β₀, β₁)")

# %%
β̂, Σ̂ = params(mvnormal_approx)
normal_approx0 = Normal(β̂[1], √Σ̂[1,1])
stephist(B0; norm=true, label="MCMC")
plot!(normal_approx0; label="normal approx", ls=:dash, lw=1.5)
title!("posterior of β₀")

# %%
β̂, Σ̂ = params(mvnormal_approx)
normal_approx1 = Normal(β̂[2], √Σ̂[2,2])
stephist(B1; norm=true, label="MCMC")
plot!(normal_approx1; label="normal approx", ls=:dash, lw=1.5)
title!("posterior of β₁")

# %% [markdown]
# ## 正規事前分布の場合

# %%
@model function poissonreg_normalprior(x, y)
    β₀ ~ Normal(0, 2)
    β₁ ~ Normal(0, 2)
    for i in eachindex(x, y)
        y[i] ~ Poisson(exp(β₀ + β₁*x[i]))
    end
end

# %%
chn_np = sample(poissonreg_normalprior(x, y), NUTS(), MCMCThreads(), 10^5, 3);

# %%
chn_np

# %%
plot(chn_np[1001:5000]; lw=0.5, leftmargin=4Plots.mm, bottommargin=4Plots.mm)

# %%
B0, B1 = vec(chn_np[:β₀]), vec(chn_np[:β₁])

# %%
mvnormal_approx = mvnormal_approx_posterior(x, y)
m = 5000
MV = rand(mvnormal_approx, m)
P1 = scatter(B0[1:m], B1[1:m]; label="MCMC", ms=2, msw=0, alpha=0.3)
plot!(xguide="β₀", yguide="β₁")
P2 = scatter(MV[1,:], MV[2,:]; label="mvnormal approx", ms=2, msw=0, alpha=0.3)
plot!(xguide="β₀", yguide="β₁")
plot(P1, P2; size=(800, 400), plot_title="posterior of (β₀, β₁)")

# %%
β̂, Σ̂ = params(mvnormal_approx)
normal_approx0 = Normal(β̂[1], √Σ̂[1,1])
stephist(B0; norm=true, label="MCMC")
plot!(normal_approx0; label="normal approx", ls=:dash, lw=1.5)
title!("posterior of β₀")

# %%
β̂, Σ̂ = params(mvnormal_approx)
normal_approx1 = Normal(β̂[2], √Σ̂[2,2])
stephist(B1; norm=true, label="MCMC")
plot!(normal_approx1; label="normal approx", ls=:dash, lw=1.5)
title!("posterior of β₁")

# %%
