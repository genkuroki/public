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

# %% [markdown]
# # Bayesian Poisson回帰
#
# * 黒木玄
# * 2020-09-16
#
# Bayesian Poisson回帰と最尤法(+ Wald検定で使う正規分布近似)によるPoisson回帰の結果が近似的に一致することを数値的に確認する.

# %% [markdown] toc=true
# <h1>目次<span class="tocSkip"></span></h1>
# <div class="toc"><ul class="toc-item"><li><span><a href="#データ生成" data-toc-modified-id="データ生成-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>データ生成</a></span></li><li><span><a href="#平坦事前分布の場合" data-toc-modified-id="平坦事前分布の場合-2"><span class="toc-item-num">2&nbsp;&nbsp;</span>平坦事前分布の場合</a></span></li><li><span><a href="#正規事前分布の場合" data-toc-modified-id="正規事前分布の場合-3"><span class="toc-item-num">3&nbsp;&nbsp;</span>正規事前分布の場合</a></span></li></ul></div>

# %%
using Distributions
using ForwardDiff
using LinearAlgebra
using Optim
using Random
using StatsPlots
default(fmt=:png, titlefontsize=10, plot_titlefontsize=10)
using Turing

# %%
# Poisson回帰の対数尤度函数の-1倍とその1階と2階の導函数達
# 1階の導函数はスコア統計量と呼ばれ, 2階の導函数の期待値はFisher情報量と呼ばれる.
# この場合に2階の導函数は確率変数の y を含まないので,
# 2階の導函数がそのままFisher情報量になる.

function negloglik_poissonreg(x, y, β)
    -sum(logpdf(Poisson(exp(β[1] + β[2]*x)), y) for (x, y) in zip(x, y))
end

function scorestat_poissonreg(x, y, β)
    ForwardDiff.gradient(β -> negloglik_poissonreg(x, y, β), β)
end

function scorestat_poissonreg_explicit(x, y, β)
    λ(i) = exp(β[1] + β[2]*x[i])
    A = sum(λ(i) - y[i] for i in eachindex(x, y))
    B = sum((λ(i) - y[i])*x[i] for i in eachindex(x, y))
    [A, B]
end

function fisherinfo_poissonreg(x, y, β)
    ForwardDiff.hessian(β -> negloglik_poissonreg(x, y, β), β)
end

function fisherinfo_poissonreg_explicit(x, y, β)
    λ(i) = exp(β[1] + β[2]*x[i])
    a = sum(λ(i) for i in eachindex(x))
    b = sum(λ(i)*x[i] for i in eachindex(x))
    c = sum(λ(i)*x[i]^2 for i in eachindex(x))
    Symmetric([a b; b c])
end

# %%
# 最尤推定 + Wald法で使う正規分布近似における
# 定数パラメータβと確率変数β̂の役割を交換したもの.
# 以下ではβ̂の側は定数とみなされており, βの側が確率変数とみなされている.
# これは本質的にBayes統計での事後分布の対数の2次近似になっている.

function mvnormal_approx_posterior_poissonreg(x, y)
    o = optimize(β->negloglik_poissonreg(x, y, β), zeros(2), LBFGS())
    β̂ = o.minimizer # maximum likelihood estimate
    Î = fisherinfo_poissonreg_explicit(x, y, β̂) # Fisher information matrix
    MvNormal(β̂, inv(Î)) # mvnormal approximation of the posterior of β
end

# %%
# data
Random.seed!(4649373)
n = 20
x = sort(rand(Uniform(0, 3), n))
y = @. rand(Poisson(exp(-0.2 + 1.0x)))

# maximum likelihood estimate
o = optimize(β->negloglik_poissonreg(x, y, β), zeros(2), LBFGS())
β̂ = o.minimizer

# %%
scorestat_poissonreg(x, y, β̂)

# %%
scorestat_poissonreg_explicit(x, y, β̂)

# %%
fisherinfo_poissonreg(x, y, β̂)

# %%
Î = fisherinfo_poissonreg_explicit(x, y, β̂)

# %%
@show β̂ inv(Î);

# %%
mvnormal_approx_posterior_poissonreg(x, y)

# %% [markdown]
# ## データ生成

# %%
Random.seed!(4649373)
n = 20
β₀, β₁ = -0.2, 1.0
x = sort(rand(Uniform(0, 3), n))
y = @. rand(Poisson(exp(β₀ + β₁*x)))
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
xlim, ylim = quantile.(Ref(B0), (0.0001, 0.9999)), quantile.(Ref(B1), (0.001, 0.999))
mvnormal_approx = mvnormal_approx_posterior_poissonreg(x, y)
m = 5000
MV = rand(mvnormal_approx, m)

P1 = scatter(B0[1:m], B1[1:m]; label="MCMC", ms=2, msw=0, alpha=0.3)
plot!(; xguide="β₀", yguide="β₁", xlim, ylim)
P2 = scatter(MV[1,:], MV[2,:]; label="mvnormal approx", ms=2, msw=0, alpha=0.3)
plot!(; xguide="β₀", yguide="β₁", xlim, ylim)
plot(P1, P2; size=(800, 400), plot_title="posterior of (β₀, β₁)  (n=$n, flat prior case)")

# %%
β̂, Σ̂ = params(mvnormal_approx)
normal_approx0 = Normal(β̂[1], √Σ̂[1,1])
stephist(B0; norm=true, label="MCMC")
plot!(normal_approx0; label="normal approx", ls=:dash, lw=1.5)
title!("posterior of β₀  (n=$n, flat prior case)")

# %%
β̂, Σ̂ = params(mvnormal_approx)
normal_approx1 = Normal(β̂[2], √Σ̂[2,2])
stephist(B1; norm=true, label="MCMC")
plot!(normal_approx1; label="normal approx", ls=:dash, lw=1.5)
title!("posterior of β₁  (n=$n, flat prior case)")

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
xlim, ylim = quantile.(Ref(B0), (0.0001, 0.9999)), quantile.(Ref(B1), (0.001, 0.999))
mvnormal_approx = mvnormal_approx_posterior_poissonreg(x, y)
m = 5000
MV = rand(mvnormal_approx, m)

P1 = scatter(B0[1:m], B1[1:m]; label="MCMC", ms=2, msw=0, alpha=0.3)
plot!(; xguide="β₀", yguide="β₁", xlim, ylim)
P2 = scatter(MV[1,:], MV[2,:]; label="mvnormal approx", ms=2, msw=0, alpha=0.3)
plot!(; xguide="β₀", yguide="β₁", xlim, ylim)
plot(P1, P2; size=(800, 400), plot_title="posterior of (β₀, β₁)  (n=$n, normal prior case)")

# %%
β̂, Σ̂ = params(mvnormal_approx)
normal_approx0 = Normal(β̂[1], √Σ̂[1,1])
stephist(B0; norm=true, label="MCMC")
plot!(normal_approx0; label="normal approx", ls=:dash, lw=1.5)
title!("posterior of β₀  (n=$n, normal prior case)")

# %%
β̂, Σ̂ = params(mvnormal_approx)
normal_approx1 = Normal(β̂[2], √Σ̂[2,2])
stephist(B1; norm=true, label="MCMC")
plot!(normal_approx1; label="normal approx", ls=:dash, lw=1.5)
title!("posterior of β₁  (n=$n, normal prior case)")

# %%
