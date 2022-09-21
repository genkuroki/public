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
# # 2つの正規分布の直積モデルのベイズ統計によるWelchのt検定の再現
#
# * 黒木玄
# * 2022-09-16, 2022-09-22
#
# 2つの正規分布の直積モデルのベイズ統計でWelchのt検定が近似的に再現されることを確認する.

# %% [markdown] toc=true
# <h1>目次<span class="tocSkip"></span></h1>
# <div class="toc"><ul class="toc-item"><li><span><a href="#分散が異なる２つの正規分布でデータを生成した場合" data-toc-modified-id="分散が異なる２つの正規分布でデータを生成した場合-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>分散が異なる２つの正規分布でデータを生成した場合</a></span></li><li><span><a href="#異なる２つのガンマ分布でデータを生成した場合" data-toc-modified-id="異なる２つのガンマ分布でデータを生成した場合-2"><span class="toc-item-num">2&nbsp;&nbsp;</span>異なる２つのガンマ分布でデータを生成した場合</a></span></li><li><span><a href="#Welchのt検定での自由度の分布" data-toc-modified-id="Welchのt検定での自由度の分布-3"><span class="toc-item-num">3&nbsp;&nbsp;</span>Welchのt検定での自由度の分布</a></span></li></ul></div>

# %%
using Distributions
using QuadGK
using Random
Random.seed!(4649373)
using StatsBase: ecdf
using StatsPlots
default(fmt=:png, titlefontsize=10, plot_titlefontsize=10,
    tickfontsize=6, guidefontsize=9)

distname(dist) = replace(string(dist), r"{[^}]*}"=>"")

# %%
# 通常のWelchのt検定

function tvalue_welch(m, x̄, sx², n, ȳ, sy²; Δμ=0)
    (x̄ - ȳ - Δμ) / √(sx²/m + sy²/n)
end

function tvalue_welch(x, y; Δμ=0)
    m, x̄, sx² = length(x), mean(x), var(x)
    n, ȳ, sy² = length(y), mean(y), var(y)
    tvalue_welch(m, x̄, sx², n, ȳ, sy²; Δμ)
end

function degree_of_freedom_welch(m, sx², n, sy²)
    (sx²/m + sy²/n)^2 / ((sx²/m)^2/(m-1) + (sy²/n)^2/(n-1))
end

function degree_of_freedom_welch(x, y)
    m, sx² = length(x), var(x)
    n, sy² = length(y), var(y)
    degree_of_freedom_welch(m, sx², n, sy²)
end

function pvalue_welch(m, x̄, sx², n, ȳ, sy²; Δμ=0)
    t = tvalue_welch(m, x̄, sx², n, ȳ, sy²; Δμ)
    ν = degree_of_freedom_welch(m, sx², n, sy²)
    2ccdf(TDist(ν), abs(t))
end

function pvalue_welch(x, y; Δμ=0)
    m, x̄, sx² = length(x), mean(x), var(x)
    n, ȳ, sy² = length(y), mean(y), var(y)
    pvalue_welch(m, x̄, sx², n, ȳ, sy²; Δμ)
end

function confint_welch(m, x̄, sx², n, ȳ, sy²; α=0.05)
    ν = degree_of_freedom_welch(m, sx², n, sy²)
    c = quantile(TDist(ν), 1-α/2)
    SEhat = √(sx²/m + sy²/n)
    [x̄-ȳ-c*SEhat, x̄-ȳ+c*SEhat]
end

function confint_welch(x, y; α=0.05)
    m, x̄, sx² = length(x), mean(x), var(x)
    n, ȳ, sy² = length(y), mean(y), var(y)
    confint_welch(m, x̄, sx², n, ȳ, sy²; α)
end

# %%
# Welchのt検定を正確に再現する Δμ に関する分布

function posterior_welch(m, x̄, sx², n, ȳ, sy²)
    ν = degree_of_freedom_welch(m, sx², n, sy²)
    SEhat = √(sx²/m + sy²/n)
    x̄ - ȳ + SEhat*TDist(ν)
end

function posterior_welch(x, y)
    m, x̄, sx² = length(x), mean(x), var(x)
    n, ȳ, sy² = length(y), mean(y), var(y)
    posterior_welch(m, x̄, sx², n, ȳ, sy²)
end

# %% [markdown]
# 次のセルの内容に関する詳しい説明については
#
# * https://github.com/genkuroki/public/blob/main/0036/InverseGammaNormal%20conjugate%20prior%20of%20normal%20distribution.ipynb
# * https://github.com/genkuroki/public/blob/main/0036/InverseGammaNormal%20conjugate%20prior%20of%20normal%20distribution.pdf
#
# を参照せよ.

# %%
# 1つの正規分布モデルの共役事前分布のベイズ更新

function bayesian_update(μ₀, v₀, κ, θ, n, x̄, σ̂²)
    μ₀_new = (μ₀/v₀ + n*x̄)/(1/v₀ + n)
    v₀_new = 1/(1/v₀ + n)
    κ_new = κ + n/2
    θ_new = θ + (n/2)*(σ̂² + ((x̄ - μ₀)^2/v₀)/(1/v₀ + n))
    μ₀_new, v₀_new, κ_new, θ_new
end

function bayesian_update(μ₀, v₀, κ, θ, y)
    n, ȳ, σ̂² = length(y), mean(y), var(y; corrected=false)
    bayesian_update(μ₀, v₀, κ, θ, n, ȳ, σ̂²)
end

# μに関する事前事後分布
posterior_μ(μ₀, v₀, κ, θ) = μ₀ + √(θ/κ*v₀)*TDist(2κ)

# 事前予測分布
preddist(μ₀, v₀, κ, θ) = μ₀ + √(θ/κ*(1 + v₀))*TDist(2κ)

# %%
# (μx, μy) に関する事後分布

function posterior_μxμy(m, x̄, σ̂x², n, ȳ, σ̂y²;
        μx₀=0, vx₀=Inf, κx=0, θx=0, μy₀=0, vy₀=Inf, κy=0, θy=0)
    μx₀_new, vx₀_new, κx_new, θx_new = bayesian_update(μx₀, vx₀, κx, θx, m, x̄, σ̂x²)
    μy₀_new, vy₀_new, κy_new, θy_new = bayesian_update(μy₀, vy₀, κy, θy, n, ȳ, σ̂y²)
    posterior_μx = posterior_μ(μx₀_new, vx₀_new, κx_new, θx_new)
    posterior_μy = posterior_μ(μy₀_new, vy₀_new, κy_new, θy_new)
    posterior_μx, posterior_μy
end

function posterior_μxμy(x, y;
        μx₀=0, vx₀=Inf, κx=0, θx=0, μy₀=0, vy₀=Inf, κy=0, θy=0)
    m, x̄, σ̂x² = length(x), mean(x), var(x; corrected=false)
    n, ȳ, σ̂y² = length(y), mean(y), var(y; corrected=false)
    posterior_μxμy(m, x̄, σ̂x², n, ȳ, σ̂y²; μx₀, vx₀, κx, θx, μy₀, vy₀, κy, θy)
end

function rand_posterior_Δμ(x, y; L=10^6,
        μx₀=0, vx₀=Inf, κx=0, θx=0, μy₀=0, vy₀=Inf, κy=0, θy=0)
    posterior_μx, posterior_μy =
        posterior_μxμy(x, y; μx₀, vx₀, κx, θx, μy₀, vy₀, κy, θy)
    [rand(posterior_μx) - rand(posterior_μy) for _ in 1:L]
end

# %%
# X ~ distx, Y ~ disty のときの X - Y が従う分布の密度函数

function pdf_diff(distx, disty, Δ)
    K(y) = pdf(distx, y+Δ) * pdf(disty, y)
    quadgk(K, extrema(disty)...)[1]
end

# %%
function plot_bayesian_welch(x, y; L=10^6, title="m=$(length(x)), n=$(length(y))")
    posterior_μx, posterior_μy = posterior_μxμy(x, y)
    pdf_Δμ(Δμ) = pdf_diff(posterior_μx, posterior_μy, Δμ)
    post_welch = posterior_welch(x, y)
    xlim = quantile.(post_welch, (0.0001, 0.9999))
    plot(legend=:outertop)
    plot!(pdf_Δμ, xlim...; norm=true, label="posterior of Δμ w.r.t. the flat prior")
    plot!(post_welch, xlim...; label="\"posterior\" reproducing Welch t-test", ls=:dash)
    plot!(xguide="Δμ")
    plot!(; title, xlim)
end

function plot_bayesian_welch(distx, m, disty, n; L=10^6,
        x = rand(distx, m), y = rand(disty, n))
    @show mean(distx) - mean(disty)
    @show var(distx) var(disty)
    @show mean(x) - mean(y)
    @show var(x) var(y)
    @show degree_of_freedom_welch(x, y)
    @show min(m,n) - 1, m + n - 2
    title="$(distname(distx)), m=$m\n$(distname(disty)), n=$n"
    plot_bayesian_welch(x, y; L, title)
end

# %% [markdown]
# ## 分散が異なる２つの正規分布でデータを生成した場合

# %%
distx = Normal(0, 4)
disty = Normal(0, 1)
plot(distx; label="distx")
plot!(disty; label="disty", ls=:dash)

# %%
plot_bayesian_welch(distx, 7, disty, 7)

# %%
plot_bayesian_welch(distx, 10, disty, 5)

# %%
plot_bayesian_welch(distx, 5, disty, 10)

# %%
plot_bayesian_welch(distx, 15, disty, 15)

# %%
plot_bayesian_welch(distx, 20, disty, 10)

# %%
plot_bayesian_welch(distx, 10, disty, 20)

# %%
plot_bayesian_welch(distx, 30, disty, 30)

# %%
plot_bayesian_welch(distx, 40, disty, 20)

# %%
plot_bayesian_welch(distx, 20, disty, 40)

# %%
plot_bayesian_welch(distx, 60, disty, 60)

# %%
plot_bayesian_welch(distx, 80, disty, 40)

# %%
plot_bayesian_welch(distx, 40, disty, 80)

# %% [markdown]
# ## 異なる２つのガンマ分布でデータを生成した場合

# %%
distx = Gamma(5, 3)
disty = Gamma(3, 1)
plot(distx; label="distx")
plot!(disty; label="disty", ls=:dash)

# %%
plot_bayesian_welch(distx, 7, disty, 7)

# %%
plot_bayesian_welch(distx, 10, disty, 5)

# %%
plot_bayesian_welch(distx, 5, disty, 10)

# %%
plot_bayesian_welch(distx, 15, disty, 15)

# %%
plot_bayesian_welch(distx, 20, disty, 10)

# %%
plot_bayesian_welch(distx, 10, disty, 20)

# %%
plot_bayesian_welch(distx, 30, disty, 30)

# %%
plot_bayesian_welch(distx, 40, disty, 20)

# %%
plot_bayesian_welch(distx, 20, disty, 40)

# %%
plot_bayesian_welch(distx, 60, disty, 60)

# %%
plot_bayesian_welch(distx, 80, disty, 40)

# %%
plot_bayesian_welch(distx, 40, disty, 80)

# %% [markdown]
# ## Welchのt検定での自由度の分布

# %%
using Random: rand!

function plot_dist_df_welch(distx, m, disty, n; L=10^6,
        title="m=$m, n=$n", kwargs...)
    tmpx = [Vector{Float64}(undef, m) for _ in 1:Threads.nthreads()]
    tmpy = [Vector{Float64}(undef, n) for _ in 1:Threads.nthreads()]
    df = Vector{Float64}(undef, L)
    Threads.@threads for i in 1:L
        x = rand!(distx, tmpx[Threads.threadid()])
        y = rand!(distx, tmpy[Threads.threadid()])
        df[i] = degree_of_freedom_welch(x, y)
    end
    
    stephist(df; norm=true, label="df")
    plot!(xlim = (min(m,n)-1-0.5, m+n-2+0.5))
    plot!(legend=:topleft)
    title!(title)
    plot!(; kwargs...)
end

function plot_dists_df_welch(distx, disty, mns; L=10^6, kwargs...)
    @assert length(mns) == 6
    PP = []
    for (m, n) in mns
        P = plot_dist_df_welch(distx, m, disty, n; L, kwargs...)
        push!(PP, P)
    end
    plot(PP...; layout=(2, 3), size=(800, 400))
    plot!(plot_title="distx=$(distname(distx)),  disty=$(distname(disty))")
end

# %%
mns = [(10-i, 5+i) for i in 0:5]
plot_dists_df_welch(Normal(0,4), Normal(0,1), mns)

# %%
mns = [(20-i, 10+i) for i in 0:2:10]
plot_dists_df_welch(Normal(0,4), Normal(0,1), mns)

# %%
mns = [(40-i, 20+i) for i in 0:4:20]
plot_dists_df_welch(Normal(0,4), Normal(0,1), mns)

# %%
mns = [(80-i, 40+i) for i in 0:8:40]
plot_dists_df_welch(Normal(0,4), Normal(0,1), mns)

# %%
mns = [(10-i, 5+i) for i in 0:5]
plot_dists_df_welch(Gamma(5,3), Gamma(3,1), mns)

# %%
mns = [(20-i, 10+i) for i in 0:2:10]
plot_dists_df_welch(Gamma(5,3), Gamma(3,1), mns)

# %%
mns = [(40-i, 20+i) for i in 0:4:20]
plot_dists_df_welch(Gamma(5,3), Gamma(3,1), mns)

# %%
mns = [(80-i, 40+i) for i in 0:8:40]
plot_dists_df_welch(Gamma(5,3), Gamma(3,1), mns)

# %%
