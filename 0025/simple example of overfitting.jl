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
#     display_name: Julia 1.7.0
#     language: julia
#     name: julia-1.7
# ---

# %%
using LinearAlgebra
using Distributions
using Plots
using Random
using Optim
using QuadGK

"""最小二乗法による多項式フィッティング"""
function polynomialfit(x, y, d)
    n = length(y)   # データサイズ
    X = x .^ (0:d)' # design matrix
    β̂ = X \ y       # Juliaでは X \ y の一発で最小二乗法が可能
    f(x) = evalpoly(x, β̂)
    ŷ = f.(x)
    r = ŷ - y
    σ̂ = norm(r)/√n  # モデルの標準偏差パラメータの推定
    dist = product_distribution([Normal(μ, σ̂) for μ in ŷ])
    loglik = logpdf(dist, y)  # 最大化された対数尤度
    AIC = -2loglik + 2(d + 2) # 所謂赤池情報量規準
    (; x, y, d, β̂, σ̂, f, loglik, AIC)
end

"""optimize函数を使った最尤法による多項式フィッティング"""
function polynomialfit_optim(x, y, d)
    negloglik(w) = begin
        β = @view w[1:end-1]
        σ = exp(w[end])
        n = length(y)
        dist = product_distribution([Normal(evalpoly(x[i], β), σ) for i in 1:n])
        -loglikelihood(dist, y)
    end
    o = optimize(negloglik, zeros(d+2), Optim.Options(store_trace=true, extended_trace = true, iterations = 10^4))
    ŵ = o.minimizer
    β̂ = ŵ[1:end-1]
    σ̂ = exp(ŵ[end])
    loglik = -o.minimum
    AIC = 2o.minimum + 2(d + 2)
    f(x) = evalpoly(x, β̂)
    (; x, y, d, β̂, σ̂, f, loglik, AIC, negloglik, o)
end

function plotpolyfit(result; legend=:topleft, titlefontsize=10, kwargs...)
    (; x, y, d, f, loglik, AIC) = result
    xs = range(extrema(x)...; length=1000)
    plot(; legend, titlefontsize)
    scatter!(x, y; label="data")
    plot!(xs, f.(xs); label="degree-$d fitting")
    title!("degree = $d,  log likelihood = $(round(loglik; digits=2)),  AIC = $(round(AIC; digits=2))")
    plot!(; kwargs...)
end

# %%
# データの生成

rng = MersenneTwister(1234567890)

n = 21
x = range(-1, 1; length=n)
σ₀ = 0.3
f₀(x) = sinpi(x)
noise = σ₀*randn(rng, length(x))
y = f₀.(x) + noise

xs = range(extrema(x)...; length=1000)
P = plot(; legend=:topleft)
scatter!(x, y; label="data")
plot!(xs, sinpi.(xs); label="sin(πx)")
title!("data and true law"; titlefontsize=10)
plot!(; size=(400, 300))

# %%
β₀ = [isodd(k) ? (-1)^(k÷2)*π^k/factorial(k) : 0.0 for k in 0:3]
@show β₀ σ₀;

# %%
(; β̂, σ̂, loglik, AIC) = polynomialfit(x, y, 3)
@show β̂ σ̂ loglik AIC;

# %%
(; β̂, σ̂, loglik, AIC) = polynomialfit_optim(x, y, 3)
@show β̂ σ̂ loglik AIC;

# %%
Q = plotpolyfit(polynomialfit(x, y, 3))
R = plotpolyfit(polynomialfit(x, y, 9))
S = plotpolyfit(polynomialfit(x, y, 15))

plot(P, Q, R, S; size=(800, 600), layout=(2, 2), titlefontsize=9)

# %% [markdown]
# 対数尤度(log likelihood)の大きさは「もっともらしさ」の度合いではない！
#
# 対数尤度の大きさはモデルのデータへのフィッティングの良さの指標の1つに過ぎない。
#
# 上の結果を見れば分かるように、フィッティングで使う多項式の次数を上げた方がフィッティングの度合いが上がって、対数尤度も大きくなるが、次数が大きな場合の結果は全然「もっともらしい」とは言えない。
#
# 尤度を最大化することによってモデルのパラメータ値を調節する最尤法は「最も尤もらしいパラメータ値」を探すのではなく、単に「モデルをデータに最もフィットさせるようなパラメータ値」を探しているだけである。
#
# 最尤法は「最も尤もらしいパラメータ値」を探す方法ではない！

# %%
result = polynomialfit.(Ref(x), Ref(y), 0:19)
d = getproperty.(result, :d)
loglik = getproperty.(result, :loglik)
AIC = getproperty.(result, :AIC)
dbest = d[argmin(AIC)]

T = plot(d, loglik; label="log likelihood", legend=:topleft)
plot!(; xtick=d, xlabel="degree", ylabel="log likelihood")
title!("log likelihoods of polynomial fitting"; titlefontsize=10)

U = plot(d, AIC; label="AIC", legend=:topright)
plot!(; xtick=d, xlabel="degree", ylabel="AIC")
title!("AICs of polynomial fitting (minimized at degree $dbest)"; titlefontsize=10)

plot(T, U; size=(800, 300), titlefontsize=9, leftmargin=3Plots.mm, bottommargin=3Plots.mm)

# %%
function animate_polyfit(x, y, d, f₀, σ₀; gifname="polyfit$d.gif", fps=20, ylim=(-8, 28))
    (; x, y, d, negloglik, o) = polynomialfit_optim(x, y, d)
    a, b = extrema(x)
    n = length(y)
    value = -getproperty.(o.trace, :value)
    centroid = getindex.(getproperty.(o.trace, :metadata), "centroid")
    
    ngenerr = Vector{Float64}(undef, length(centroid))
    Threads.@threads for t in keys(centroid)
        β = centroid[t][1:end-1]
        f(x) = evalpoly(x, β)
        σ = exp(centroid[t][end])
        GE = (1/2)*(log(2π*σ^2) + σ₀^2/σ^2 + 1/(b-a)*quadgk(x -> (f(x) - f₀(x))^2, a, b)[1])
        ngenerr[t] = n*GE
    end

    xs = range(extrema(x)...; length=1000)
    ylim = isnothing(ylim) ? extrema([-value; ngenerr]) .+ (-1, 1) : ylim
    
    L = length(value)
    tstep = max(1, L ÷ 400)

    anim = @animate for t in [fill(1, fps); 1:tstep:L; fill(L, 3fps)]
        β = centroid[t][1:end-1]
        f(x) = evalpoly(x, β)
        P = plot(; legend=false)
        scatter!(x, y)
        plot!(xs, f.(xs); ylim=(-1.5, 1.5))
        title!("degree-$d fitting of size-$n data"; titlefontsize=10)

        Q = plot(; legend=:topright)
        plot!(-value[1:t]; label="(-1)×(log likelihood)")
        plot!(ngenerr[1:t]; label="n×(generalization error)")
        plot!(;xlim=extrema(axes(value, 1)).+(-0.05L, 0.05L), ylim)
        title!("t = $t"; titlefontsize=10)

        plot(P, Q; size=(800, 300))    end
    gif(anim, gifname; fps)
end

# %%
animate_polyfit(x, y, 1, f₀, σ₀)

# %%
animate_polyfit(x, y, 3, f₀, σ₀)

# %%
animate_polyfit(x, y, 9, f₀, σ₀)

# %%
animate_polyfit(x, y, 15, f₀, σ₀)

# %%
GE₀ = (1/2)*(log(2π*σ₀^2) + σ₀^2/σ₀^2)
n * GE₀

# %%
