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
#     display_name: Julia 1.10.0
#     language: julia
#     name: julia-1.10
# ---

# %% [markdown]
# # 統計力学におけるカノニカル分布の最も簡単な場合
#
# 黒木玄
#
# 2020-09-24～2020-09-26, 2024-01-20
#
# Reference: [Kullback-Leibler 情報量と Sanov の定理](https://genkuroki.github.io/documents/20160616KullbackLeibler.pdf)

# %% [markdown] toc=true
# <h1>目次<span class="tocSkip"></span></h1>
# <div class="toc"><ul class="toc-item"><li><span><a href="#カノニカル分布が指数分布になる場合" data-toc-modified-id="カノニカル分布が指数分布になる場合-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>カノニカル分布が指数分布になる場合</a></span><ul class="toc-item"><li><span><a href="#成分の分布が指数分布で近似されるランダムベクトルを直接的に生成" data-toc-modified-id="成分の分布が指数分布で近似されるランダムベクトルを直接的に生成-1.1"><span class="toc-item-num">1.1&nbsp;&nbsp;</span>成分の分布が指数分布で近似されるランダムベクトルを直接的に生成</a></span></li><li><span><a href="#近似的に指数分布が得られるMCMC法(1)" data-toc-modified-id="近似的に指数分布が得られるMCMC法(1)-1.2"><span class="toc-item-num">1.2&nbsp;&nbsp;</span>近似的に指数分布が得られるMCMC法(1)</a></span></li><li><span><a href="#近似的に指数分布が得られるMCMC法(2)" data-toc-modified-id="近似的に指数分布が得られるMCMC法(2)-1.3"><span class="toc-item-num">1.3&nbsp;&nbsp;</span>近似的に指数分布が得られるMCMC法(2)</a></span></li></ul></li><li><span><a href="#カノニカル分布が正規分布になる場合-(Maxwell–Boltzmann分布)" data-toc-modified-id="カノニカル分布が正規分布になる場合-(Maxwell–Boltzmann分布)-2"><span class="toc-item-num">2&nbsp;&nbsp;</span>カノニカル分布が正規分布になる場合 (Maxwell–Boltzmann分布)</a></span><ul class="toc-item"><li><span><a href="#成分の分布が正規分布で近似されるランダムベクトルを直接的に生成" data-toc-modified-id="成分の分布が正規分布で近似されるランダムベクトルを直接的に生成-2.1"><span class="toc-item-num">2.1&nbsp;&nbsp;</span>成分の分布が正規分布で近似されるランダムベクトルを直接的に生成</a></span></li><li><span><a href="#近似的に正規分布が得られるMCMC法(2)" data-toc-modified-id="近似的に正規分布が得られるMCMC法(2)-2.2"><span class="toc-item-num">2.2&nbsp;&nbsp;</span>近似的に正規分布が得られるMCMC法(2)</a></span></li></ul></li><li><span><a href="#カノニカル分布がガンマ分布になる場合" data-toc-modified-id="カノニカル分布がガンマ分布になる場合-3"><span class="toc-item-num">3&nbsp;&nbsp;</span>カノニカル分布がガンマ分布になる場合</a></span><ul class="toc-item"><li><span><a href="#近似的にガンマ分布が得られるMCMC法(2)" data-toc-modified-id="近似的にガンマ分布が得られるMCMC法(2)-3.1"><span class="toc-item-num">3.1&nbsp;&nbsp;</span>近似的にガンマ分布が得られるMCMC法(2)</a></span></li></ul></li><li><span><a href="#特に名前がないと思われる分布に収束する場合" data-toc-modified-id="特に名前がないと思われる分布に収束する場合-4"><span class="toc-item-num">4&nbsp;&nbsp;</span>特に名前がないと思われる分布に収束する場合</a></span><ul class="toc-item"><li><span><a href="#NLsolve.jl-の使い方の確認" data-toc-modified-id="NLsolve.jl-の使い方の確認-4.1"><span class="toc-item-num">4.1&nbsp;&nbsp;</span>NLsolve.jl の使い方の確認</a></span></li><li><span><a href="#特に名前がないと思われる分布への収束の確認" data-toc-modified-id="特に名前がないと思われる分布への収束の確認-4.2"><span class="toc-item-num">4.2&nbsp;&nbsp;</span>特に名前がないと思われる分布への収束の確認</a></span></li><li><span><a href="#特に名前がないと思われる分布への収束のアニメーション" data-toc-modified-id="特に名前がないと思われる分布への収束のアニメーション-4.3"><span class="toc-item-num">4.3&nbsp;&nbsp;</span>特に名前がないと思われる分布への収束のアニメーション</a></span></li></ul></li></ul></div>

# %% [markdown]
# ## カノニカル分布が指数分布になる場合

# %% [markdown]
# ### 成分の分布が指数分布で近似されるランダムベクトルを直接的に生成

# %%
# 総和がNのN個の正の実数を成分とするランダムベクトルの生成

using Random: AbstractRNG, SamplerTrivial, rand!, default_rng

struct Simplex{T} N::Int end
Simplex(N) = Simplex{Float64}(N)

function Base.rand(rng::AbstractRNG, d::SamplerTrivial{Simplex{T}}) where T
    N = d[].N
    X = Vector{T}(undef, N)
    rand!(rng, @view(X[begin:end-1]))
    X[end] = 1
    X .*= N
    sort!(X)
    @. @views X[begin+1:end] -= X[begin:end-1]
    X
end

# %%
rand(Simplex(3))

# %%
rand(Simplex(10)) |> sum

# %%
rand(Simplex(3), 5)

# %%
# 総和がNのN個の正の実数を成分とするランダムベクトルの成分の分布は
# Nを大きくすると期待値1の指数分布に近付く.
# これは「統計力学」の最も易しい場合になっている.

using Plots; default(fmt=:png)

N = 10^4
histogram(rand(Simplex(N)); norm=true, alpha=0.3, label="rand(Simplex($N))")
plot!(x->exp(-x), 0, 8; label="exp(-x)", lw=2)
plot!(xlim=(0, 8))

# %% [markdown]
# ### 近似的に指数分布が得られるMCMC法(1)
#
# N = length(X)が大きなとき, Xのすべての成分を正に保ったままで, N == sum(X)が不変になるようなランダムウォークを繰り返すと, Xの成分の分布は期待値1の指数分布で近似されるようになる.

# %%
# 富のランダム再分配

function randomly_redistribute!(rng::AbstractRNG, X, h)
    N = length(X)
    i = rand(rng, 1:N)
    j = rand(rng, 1:N-1)
    j = ifelse(j == i, N, j)
    dx = h*abs(randn())
    if dx < X[i]
        X[i] -= dx
        X[j] += dx
    end
    X
end

function iter_redistribute!(X, L, h; rng = default_rng())
    for _ in 1:L
        randomly_redistribute!(rng, X, h)
    end
    X
end

# %%
N = 10^4
X = ones(N)
L = 10^6
h = 0.02
tmax = 200

anim = @animate for t in [fill(0, 10); 1:tmax-1; fill(tmax, 20)]
    t == 0 || t == tmax || iter_redistribute!(X, L, h)
    histogram(X; norm=true, alpha=0.3, label="MCMC")
    plot!(x -> exp(-x), 0, 8; label="exp(-x)", lw=2)
    plot!(xlim=(-0.1, 8), ylim=(-0.02, 1.05))
    title!("t = $t"; titlefontsize=12)
end

gif(anim, "exp1.gif", fps=20)

# %% [markdown]
# ### 近似的に指数分布が得られるMCMC法(2)
#
# N = length(X)が大きなとき, Xのすべての成分が正でかつ, sum(X)がN以下であるという条件を保ったままランダムウォークを繰り返すと, Xの成分の分布は期待値1の指数分布で近似されるようになる.  sum(X)を**ぴったりNに固定**せずに, **N以下**に広げても, 収束する先の分布は変わらない.
#
# 詳しい一般論については, [『Kullback-Leibler情報量とSanovの定理』の第4.2節](https://genkuroki.github.io/documents/20160616KullbackLeibler.pdf#page23)を参照せよ. そこではエネルギー $E_i$ の平均値をぴったり固定せずに不等号で制限している.

# %%
function iter_randwalk!(X, L, h; rng = default_rng())
    N = length(X)
    S = sum(X)
    if S > N
        X ./= 1.1S
    end
    for _ in 1:L
        i = rand(rng, 1:N)
        dx = h*randn()
        if X[i] + dx > 0 && S + dx ≤ N
            X[i] += dx
            S += dx
        end
    end
    X
end

# %%
N = 10^4
X = 0.8ones(N) # テキトーな初期条件
L = 10^6
h = 0.02
tmax = 200

anim = @animate for t in [fill(0, 10); 1:tmax-1; fill(tmax, 20)]
    t == 0 || t == tmax || iter_randwalk!(X, L, h)
    histogram(X; norm=true, alpha=0.3, label="MCMC")
    plot!(x -> exp(-x), 0, 8; label="exp(-x)", lw=2)
    plot!(xlim=(-0.1, 8), ylim=(-0.02, 1.05))
    title!("t = $t"; titlefontsize=12)
end

gif(anim, "exp2.gif", fps=20)

# %% [markdown]
# ## カノニカル分布が正規分布になる場合 (Maxwell–Boltzmann分布)
#
# * [『Kullback-Leibler 情報量と Sanov の定理』第2.4節](https://genkuroki.github.io/documents/20160616KullbackLeibler.pdf#page23)を参照せよ.

# %% [markdown]
# ### 成分の分布が正規分布で近似されるランダムベクトルを直接的に生成

# %%
# Euclidノルムの二乗がNのランダムベクトルを生成

using Random: AbstractRNG, SamplerTrivial, rand!, default_rng
using LinearAlgebra: norm, dot
norm2(x) = dot(x, x)

struct Sphere{T} N::Int end
Sphere(N) = Sphere{Float64}(N)

function Base.rand(rng::AbstractRNG, d::SamplerTrivial{Sphere{T}}) where T
    N = d[].N
    X = randn(rng, T, N)
    X .*= √N/norm(X)
    X
end

# %%
X = rand(Sphere(4))

# %%
X = rand(Sphere(10)) |> norm2

# %%
rand(Sphere(3), 5)

# %%
# Euclidノルムの二乗がNのランダムベクトルの成分の分布は
# Nを大きくすると標準正規分布に近付く.
# これは「統計力学」におけるMaxwell–Boltzmann分布の場合になっている.

using Plots
N = 10^4
histogram(rand(Sphere(N)); norm=true, alpha=0.3, label="rand(Sphere($N))")
plot!(x->exp(-x^2/2)/√(2π), -4, 4; label="standard normal dist.", lw=2)

# %% [markdown]
# ### 近似的に正規分布が得られるMCMC法(2)
#
# N = length(X)が大きなとき, XのEuclidノルムの2乗がN以下であるという条件を保ったままランダムウォークを繰り返すと, Xの成分の分布は標準正規分布で近似されるようになる. XのEuclidノルムの2乗を**ぴったりNに固定**せずに, **N以下**に広げても, 収束する先の分布は変わらない.
#
# 詳しい一般論については, [『Kullback-Leibler情報量とSanovの定理』の第4.2節](https://genkuroki.github.io/documents/20160616KullbackLeibler.pdf#page31)を参照せよ. そこではエネルギー $E_i$ の平均値をぴったり固定せずに不等号で制限している.

# %%
function iter_randwalk_in_ball!(X, L, h; rng = default_rng())
    N = length(X)
    E = norm2(X) # E stands for the total Energy
    if E > N
        X ./= 1.1*√E
    end
    for _ in 1:L
        i = rand(rng, 1:N)
        dx = h*randn()
        dE = 2*X[i]*dx + dx^2
        if E + dE ≤ N
            X[i] += dx
            E += dE
        end
    end
    X
end

# %%
N = 10^4
X = [0.8 .+ 0.5rand(N ÷ 2); -0.9 .+ 0.2randn(N ÷ 2)] # テキトーな初期条件
L = 10^6
h = 0.01
tmax = 200

anim = @animate for t in [fill(0, 10); 1:tmax-1; fill(tmax, 20)]
    t == 0 || t == tmax || iter_randwalk_in_ball!(X, L, h)
    histogram(X; norm=true, alpha=0.3, label="MCMC")
    plot!(x->exp(-x^2/2)/√(2π), -4, 4; label="standard normal dist.", lw=2)
    plot!(xlim=(-4, 4), ylim=(-0.02, 0.5))
    title!("t = $t"; titlefontsize=12)
end

gif(anim, "normal2.gif", fps=20)

# %% [markdown]
# ## カノニカル分布がガンマ分布になる場合
#
# *  [『Kullback-Leibler情報量とSanovの定理』の第4.2節](https://genkuroki.github.io/documents/20160616KullbackLeibler.pdf#page31)を参照せよ.

# %% [markdown]
# ### 近似的にガンマ分布が得られるMCMC法(2)
#
# N = length(X)が大きなとき, Xの各成分が正で加法平均がμ以下であるという条件だけではなく, 対数平均がある値ν以上になるという条件を保ったままランダムウォークを繰り返すと, Xの成分の分布はガンマ分布で近似されるようになる.

# %%
using Random: default_rng
using Plots; default(fmt=:png)
using SpecialFunctions
using Distributions
using Printf
rd(x, d=2) = round(x; digits=d)
using ProgressMeter

function iter_randwalk_gamma!(X, μ, ν, L, h; rng = default_rng())
    N = length(X)
    S = mean(X)
    U = mean(log, X)
    for _ in 1:L
        i = rand(rng, 1:N)
        dx = h*randn()
        X_i_new = X[i] + dx
        if X_i_new > 0
            S_new = S + dx/N
            U_new = U + (-log(X[i]) + log(X_i_new))/N
            if S_new ≤ μ && U_new ≥ ν
                X[i] = X_i_new
                S = S_new
                U = U_new
            end
        end
    end
    X
end

Delta(α) = log(α) - digamma(α)

# %%
plot(Delta, 0.2, 10; label="Delta(α)", size=(400, 250))

# %%
function anim_gamma!(μ, α, X, L, h, tmax, xlim=(0, 10), ylim=(0, 0.6))
    θ = μ/α
    Δ = Delta(α)
    ν = log(μ) - Δ
    
    S = [mean(X)]
    U = [mean(log, X)]
    
    prog = Progress(tmax+30, 0)
    anim = @animate for t in [fill(0, 10); 1:tmax; fill(tmax+1, 20)]
        if 1 ≤ t ≤ tmax
            iter_randwalk_gamma!(X, μ, ν, L, h)
            push!(S, mean(X))
            push!(U, mean(log, X))
        end

        P = histogram(X; norm=true, alpha=0.3, label="MCMC")
        title!("μ = $(rd(μ)),  ν = log(μ) - $(rd(Δ,4)),  t = $t")
        f_label = @sprintf "Gamma(%.2f, %.2f)" α θ
        plot!(x -> pdf(Gamma(α, θ), x), xlim...; label=f_label, color=:red, lw=2)
        plot!(; xlim, ylim)

        Q = hline([μ]; label="", color=:red, lw=2)
        plot!(0:min(tmax,t), S; label="", color=:blue)
        title!("mean(X)")
        plot!(xlim=(0, tmax), ylim=(minimum(S), μ + 0.05*(μ - minimum(S))))

        R = hline([ν]; label="", color=:red, lw=2)
        plot!(0:min(tmax,t), U; label="", color=:blue)
        title!("mean(log, X)")
        plot!(xlim=(0, tmax), ylim=(ν - 0.05*(maximum(U) - ν), maximum(U)))

        layout = @layout[a [b; c]]
        plot(P, Q, R; layout, size=(800, 300))
        plot!(legendfontsize=8, titlefontsize=10)
        
        next!(prog)
    end

    gif(anim, @sprintf("Gamma(%.2f, %.2f).gif", α, θ), fps=20)
end

# %%
μ = 2.0
α = 1.0

N = 10^4
X = 0.7μ*ones(N)
L = 10^6
h = 0.05
tmax = 200

xlim=(0, 10)
ylim=(0, 0.6)

anim_gamma!(μ, α, X, L, h, tmax, xlim, ylim)

# %%
μ = 2.0
α = 5.0

N = 10^4
X = 0.95μ*ones(N)
L = 8*10^5
h = 0.01
tmax = 200

xlim = (0, 8)
ylim = (0, 0.6)

anim_gamma!(μ, α, X, L, h, tmax, xlim, ylim)

# %%
μ = 2.0
α = 25.0

N = 10^4
X = 0.99μ*ones(N)
L = 5*10^5
h = 0.005
tmax = 200

xlim = (0, 5)
ylim = (0, 1.2)

anim_gamma!(μ, α, X, L, h, tmax, xlim, ylim)

# %% [markdown]
# ## 特に名前がないと思われる分布に収束する場合

# %% [markdown]
# ### NLsolve.jl の使い方の確認

# %%
using QuadGK
using NLsolve

bf(x, a, f=√) = exp(-a[1]*x + (a[2]-1)*f(x))
Zf(a, f=√, xmax=1e3) = quadgk(x -> bf(x, a, f), 0, xmax)[1]

function Eq!(p, a, f=√, μ=2.0, ν=√μ - 0.05, xmax=1e3)
    Z = Zf(a, f, xmax)
    p[1] = quadgk(x ->    x*bf(x, a, f), 0, xmax)[1]/Z - μ
    p[2] = quadgk(x -> f(x)*bf(x, a, f), 0, xmax)[1]/Z - ν
    p
end

function solve_Eq(f=√, μ=2.0, ν=√μ - 0.05, xmax=1e3)
    F!(p, a) = Eq!(p, a, f, μ, ν, xmax)
    sol = nlsolve(F!, ones(2))
    a = sol.zero
    Z = Zf(a, f, xmax)
    p(x) = bf(x, a, f)/Z
    (a = a, p = p)
end

# %%
@time sol = nlsolve(Eq!, ones(2))

# %% [markdown]
# ↑コンパイルによる遅延時間が含まれている.
#
# 数値積分を3回実行して計算される函数の零点を見付ける問題を解かせることになるので, 実行時間が長くなり過ぎることを心配したが, Iterations: 8 なので大丈夫そうである.

# %%
dump(sol)

# %%
@time solve_Eq()

# %% [markdown]
# ↑二度目の実行はそこそこ速い.

# %% [markdown]
# ### 特に名前がないと思われる分布への収束の確認
#
# 以下, $N$ は十分に大きいと仮定する.
#
# $f(x)$ は正の実数 $x$ について定義された, 上に凸な狭義単調増加函数であるとする(所謂効用函数).  パラメーター $a=(a_1,a_2)$ を持つ $x$ の確率密度函数 $p(x|a)$ を
#
# $$
# p(x|a) = Z(a,b)^{-1}e^{-a_1 x + (a_2 - 1)f(x)}, \quad
# Z(a,b) = \int_0^\infty e^{-a_1 x + (a_2 - 1)f(x)}\,dx
# $$
#
# と定める.  $\mu$ は正の実数であるとし, $\nu$ は
#
# $$
# \frac{1}{\mu}\int_0^\infty f(x)e^{-x/\mu}\,dx =
# \int_0^\infty f(x)p(x|1/\mu, 1)\,dx \le \nu < \log\mu
# $$
#
# を満たしていると仮定する. $a=(a_1,a_2)$ をこのような $\mu, \nu$ から
#
# $$
# \int_0^\infty  x\,p(x|a)\,dx = \mu, \quad
# \int_0^\infty f(x)p(x|a)\,dx = \nu
# $$
#
# という条件によって定める.
#
# 実数達 $X_i$ ($i=1,2,\ldots,N$) の範囲に以下の制限を課す:
#
# $$
# X_i > 0, \quad
# \frac{1}{N}\sum_{i=1}^N X_i \le \mu, \quad
# \frac{1}{N}\sum_{i=1}^N f(X_i) \ge \nu.
# $$
#
# この条件のもとで, $X=(X_1,\ldots,X_N)$ をランダムウォークさせると, $X_i$ 達のなす分布(大雑把には $X_i$ 達のヒストグラムだとみなしてよい)は, 上のように $\mu, \nu$ から $a=(a_1,a_2)$ を定めたときの確率密度函数 $p(x|a)$ を持つ分布で近似される. 
#
# このことを $f(x) = \sqrt{x}$ の場合について数値的に確認しよう.

# %%
using Random: default_rng
using Plots; default(fmt=:png)
using Printf
rd(x, d=2) = round(x; digits=d)

function iter_randwalk_lbfun!(f, X, μ, ν, L, h; rng = default_rng())
    N = length(X)
    S = mean(X)
    U = mean(f, X)
    for _ in 1:L
        i = rand(rng, 1:N)
        dx = h*randn()
        X_i_new = X[i] + dx
        if X_i_new > 0
            S_new = S + dx/N
            U_new = U + (-f(X[i]) + f(X_i_new))/N
            if S_new ≤ μ && U_new ≥ ν
                X[i] = X_i_new
                S = S_new
                U = U_new
            end
        end
    end
    X
end

# %%
μ = 2.0
ν = √μ - 1e8
νmin = quadgk(x -> √x*bf(x, [1/μ, 1], √), 0, 1e3)[1]/Zf([1/μ, 1], √)
@show νmin
@show Δmax = √μ - νmin 

N = 10^4
X = 0.99μ*ones(N)
L = 10^6
h = 0.05
tmax = 100

@time iter_randwalk_lbfun!(√, X, μ, ν, tmax*L, h)
@show mean(√, X)
@show √μ - mean(√, X)

histogram(X; norm=true, alpha=0.3, label="")
plot!(x -> exp(-x/μ)/μ, 0, 15; label="", lw=2)

# %%
μ = 2.0
ν = √μ - Δmax

N = 10^4
X = 0.99μ*ones(N)
L = 10^6
h = 0.05
tmax = 100

@time iter_randwalk_lbfun!(√, X, μ, ν, tmax*L, h)
@show mean(√, X)
@show Δmax = √μ - mean(√, X)

histogram(X; norm=true, alpha=0.3, label="")
plot!(x -> exp(-x/μ)/μ, 0, 15; label="", lw=2)

# %%
μ = 2.0
ν = √μ - 0.1

@time a, p = solve_Eq(√, μ, ν)

N = 10^4
X = 0.99μ*ones(N)
L = 10^6
h = 0.05
tmax = 100

@time iter_randwalk_lbfun!(√, X, μ, ν, tmax*L, h)

histogram(X; norm=true, alpha=0.3, label="")
plot!(p, 0, 12; label="", lw=2)

# %%
μ = 2.0
ν = √μ - 0.05

@time a, p = solve_Eq(√, μ, ν)

N = 10^4
X = 0.99μ*ones(N)
L = 10^6
h = 0.05
tmax = 100

@time iter_randwalk_lbfun!(√, X, μ, ν, tmax*L, h)

histogram(X; norm=true, alpha=0.3, label="")
plot!(p, 0, 8; label="", lw=2)

# %%
μ = 2.0
ν = √μ - 0.02

@time a, p = solve_Eq(√, μ, ν)

N = 10^4
X = 0.99μ*ones(N)
L = 10^6
h = 0.05
tmax = 100

@time iter_randwalk_lbfun!(√, X, μ, ν, tmax*L, h)

histogram(X; norm=true, alpha=0.3, label="")
plot!(p, 0, 5; label="", lw=2)

# %%
μ = 2.0
ν = √μ - 0.01

@time a, p = solve_Eq(√, μ, ν)

N = 10^4
X = 0.99μ*ones(N)
L = 10^6
h = 0.05
tmax = 100

@time iter_randwalk_lbfun!(√, X, μ, ν, tmax*L, h)

histogram(X; norm=true, alpha=0.3, label="")
plot!(p, 0, 5; label="", lw=2)

# %% [markdown]
# 以上を見ると理論通りの分布に収束していることが分かる.

# %% [markdown]
# ###  特に名前がないと思われる分布への収束のアニメーション

# %%
using ProgressMeter

function anim_lbfun!(f, μ, Δ, X, L, h, tmax, xlim=(0, 10), ylim=(0, 0.6))
    ν = f(μ) - Δ
    a, p = solve_Eq(f, μ, ν)
    
    M = [mean(X)]
    U = [mean(f, X)]
    
    prog = Progress(tmax+30, 0)
    anim = @animate for t in [fill(0, 10); 1:tmax; fill(tmax+1, 20)]
        if 1 ≤ t ≤ tmax
            iter_randwalk_lbfun!(f, X, μ, ν, L, h)
            push!(M, mean(X))
            push!(U, mean(f, X))
        end

        P = histogram(X; norm=true, alpha=0.3, label="MCMC")
        title!("μ = $(rd(μ)),  ν = $f(μ) - $(rd(Δ,4)),  t = $t")
        f_label = @sprintf "Dist(%s, %.2f, %.2f)" f a[1] a[2]
        plot!(p, xlim...; label=f_label, color=:red, lw=2)
        plot!(; xlim, ylim)

        Q = hline([μ]; label="", color=:red, lw=2)
        plot!(0:min(tmax,t), M; label="", color=:blue)
        title!("mean(X)")
        plot!(xlim=(0, tmax), ylim=(minimum(M), μ + 0.05*(μ - minimum(M))))

        R = hline([ν]; label="", color=:red, lw=2)
        plot!(0:min(tmax,t), U; label="", color=:blue)
        title!("mean($f, X)")
        plot!(xlim=(0, tmax), ylim=(ν - 0.05*(maximum(U) - ν), maximum(U)))

        layout = @layout[a [b; c]]
        plot(P, Q, R; layout, size=(800, 300))
        plot!(legendfontsize=8, titlefontsize=10)
        
        next!(prog)
    end

    gif(anim, @sprintf("Dist(%s, %.2f, %.2f).gif", f, a[1], a[2]), fps=20)
end

# %%
μ = 2.0
Δ = 0.16

N = 10^4
X = 0.99μ*ones(N)
L = 5*10^5
h = 0.05
tmax = 200

xlim = (0, 15)
ylim = (0, 0.6)

anim_lbfun!(√, μ, Δ, X, L, h, tmax, xlim, ylim)

# %%
μ = 2.0
Δ = 0.1

N = 10^4
X = 0.99μ*ones(N)
L = 10^5
h = 0.05
tmax = 200

xlim = (0, 12)
ylim = (0, 0.5)

anim_lbfun!(√, μ, Δ, X, L, h, tmax, xlim, ylim)

# %%
μ = 2.0
Δ = 0.05

N = 10^4
X = 0.99μ*ones(N)
L = 10^5
h = 0.05
tmax = 200

xlim = (0, 8)
ylim = (0, 0.5)

anim_lbfun!(√, μ, Δ, X, L, h, tmax, xlim, ylim)

# %%
μ = 2.0
Δ = 0.02

N = 10^4
X = 0.99μ*ones(N)
L = 2*10^4
h = 0.05
tmax = 200

xlim = (0, 5)
ylim = (0, 0.7)

anim_lbfun!(√, μ, Δ, X, L, h, tmax, xlim, ylim)

# %%
μ = 2.0
Δ = 0.01

N = 10^4
X = 0.99μ*ones(N)
L = 10^4
h = 0.05
tmax = 200

xlim = (0, 5)
ylim = (0, 0.9)

anim_lbfun!(√, μ, Δ, X, L, h, tmax, xlim, ylim)

# %%
