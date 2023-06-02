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
#     display_name: Julia 1.9.0
#     language: julia
#     name: julia-1.9
# ---

# %%
using BenchmarkTools
using Combinatorics
using Distributions
using HypothesisTests
using Memoization
using StatsBase
using StatsPlots
default(fmt=:png)

# %% [markdown]
# $\{1,2,\ldots,N\}$ から重複無しに $m$ 個の組み合わせ $1\le i_1<\cdots<i_m\le N$ を取るとき, その和が $r$ になる場合の数を $f(N, m, r)$ と書くと, 以下が成立することがわかる:
#
# $$
# \begin{aligned}
# &
# f(N, m, r) = 0\quad\text{unless $0\le m \le N$ and $0\le r - m(m+1)/2 \le m(N-m)$},
# \\ &
# f(N, 0, 0) = 1,
# \\ &
# f(N, 1, r) = \text{if $1\le r\le N$ then $1$ else $0$},
# \\ &
# f(N, m, r) = f(N-1, m, r) + f(N-1, m-1, r-N).
# \end{aligned}
# $$
#
# ただし, この方法で計算するにはメモ化がほぼ必須である.
#
# 以下のセルのコードは上の方法による $f(N, m, r)$ の計算法をそのままコードに翻訳したものである.

# %%
function f_naive(N, m, r)
    0 ≤ m ≤ N || return zero(N)
    0 ≤ r - m*(m+1)÷2 ≤ m*(N-m) || return zero(N)
    m == 0 && return one(N)
    m == 1 && return 1 ≤ r ≤ N ? one(N) : zero(N)
    f_naive(N-1, m, r) + f_naive(N-1, m-1, r-N)
end

@memoize function f_memoize(N, m, r)
    0 ≤ m ≤ N || return zero(N)
    0 ≤ r - m*(m+1)÷2 ≤ m*(N-m) || return zero(N)
    m == 0 && return one(N)
    m == 1 && return 1 ≤ r ≤ N ? one(N) : zero(N)
    f_memoize(N-1, m, r) + f_memoize(N-1, m-1, r-N)
end

# %% [markdown]
# メモ化していない `f_naive` はひどく遅い.

# %%
for m in 10:18
    N = 2m
    r = m*(m+1)÷2 + m*(N-m)÷2
    print("m = $m, N = $N, r = $r: ")
    @time f_naive(N, m, r)
end

# %% [markdown]
# メモ化された `f` はこれよりかなり速くなっている.

# %%
for m in 10:18
    N = 2m
    r = m*(m+1)÷2 + m*(N-m)÷2
    print("m = $m, N = $N, r = $r: ")
    @time f_memoize(N, m, r)
end

# %% [markdown]
# 念のために正しく計算されていることを確認しておこう.
#
# 以下では $N = n + m$ とおく.

# %%
function F_memoize(n, m)
    rmin = m*(m+1)÷2
    rmax = rmin + m*n
    rs = rmin:rmax
    fs = [f_memoize(n+m, m, r) for r in rs]
    rs, fs
end

function F_comb(n, m)
    rank_sums = [sum(comb) for comb in combinations(1:(n+m), m)]
    c = countmap(rank_sums)
    rmin, rmax = extrema(keys(c))
    rs = rmin:rmax
    fs = [c[k] for k in rs]
    rs, fs
end

# %%
N = 10
F1 = F_memoize
F2 = F_comb
for m in 0:N
    @show N, m, F1(N-m, m)
    @show N, m, F2(N-m, m)
    @show F1(N-m, m) == F2(N-m, m)
    println()
end

# %%
[F_memoize(n, m) == F_comb(n, m) for n in 0:10, m in 0:10]

# %%
for m in 25:5:50
    n = m
    print("n = $n, m = $m: ")
    @time F_memoize(m, m)
end

# %% [markdown]
# $$
# \begin{aligned}
# &
# f(N, m, r) = 0\quad\text{unless $0\le m \le N$ and $0\le r - m(m+1)/2 \le m(N-m)$},
# \\ &
# f(N, 0, 0) = 1,
# \\ &
# f(N, 1, r) = \text{if $1\le r\le N$ then $1$ else $0$},
# \\ &
# f(N, m, r) = f(N-1, m, r) + f(N-1, m-1, r-N).
# \end{aligned}
# $$
#
# より,
#
# $$
# N = m + n, \quad
# r = m(m+1)/2 + u, \quad
# g(n, m, u) = f(m+n, m, m(m+1)/2 + u)
# $$
#
# とおくと,
#
# $$
# \begin{aligned}
# &
# g(n, m, u) = 0\quad\text{unless $0 \le m$, $0\le n$, and $0\le u \le mn$},
# \\ &
# g(n, 0, 0) = 1,
# \\ &
# g(n, 1, u) = \text{if $0\le u \le m+n-1$ then $1$ else $0$},
# \\ &
# g(n, m, u) = g(n-1, m, u) + g(n, m-1, u-n).
# \end{aligned}
# $$

# %% [markdown]
# 特に,
#
# $$
# g(0, 1, u) = \text{if $u=0$ then $1$ else $0$}.
# $$
#
# $g(n, m, u)$ の値は以下を帰納的に計算すれば得られる:
#
# $$
# g(j, i, u-(n-j)n), \quad 0\le j\le n, \quad 0\le i \le m.
# $$

# %%
function G(n, m,
        prev = Matrix{Int}(undef, m*n+1, m+1),
        next = similar(prev),
    )
    # g(0, i, u) = δ_{i0} δ_{u0}
    @. next = 0
    next[1+0, 1+0] = 1
    for k in 1:m+n
        @. prev = next
        for i in max(0, k-n):min(m, k)
            j = k - i
            # g(j, i, u) = g(j-1, i, u) + g(j, i-1, u-j)
            @views @. next[1:1+i*(j-1), 1+i] = prev[1:1+i*(j-1), 1+i]
            if i ≥ 1
                @views @. next[1+j:1+i*j, 1+i] += prev[1:1+(i-1)*j, 1+i-1]
            end
            # g(0, i, u) = δ_{i0} δ_{u0}
            if j == 0
                next[1+0, 1+i] = 1
            end
        end
    end
    @view next[:, 1+m]
end

function F_recursive(n, m,
        prev = Matrix{Int}(undef, m*n+1, m+1),
        next = similar(prev),
    )
    rmin = m*(m+1)÷2
    rmax = rmin + m*n
    rs = rmin:rmax
    fs = G(n, m, prev, next)
    rs, fs
end

# %%
N = 10
F1 = F_memoize
F2 = F_recursive
for m in 0:N
    @show N, m, F1(N-m, m)
    @show N, m, F2(N-m, m)
    @show F1(N-m, m) == F2(N-m, m)
    println()
end

# %%
[F_recursive(n, m) == F_memoize(n, m) for n in 0:20, m in 0:20]

# %%
for m in 25:5:100
    n = m
    print("n = $n, m = $m: ")
    prev = Matrix{Int}(undef, m*n+1, m+1)
    next = similar(prev)
    @time G(n, m, prev, next)
end

# %%
n, m = 100, 100
prev = Matrix{Int}(undef, m*n+1, m+1)
next = similar(prev)
@btime F_recursive($n, $m, $prev, $next);

# %%
n, m = 100, 100
x, y = rand(m), rand(n)
@time ExactMannWhitneyUTest(x, y)
@time ExactMannWhitneyUTest(x, y)
@time ExactMannWhitneyUTest(x, y)

# %%
n, m = 100, 100
prev = Matrix{Int}(undef, m*n+1, m+1)
next = similar(prev)
rs, fs = @time F_recursive(n, m, prev, next)
@show rs length(fs) typeof(fs)

ps = fs/sum(fs)
@show mu = m*(m+n+1)/2
@show s2 = m*n*(m+n+1)/12

plot(rs[begin:10:end], ps[begin:10:end]; label="")
plot!(Normal(mu, √s2), extrema(rs)...; ls=:dash, label="")

# %%
extrema(fs)

# %%
n, m = 100, 100
prev = Matrix{Int128}(undef, m*n+1, m+1)
next = similar(prev)
rs, fs = @time F_recursive(n, m, prev, next)
@show rs length(fs) typeof(fs)

ps = fs/sum(fs)
@show mu = m*(m+n+1)/2
@show s2 = m*n*(m+n+1)/12

plot(rs[begin:10:end], ps[begin:10:end]; label="")
plot!(Normal(mu, √s2), extrema(rs)...; ls=:dash, label="")

# %%
extrema(fs)

# %%
n, m = 100, 100
prev = Matrix{BigInt}(undef, m*n+1, m+1)
next = similar(prev)
rs, fs = @time F_recursive(n, m, prev, next)
@show rs length(fs) typeof(fs)

ps = fs/sum(fs)
@show mu = m*(m+n+1)/2
@show s2 = m*n*(m+n+1)/12

plot(rs[begin:10:end], ps[begin:10:end]; label="")
plot!(Normal(mu, √s2), extrema(rs)...; ls=:dash, label="")

# %%
extrema(fs)

# %%
n, m = 100, 100
prev = Matrix{Float64}(undef, m*n+1, m+1)
next = similar(prev)
rs, fs = @time F_recursive(n, m, prev, next)
@show rs length(fs) typeof(fs)

ps = fs/sum(fs)
@show mu = m*(m+n+1)/2
@show s2 = m*n*(m+n+1)/12

plot(rs[begin:10:end], ps[begin:10:end]; label="")
plot!(Normal(mu, √s2), extrema(rs)...; ls=:dash, label="")

# %%
extrema(fs)

# %%
sum(fs)

# %%
binomial(big(200), 100)

# %%
binomial(big(200), 100) |> float

# %%
