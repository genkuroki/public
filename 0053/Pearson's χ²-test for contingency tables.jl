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
#     display_name: Julia 1.11.4
#     language: julia
#     name: julia-1.11
# ---

# %%
# using Pkg
# Pkg.add("Distributions")
# Pkg.add("StatsPlots")

# %%
using Random
using Distributions
using StatsPlots
default(fmt=:png)

_ecdf(A, x) = count(≤(x), A) / length(A)

safediv(x, y) = x == 0 ? zero(x/y) : x/y

function pearson_chi_squared(X::Matrix)
    N = sum(X)
    function Z(i, j)
        Xtilde_ij = sum(X[i,j] for j in axes(X, 2)) * sum(X[i,j] for i in axes(X, 1)) / N
        safediv((X[i,j] - Xtilde_ij)^2, Xtilde_ij)
    end
    sum(Z(i, j) for i in axes(X, 1), j in axes(X, 2))
end

"""
    rand_multinomials!(n, p, X)

例えば, `n = [10, 20]`, `p = [0.1, 0.2, 0.7]` であるとする.
そのとき, `X` が整数の2×3行列ならば,
`rand_multinomials!(n, p, X)` は行列`X`の第`i`行に
多項分布`Multinomial(n[i], p)`の乱数を格納する.
その返り値は変更後の`X`になる.

この関数で生成されるサンプルが従う確率分布はPearsonのχ²検定の帰無仮説を満たしている.

Example:
```
julia> X = zeros(Int, 2, 3);
julia> rand_multinomials!([10, 20], [0.1, 0.2, 0.7], X)
2×3 Matrix{Int64}:
 1  2   7
 2  4  14
```
"""
function rand_multinomials!(n, p, X)
    for i in eachindex(n)
        rand!(Multinomial(n[i], p), @view(X[i, :]))
    end
    X
end

function sim_chisq_stat(n, p; niters=10)
    chisq = zeros(niters)
    Xtmp = [zeros(Int, length(n), length(p)) for _ in 1:Threads.nthreads()]
    Threads.@threads for k in 1:niters
        tid = Threads.threadid()
        X = Xtmp[tid]
        rand_multinomials!(n, p, X)
        chisq[k] = pearson_chi_squared(X)
    end
    chisq
end

pvalue_chisq(chisq, df) = ccdf(Chisq(df), chisq)

"""
    plot_sim_chisq_stat(; n=[10, 15], p=[0.1, 0.2, 0.7], niters=10^5)

`rand_multinomials!(n, p, X)` によって繰り返しテストサンプルを生成して, Pearsonのχ²統計量の計算結果を記録し, Pearsonのχ²統計量の経験累積分布関数(ecdf)を計算してプロットする.  そのとき, その分布を近似するχ²分布の累積分布関数も同時プロットする.

その右にP値がα以下になる確率をプロットする.  適切なP値を使っていれば, 帰無仮説の下でのモデルの確率分布についてP値がα以下になる確率(実質アルファエラー率)はαで近似される.
"""
function plot_sim_chisq_stat(; n=[10, 15], p=[0.1, 0.2, 0.7], niters=10^6)
    df = (length(n) - 1) * (length(p) - 1)
    xmax = quantile(Chisq(df), 0.995)
    @show n
    @show p
    @show df
    @time chisq = sim_chisq_stat(n, p; niters)
    pval = pvalue_chisq.(chisq, df)
    
    P = plot(x -> _ecdf(chisq, x), 0, xmax; label="ecdf of Pearson's χ²-statistics")
    plot!(x -> cdf(Chisq(df), x), 0, xmax; label="cdf of χ²-distribution with df=$df", ls=:dash)
    plot!(legend=:bottomright, legendfontsize=12)
    plot!(xguide="Pearson's χ²")
    plot!(ytick=0:0.05:1)
    
    Q = plot(alpha -> _ecdf(pval, alpha), 0, 1; label="")
    plot!(identity; label="", ls=:dash)
    plot!(xguide="nominal significance level α", yguide="probability of P-value ≤ α")
    plot!(xtick=0:0.05:1, ytick=0:0.05:1, xrotation=90)
    
    plot(P, Q; size=(960, 400), layout=@layout[a b{0.4w}])
    plot!(bottommargin=10Plots.mm)
end

# %%
@doc rand_multinomials!

# %%
# このセルを繰り返し実行してみよ.

X = zeros(Int, 2, 3);
rand_multinomials!([10, 20], [0.1, 0.2, 0.7], X)

# %%
@doc plot_sim_chisq_stat

 # %%
 plot_sim_chisq_stat(; n=[6, 10], p=[0.2, 0.8])

 # %%
 plot_sim_chisq_stat(; n=[12, 20], p=[0.2, 0.8])

 # %%
 plot_sim_chisq_stat(; n=[24, 40], p=[0.2, 0.8])

 # %%
 plot_sim_chisq_stat(; n=[48, 80], p=[0.2, 0.8])

 # %%
 plot_sim_chisq_stat(; n=[100, 120], p=[0.05, 0.95])

 # %%
 plot_sim_chisq_stat(; n=[10, 20], p=[0.1, 0.2, 0.7])

 # %%
 plot_sim_chisq_stat(; n=[20, 40], p=[0.1, 0.2, 0.7])

# %%
plot_sim_chisq_stat(; n=[8, 10, 12], p=[0.1, 0.2, 0.3, 0.4])

# %%
plot_sim_chisq_stat(; n=[16, 20, 24], p=[0.1, 0.2, 0.3, 0.4])

# %%
