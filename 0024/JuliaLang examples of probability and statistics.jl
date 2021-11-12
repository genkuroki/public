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
#     display_name: Julia 1.6.3
#     language: julia
#     name: julia-1.6
# ---

# %%
using Distributions
using StatsPlots

# warm up (少し時間を取られる。その間に下のセルでコードを入力する)
histogram(randn(10^5); norm=true, alpha=0.3, size=(300, 200))
plot!(Normal())

# %%
# コインを20回投げたときの表が出る回数の分布は正規分布で近似される。

n = 20
p = 0.5
sample_of_means = rand(Binomial(n, p), 10^5)
histogram(sample_of_means; norm=true, alpha=0.3, bin=-0.5:n+0.5, label="rand(Binomial($n, $p))")
plot!(Normal(n*p, √(n*p*(1-p))), 0, 20; label="normal approx.", lw=2)

# %%
# 回数n、確率パラメータpの二項分布オブジェクトを次のように作れる。

bin = Binomial(n, p)

# %%
# 確率分布オブジェクトの乱数

rand(bin)

# %%
# 確率分布オブジェクトの乱数を複数個生成

@show rand(bin, 20);

# %%
# サンプルのヒストグラム
# norm = true で確率密度函数と比較できるようにする。
# alpha = 0.3 で半透明化pp
# bin = -0.5:n+0.5 で 1 刻みの適切なビンを設定

histogram(sample_of_means; norm=true, alpha=0.3, bin=-0.5:n+0.5, label="sample of means")

# %%
# -0.5:n+0.5 は [-0.5, 0.5, 1.5, 2.5, ..., 19.5, 20.5] とほぼ等価

-0.5:n+0.5

# %%
@show collect(-0.5:n+0.5);

# %%
# 平均μ、標準偏差σの正規分布オブジェクトは次にようにして作れる。

normal = Normal(n*p, √(n*p*(1-p)))

# %%
# 確率分布オブジェクト dist は plot(dist) でプロットできる。
# lw = 2 で2倍の太さでプロット

plot(normal, 0, 20; label="normal dist.", lw=2)

# %%
# 以上を重ねてプロット

histogram(sample_of_means; norm=true, alpha=0.3, bin=-0.5:n+0.5, label="rand(Binomial($n, $p))")
plot!(Normal(n*p, √(n*p*(1-p))), 0, 20; label="normal approx.", lw=2)

# %%
# サイコロを10回振って出た目の平均値の分布は正規分布でよく近似される。

n = 10
M = [mean(rand(1:6, n)) for _ in 1:10^5]
histogram(M; norm=true, alpha=0.3, bin=0.95:0.1:6.05, label="mean(rand(1:6, 10))")
plot!(Normal(mean(1:6), std(1:6; corrected=false)/√n), 1, 6; label="normal approx.", lw=2)

# %%
# 2つの確率パラメータ(比率パラメータ)が等しい二項分布で生成されたサンプルのZ統計量は
# 近似的に標準正規分布に従う。

function gensample(bin1, bin2)
    m, n = ntrials(bin1), ntrials(bin2)
    a, b = rand(bin1), rand(bin2)
    c, d = m - a, n - b
    a, b, c, d
end

function zstat(a, b, c, d)
    a*d - b*c == 0 && return zero(inv(a))
    m, n = a + c, b + d
    N = a + b + c + d
    p = (a + b)/N
    Z = (a/m - b/n) / √(p*(1-p)*(1/m + 1/n))
    Z
end

m, n, p = 30, 20, 0.3
bin1 = Binomial(m, p)
bin2 = Binomial(n, p)
Z = [zstat(gensample(bin1, bin2)...) for _ in 1:10^5]

histogram(Z; norm=true, alpha=0.3, bin=-4.25:0.5:4.25, label="Z statistics")
plot!(Normal(), -4, 4; label="std normal dist.", lw=2)

# %%
@code_warntype zstat(1, 2, 3, 4)

# %%
# 数式処理を使って、Z統計量の2乗がPearsonのχ²統計量に一致することを確認。
# 上の数値計算で使ったzstat函数をそのまま数式処理で利用する。

using SymPy

@syms a b c d
Z = zstat(a, b, c, d)

# %%
# Z^2 |> factor は factor(Z^2) に等価
# 以下の計算結果はPearsonのχ²統計量の有名な表示に一致

Z^2 |> factor

# %%
# (O - E)²/E の和で定義されたPearsonのχ²統計量

function chisqstat(a, b, c, d)
    a*d - b*c == 0 && return zero(inv(a))
    N = a + b + c + d
    s, f = a + b, c + d
    m, n = a + c, b + d
    Ea, Eb, Ec, Ed = m*s/N, n*s/N, m*f/N, n*f/N
    X² = (a - Ea)^2/Ea + (b - Eb)^2/Eb + (c - Ec)^2/Ec + (d - Ed)^2/Ed
end

X² = chisqstat(a, b, c, d)

# %%
# 以下の計算結果は2×2の分割表のPearson χ²統計量に関する有名な公式

X² |> factor

# %%
# 正確二項検定のP値の定義

# \lessapprox TAB → ⪅
# \approx TAB → ≈
x ⪅ y = x < y || x ≈ y

pval_exact(dist, k) = sum(pdf(dist, j) for j in support(dist) if pdf(dist, j) ⪅ pdf(dist, k))

# %%
# 正確二項検定のP値の例　(パラメータpを固定した場合)

n, p = 20, 0.3
bin = Binomial(n, p)

k = support(bin)
y = pval_exact.(bin, k)
plot(k, y; label="")
plot!(; xtick=0:20, ytick=0:0.1:1)
plot!(; xlabel="data k", ylabel="p-value for parameter n = $n, p = $p")

# %%
# 正確二項検定のP値の例　(データkを固定した場合)

n, k = 20, 6

p = 0:0.002:1
y = @. pval_exact(Binomial(n, p), k)
plot(p, y; label="")
plot!(; xtick=0:0.1:1, ytick=0:0.1:1)
plot!(; xlabel="parameter p", ylabel="p-value for data n = $n, k = $k")

# %%
# 信頼区間函数の定義

using Roots

function confint(pvalfunc, mlefunc, n, k, pmin, pmax; α = 0.05)
    f(p) = pvalfunc(n, p, k) - α
    ci = find_zeros(f, pmin, pmax)
    length(ci) ≥ 2 && return ci
    pmle = mlefunc(n, k)
    first(ci) > pmle ? [pmin, first(ci)] : [first(ci), pmax]
end

# %%
# 信頼区間の例

pvalfunc(n, p, k) = pval_exact(Binomial(n, p), k)
mlefunc(n, k) = k/n
pmin, pmax = 0.0, 1.0

n = 20
k = 0:20
α = 0.05
ci = confint.(pvalfunc, mlefunc, n, k, pmin, pmax; α)
ciL, ciR = first.(ci), last.(ci)

plot(; legend=:topleft)
plot!(k, ciL; label="min. of conf. int.")
plot!(k, ciR; label="max. of conf. int.")
plot!(k, ciL; label="$(100(1 - α))% confidence interval", lw=0, frange=ciR, fa=0.1)
plot!(; xtick=0:n, ytick=0:0.1:1)
plot!(; xlabel="data", ylabel="parameter")

# %%
