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
#     display_name: Julia 1.11.6
#     language: julia
#     name: julia-1.11
# ---

# %% [markdown]
# # Agresti-Coullの信頼区間
#
# * 黒木玄
# * Created at 2025-02-21, Updated at 2025-08-10
#
# 有名なAgresti-Coullの信頼区間は論文
#
# * A. Agresti and B.A. Coull, Approximate is better than “exact” for interval estimation of binomial proportions, The American Statistician, 1998.  https://scholar.google.co.jp/scholar?cluster=5129299358902170657
#
# で提案された。Agresti-Coullの信頼区間両端の値は次の公式で計算される:
#
# $$
# p^{\text{AC}}_\pm = \tilde{p} + z \sqrt{\frac{\tilde{p}(1-\tilde{p})}{\tilde{n}}}, \quad
# z = \operatorname{cquantile}\left(\operatorname{Normal}(0,1), \frac{\alpha}{2}\right), \quad
# \tilde{k}=k+\frac{z^2}{2}, \quad
# \tilde{n}=n+z^2, \quad
# \tilde{p}=\frac{\tilde{k}}{\tilde{n}}.
# $$
#
# Agresti-Coullの信頼区間はデータの値 $(k, n-k)$ を $\left(k+\frac{z^2}{2}, n-k+\frac{z^2}{2}\right)$ に修正した後のWaldの信頼区間に等しい.
#
# このノートでは以下の2つを実行する:
#
# * Agresti-Coullの信頼区間がWilsonのスコア信頼区間をよく近似していることの確認
# * Agresti-Coullの信頼区間を与えるP値関数の実装

# %% [markdown]
# <img src="https://raw.githubusercontent.com/genkuroki/public/refs/heads/main/0053/Agresti-Coull-1.png">

# %% [markdown]
# <img src="https://raw.githubusercontent.com/genkuroki/public/refs/heads/main/0053/Agresti-Coull-2.png">

# %%
using Distributions
using LaTeXStrings
using Roots
using StatsPlots
default(fmt=:png, legendfontsize=12, guidefontsize=12)
using SymPy

# %%
function confint_score(k, n, α = 0.05)
    z = cquantile(Normal(), α/2)
    p̂ = k/n
    sehat² = p̂ * (1 - p̂) / n
    a, b = 1 + z^2/n, p̂ + z^2/(2n)
    sqrtD = z * √(sehat² + z^2/(4n^2))
    p_L = (b - sqrtD) / a
    p_U = (b + sqrtD) / a
    [p_L, p_U]
end

function confint_wald(k, n, α = 0.05)
    z = cquantile(Normal(), α/2)
    p̂ = k/n
    sehat = √(p̂ * (1 - p̂) / n)
    p_L = p̂ - z*sehat
    p_U = p̂ + z*sehat
    [p_L, p_U]
end

function confint_agresti_coull(k, n, α = 0.05)
    z = cquantile(Normal(), α/2)
    k̃, ñ = k + z^2/2, n + z^2
    confint_wald(k̃, ñ, α)
end

# %%
safediv(x, y) = x == 0 ? zero(x/y) : x/y

function pvalue_score(k, n, p)
    0 ≤ p ≤ 1 || return zero(p)
    bin = Binomial(n, p)
    z = safediv(k - mean(bin), std(bin))
    2ccdf(Normal(), abs(z))
end

function pvalue_wald(k, n, p)
    p̂ = k/n
    z = safediv(k - n*p, √(n*p̂*(1-p̂)))
    2ccdf(Normal(), abs(z))
end

# z²に関する3次方程式
function _eq_agresti_coull(z², k, n, p)
    k̃, ñ = k + z²/2, n + z²
    ñ * (k̃ - ñ*p)^2 - z² * k̃ * (ñ - k̃)
end

# z²に関する3次方程式を解く
function z²_agresti_coull(k, n, p)
    f(t) = _eq_agresti_coull(exp(t), k, n, p)
    (1 - p ≈ 1 || p ≈ 1) && return oftype(p, Inf)
    k ≈ n*p && return zero(p)
    z² = exp(find_zero(f, (-1e2, 1e2))) # 手抜き
end

function pvalue_agresti_coull(k, n, p)
    0 ≤ p ≤ 1 || return zero(p)
    2ccdf(Normal(), √z²_agresti_coull(k, n, p))
end

# %%
@syms x, k, n, p
print("cubic equation for x = z²: ")
Eq(_eq_agresti_coull(x, k, n, p), 0) |> display
println()

@show k, n, p = 3, 10, 0.6
@show pvalue_score(k, n, p) pvalue_agresti_coull(k, n, p) pvalue_wald(k, n, p)
println()
stack(([p, pvalue_score(k, n, p), pvalue_agresti_coull(k, n, p), pvalue_wald(k, n, p)] for p in 0:0.1:1); dims=1) |> display
println()

pgfplotsx()
#gr()
plot(z² -> _eq_agresti_coull(z², k, n, p), -17, 5; label="")
scatter!([z²_agresti_coull(k, n, p)], [0.0]; label="positive solution")
hline!([0]; label="", c=:red)
vline!([0]; label="", c=:red)
plot!(legend=:top)
plot!(xguide=L"x = z^2")
title!(L"cubic equation for $x = z^2$: $n=%$n$, $k=%$k$, $p=%$p$")
plot!(size=(500, 350)) |> display

# %%
pgfplotsx()
#gr()
k, n = 3, 10

r(x) = round(x; sigdigits=3)
@show confint_score(k, n) .|> r
@show confint_agresti_coull(k, n) .|> r
@show confint_wald(k, n) .|> r

plot()
plot!(p -> pvalue_score(k, n, p), -0.2, 1; label="Wilson score", c=1)
plot!(confint_score(k, n), fill(0.06, 2); label="", c=1)
plot!(p -> pvalue_agresti_coull(k, n, p); label="Agresti-Coull", ls=:dash, c=2)
plot!(confint_agresti_coull(k, n), fill(0.05, 2); label="", c=2)
plot!(p -> pvalue_wald(k, n, p); label="Wald", ls=:dot, c=3)
plot!(confint_wald(k, n), fill(0.04, 2); label="", c=3)
plot!(legend=:topright)
title!(L"P-value function for $n=%$n$, $k=%$k$")
plot!(xtick=-0.2:0.1:1.2, ytick=0:0.1:1)
plot!(xguide=L"p", yguide="P-value")
plot!(size=(500, 350))

# %%
pgfplotsx()
#gr()
k, n, = 7, 10

r(x) = round(x; sigdigits=3)
@show confint_score(k, n) .|> r
@show confint_agresti_coull(k, n) .|> r
@show confint_wald(k, n) .|> r

plot()
plot!(p -> pvalue_score(k, n, p), 0, 1.2; label="Wilson score", c=1)
plot!(confint_score(k, n), fill(0.06, 2); label="", c=1)
plot!(p -> pvalue_agresti_coull(k, n, p); label="Agresti-Coull", ls=:dash, c=2)
plot!(confint_agresti_coull(k, n), fill(0.05, 2); label="", c=2)
plot!(p -> pvalue_wald(k, n, p); label="Wald", ls=:dot, c=3)
plot!(confint_wald(k, n), fill(0.04, 2); label="", c=3)
plot!(legend=:topleft)
title!(L"P-value function for $n=%$n$, $k=%$k$")
plot!(xtick=-0.2:0.1:1.2, ytick=0:0.1:1)
plot!(xguide=L"p", yguide="P-value")
plot!(size=(500, 350))

# %%
pgfplotsx()
#gr()
k, n, = 70, 100

r(x) = round(x; sigdigits=3)
@show confint_score(k, n) .|> r
@show confint_agresti_coull(k, n) .|> r
@show confint_wald(k, n) .|> r

plot()
plot!(p -> pvalue_score(k, n, p), 0, 1; label="Wilson score", c=1)
plot!(confint_score(k, n), fill(0.06, 2); label="", c=1)
plot!(p -> pvalue_agresti_coull(k, n, p); label="Agresti-Coull", ls=:dash, c=2)
plot!(confint_agresti_coull(k, n), fill(0.05, 2); label="", c=2)
plot!(p -> pvalue_wald(k, n, p); label="Wald", ls=:dot, c=3)
plot!(confint_wald(k, n), fill(0.04, 2); label="", c=3)
plot!(legend=:topleft)
title!(L"P-value function for $n=%$n$, $k=%$k$")
plot!(xtick=-0.2:0.1:1.2, ytick=0:0.1:1)
plot!(xguide=L"p", yguide="P-value")
plot!(size=(500, 350))

# %%
pgfplotsx()
#gr()
k, n = 3, 10

P = plot()
plot!(p -> pvalue_score(k, n, p), -0.2, 1; label="", c=1)
for α in [0.001; 0.01:0.01:1]
    ci = confint_score(k, n, α)
    plot!(ci, fill(α, 2); label=(α==1 ? "CIs" : ""), c=4)
end
title!(L"Wilson score: $n=%$n$, $k=%$k$")

Q = plot()
plot!(p -> pvalue_agresti_coull(k, n, p), -0.2, 1; label="", c=2)
for α in [0.001; 0.01:0.01:1]
    ci = confint_agresti_coull(k, n, α)
    plot!(ci, fill(α, 2); label=(α==1 ? "CIs" : ""), c=5)
end
title!(L"Agresti-Coull: $n=%$n$, $k=%$k$")

R = plot()
plot!(p -> pvalue_wald(k, n, p), -0.2, 1; label="", c=3)
for α in [0.001; 0.01:0.01:1]
    ci = confint_wald(k, n, α)
    plot!(ci, fill(α, 2);  label=(α==1 ? "CIs" : ""), c=6)
end
title!(L"Wald: $n=%$n$, $k=%$k$")

plot(P, Q, R; size=(1000, 700), layout=(2, 2))
plot!(xtick=-0.2:0.1:1, ytick=0:0.1:1)
plot!(xguide=L"p", yguide="P-value")
plot!(legend=:topright)
#plot!(leftmargin=4Plots.mm, bottommargin=4Plots.mm)p

# %%
