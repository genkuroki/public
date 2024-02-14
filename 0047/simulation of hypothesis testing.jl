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

# %%
using Distributions

function ci_wald(n, k; α=0.05)
    c = cquantile(Normal(), α/2)
    p̂ = k/n
    sehat = √(p̂*(1-p̂)/n)
    p̂ - c*sehat, p̂ + c*sehat
end

"""
kは確率1-qで二項分布Binomial(n, p0)に従う乱数であり,
確率qで二項分布Binomial(n, p1)に従う乱数であるとする.
n, kからWaldの100(1-α)%信頼区間[L, U]を計算する.
以上をniters回実行する.

条件A, Bの定義:

* A = (L ≤ p0 ≤ U)
* B = (L ≤ p1 ≤ U)

kが二項分布Binomial(n, p0)に従う乱数の場合について

* a = AかつBとなった回数
* b = Aかつnot Bとなった回数
* c = not AかつBとなった回数
* d = not Aかつnot Bとなった回数

kが二項分布Binomial(n, p1)に従う乱数の場合について

* e = AかつBとなった回数
* f = Aかつnot Bとなった回数
* g = not AかつBとなった回数
* h = not Aかつnot Bとなった回数

このとき

* g/(c+d+g+h)はp0が信頼区間に含まれない場合にkの生成に使われたp=p0, p1が信頼区間に含まれる条件付き確率の近似値になる.
* (a+b+e)/(a+b+e+f)はp0が信頼区間に含まれる場合にkの生成に使われたp=p0, p1が信頼区間に含まれる条件付き確率の近似値になる.
"""
function sim(; p0=0.5, p1=0.6, n=200, q=1e-4, α=0.05, niters=10^8)
    a, b, c, d, e, f, g, h = 0, 0, 0, 0, 0, 0, 0, 0
    for i in 1:niters
        if rand() > q
            k = rand(Binomial(n, p0))
            L, U = ci_wald(n, k; α)
            A = L ≤ p0 ≤ U
            B = L ≤ p1 ≤ U
            a += A & B
            b += A & !B
            c += !A & B
            d += !A & !B
        else
            k = rand(Binomial(n, p1))
            L, U = ci_wald(n, k; α)
            A = L ≤ p0 ≤ U
            B = L ≤ p1 ≤ U
            e += A & B
            f += A & !B
            g += !A & B
            h += !A & !B
        end
    end
    a, b, c, d, e, f, g, h
end

# %%
@doc sim

# %%
@show a, b, c, d, e, f, g, h = sim(p0=0.5, p1=0.6, n=200, q=1e-4, α=0.05)
@show g/(c+d+g+h) (a+b+e)/(a+b+e+f) (a+b)/(a+b+c+d) (g+h)/(e+f+g+h);

# %%
@show a, b, c, d, e, f, g, h = sim(p0=0.5, p1=0.6, n=200, q=0.001, α=0.05)
@show g/(c+d+g+h) (a+b+e)/(a+b+e+f) (a+b)/(a+b+c+d) (g+h)/(e+f+g+h);

# %%
@show a, b, c, d, e, f, g, h = sim(p0=0.5, p1=0.6, n=200, q=0.01, α=0.05)
@show g/(c+d+g+h) (a+b+e)/(a+b+e+f) (a+b)/(a+b+c+d) (g+h)/(e+f+g+h);

# %%
@show a, b, c, d, e, f, g, h = sim(p0=0.5, p1=0.6, n=300, q=0.01, α=0.01)
@show g/(c+d+g+h) (a+b+e)/(a+b+e+f) (a+b)/(a+b+c+d) (g+h)/(e+f+g+h);

# %%
@show a, b, c, d, e, f, g, h = sim(p0=0.5, p1=0.6, n=420, q=0.01, α=0.001)
@show g/(c+d+g+h) (a+b+e)/(a+b+e+f) (a+b)/(a+b+c+d) (g+h)/(e+f+g+h);

# %%
@show a, b, c, d, e, f, g, h = sim(p0=0.5, p1=0.6, n=480, q=0.01, α=0.0004)
@show g/(c+d+g+h) (a+b+e)/(a+b+e+f) (a+b)/(a+b+c+d) (g+h)/(e+f+g+h);

# %%
@show a, b, c, d, e, f, g, h = sim(p0=0.5, p1=0.6, n=200, q=0.5, α=0.05)
@show g/(c+d+g+h) (a+b+e)/(a+b+e+f) (a+b)/(a+b+c+d) (g+h)/(e+f+g+h);

# %%
@show a, b, c, d, e, f, g, h = sim(p0=0.5, p1=0.6, n=200, q=0.8, α=0.05)
@show g/(c+d+g+h) (a+b+e)/(a+b+e+f) (a+b)/(a+b+c+d) (g+h)/(e+f+g+h);

# %%
