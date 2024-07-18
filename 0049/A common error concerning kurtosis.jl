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
#     display_name: Julia 1.10.4
#     language: julia
#     name: julia-1.10
# ---

# %% [markdown]
# * A Common Error Concerning Kurtosis
# * Irving Kaplansky
# * Journal of the American Statistical Association, Vol. 40, No. 230 (Jun., 1945), p. 259
# * https://psycnet.apa.org/doi/10.2307/2280139
#
# <img src="IMG_5004.jpeg" width=70%>

# %% tags=[]
using Distributions
using QuadGK
using StatsPlots
default(fmt=:png)

meanf(f) = quadgk(x -> x*f(x), -Inf, Inf)[1]

function varf(f)
    μ = meanf(f)
    quadgk(x -> (x-μ)^2*f(x), -Inf, Inf)[1]
end

function muf(f, k)
    μ = meanf(f)
    σ = √varf(f)
    quadgk(x -> (x-μ)^k*f(x), -Inf, Inf)[1]/σ^k
end

kurtosisf(f) = muf(f, 4) - 3

P(x) = 1/(3√π) * (9/4 + x^4)*exp(-x^2)
Q(x) = 3/(2√(2π)) * exp(-x^2/2) - 1/2 * P(x)
R(x) = 1/(6√π) * (exp(-x^2/4) + 4exp(-x^2))
S(x) = 3√3/(16√π) * (2 + x^2) * exp(-3x^2/4)

for f in (P, Q, R, S)
    @eval (
        f = $f,
        #one = round(quadgk($f, -Inf, Inf)[1]; digits=3),
        μ = round(meanf($f); digits=3),
        σ² = round(varf($f); digits=3),
        μ₄ = round(muf($f, 4); digits=3),
        f₀ = round($f(0); digits=3),
        η = round(kurtosisf($f); digits=3), 
    ) |> println
end

plot(
    plot(P; label="P(x)", c=1), plot(Q; label="Q(x)", c=2),
    plot(R; label="R(x)", c=3), plot(S; label="S(x)", c=4);
    size=(800, 600), layout=(2,2)
)

# %%
function makef(dist, w=1/2)
    g(x) = (1 - w) * pdf(dist, x)
    g₀ = g(0)
    a = w/(2g₀)
    _f(x) = abs(x) < a ? g₀ : g(sign(x)*(abs(x) - a))
    _f
end

f = makef(Normal(), 0.1)
@show quadgk(f, -Inf, Inf)
@show kurtosisf(f)
plot(f, -5, 5)

# %%
f = makef(TDist(4), 0.8)
@show quadgk(f, -Inf, Inf)
@show kurtosisf(f)
plot(f, -10, 10; label="")

# %%
f = makef(TDist(4.5), 0.0)
@show quadgk(f, -Inf, Inf)
@show kurtosisf(f)
plot(f, -10, 10; label="")
