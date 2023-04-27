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
using Plots
default(fmt=:png)
using SpecialFunctions
fact(x) = exp(logfactorial(x))
binom(x, k) = exp(logfactorial(x) - logfactorial(k) - logfactorial(x-k))

f(x, m, n=m) = sum(binom(m, k)/binom(m+n, k) * x^k/fact(k) for k in 0:m)
pade_exp(x, m, n=m) = f(x, m, n)/ f(-x, m, n)
taylor_exp(x, m) = sum(x^k/fact(k) for k in 0:m)

# %%
a, b = -10, 10
x = range(a, b, 1000)
safe_log(x) = x ≥ 0 ? log(x) : oftype(x, NaN)

anim = @animate for m in 1:20
    P = plot(x, exp; label="exp(x)", ylim=(-10, 1000))
    plot!(x, x -> pade_exp(x, m); label="pade_exp(x, $m)", ls=:dash)
    plot!(x, x -> taylor_exp(x, 2m); label="taylor_exp(x, $(2m))", ls=:dashdot)
    title!("m = $m")
    
    Q = plot(x, identity; label="", ylim=(a, b))
    plot!(x, x -> safe_log(pade_exp(x, m)); label="log(pade_exp(x, $m))", ls=:dash)
    plot!(x, x -> safe_log(taylor_exp(x, 2m)); label="log(taylor_exp(x, $(2m)))", ls=:dashdot)
    title!("m = $m")
    
    plot(P, Q; size=(600, 600), layout=(2, 1), legend=:topleft)
end

gif(anim, "pade_exp.gif", fps=2)

# %%
