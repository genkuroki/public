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
#     display_name: Julia 1.7.3
#     language: julia
#     name: julia-1.7
# ---

# %%
using Combinatorics
using Distributions
using StatsBase: countmap
using StatsPlots
using SymPy

# %%
n, k = 21, 10

# %%
rank_sums = [sum(comb) for comb in combinations(1:n, k)]
c = countmap(rank_sums)
xmin, xmax = extrema(keys(c))
@show xs = xmin:xmax
ys = (k -> c[k]).(xs) / length(rank_sums)
plot(xs, ys; label="dist. of rank sums")

# %%
@vars q x
F = sympy.Poly(prod(1 + q^i*x for i in 1:n), x)
f = F.coeffs()
f_k = sympy.Poly(f[n+1 - k], q)
Ys = float(f_k.coeffs())
Ys ./= sum(Ys)
Xmax = f_k.degree()
Xmin = Xmax - length(Ys) + 1
@show Xs = Xmin:Xmax
plot(Xs, Ys; label="dist. of coeffs. of q-binom")

# %%
xs == Xs

# %%
ys == Ys

# %%
mu = mean(rank_sums)

# %%
k*(n-k)/2 + k*(k+1)/2 == mu

# %%
s2 = var(rank_sums; corrected=false)

# %%
k*(n-k)*(n+1)/12 == s2

# %%
plot(Xs, Ys; label="dist. of coeffs. of q-binom")
plot!(Normal(mu, √s2); label="normal approx.", ls=:dash)

# %%
f_k

# %%
# (x)_q = (1 - q^x)/(1 - q)
# (n)_q! = (1)_q (2)_q … (n)_q
# g_k = q^{k(k+1)/2} (n)_q!/((k)_q! (n-k)_q!), essentially a q-binomial coefficient
g_k = sympy.Poly(q^(k*(k+1)÷2) * prod((1 - q^(n-k+j))/(1 - q^j) for j in 1:k; init=Sym(1)).factor(), q)

# %%
f_k == g_k

# %%
