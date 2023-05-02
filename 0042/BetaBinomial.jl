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
using Distributions
using QuadGK
using StatsPlots
default(fmt=:png)

# %%
function pdf_bb2d(bb::BetaBinomial, x, p)
    n, α, β = params(bb)
    pdf(Binomial(n, p), round(Int, x)) * pdf(Beta(α, β), p)
end

function pdf_bb(bb::BetaBinomial, x)
    quadgk(p -> pdf_bb2d(bb, x, p), 0, 1)[1]
end

# %%
n, α, β = 10, 0.5, 0.5
bb = BetaBinomial(n, α, β)
x = range(-1, n+1, 1001)
plot(x, x -> pdf(bb, round(Int, x)); label="BetaBinomial($n, $α, $β)")
plot!(x, x -> pdf_bb(bb, x); label="", ls=:dash)
#plot!(x, x -> pdf(Binomial(n, α/(α+β)), round(Int, x)); label="Binomial($n, $α/($α+$β)")
plot!(xtick=0:n)

# %%
n, α, β = 10, 5, 5
bb = BetaBinomial(n, α, β)
x = range(-0.5, n+0.5, 201)
p = range(0, 1, 201)
heatmap(x, p, (x, p) -> pdf_bb2d(bb, x, p))
plot!(xtick=0:n, ytick=0:0.1:1)

# %%
n, α, β = 10, 5, 5
bb = BetaBinomial(n, α, β)
x = range(-0.5, n+0.5, 201)
p = range(0, 1, 201)
surface(x, p, (x, p) -> pdf_bb2d(bb, x, p))
plot!(xtick=0:n, ytick=0:0.1:1)

# %%
n, α, β = 10, 5, 5
bb = BetaBinomial(n, α, β)
x = range(-1, n+1, 1001)
plot(x, x -> pdf(bb, round(Int, x)); label="BetaBinomial($n, $α, $β)")
plot!(x, x -> pdf(Binomial(n, α/(α+β)), round(Int, x)); label="Binomial($n, $α/($α+$β)")
plot!(xtick=0:n)

# %%
n, α, β = 100, 5, 5
bb = BetaBinomial(n, α, β)
x = range(-1, n+1, 1001)
plot(x, x -> pdf(bb, round(Int, x)); label="BetaBinomial($n, $α, $β)")
plot!(x, x -> pdf(Binomial(n, α/(α+β)), round(Int, x)); label="Binomial($n, $α/($α+$β)")
plot!(xtick=0:10:n)

# %%
n, α, β = 10, 0.5, 0.5
bb = BetaBinomial(n, α, β)
x = range(-0.5, n+0.5, 201)
p = range(0, 1, 201)
heatmap(x, p, (x, p) -> pdf(Binomial(n, p), round(Int, x)))
plot!(xtick=0:n, ytick=0:0.1:1)

# %%
n, α, β = 10, 0.5, 0.5
bb = BetaBinomial(n, α, β)
x = range(-0.5, n+0.5, 201)
p = range(0, 1, 201)
surface(x, p, (x, p) -> pdf(Binomial(n, p), round(Int, x)))
plot!(xtick=0:n, ytick=0:0.1:1)

# %%
n, α, β = 10, 0.5, 0.5
bb = BetaBinomial(n, α, β)
x = range(-1, n+1, 1001)
plot(x, x -> pdf(bb, round(Int, x)); label="BetaBinomial($n, $α, $β)")
plot!(x, x -> pdf(Binomial(n, α/(α+β)), round(Int, x)); label="Binomial($n, $α/($α+$β)")
plot!(xtick=0:n)

# %%
n, α, β = 100, 0.5, 0.5
bb = BetaBinomial(n, α, β)
x = range(-1, n+1, 1001)
plot(x, x -> pdf(bb, round(Int, x)); label="BetaBinomial($n, $α, $β)")
plot!(x, x -> pdf(Binomial(n, α/(α+β)), round(Int, x)); label="Binomial($n, $α/($α+$β)")
plot!(xtick=0:10:n)

# %%
