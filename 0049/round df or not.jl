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

# %%
using Distributions
using StatsPlots
default(fmt=:png)

pvalue_welch_t(X, Y; Δμ=0.0, roundfunc=identity) =
    pvalue_welch_t(length(X), mean(X), var(X), length(Y), mean(Y), var(Y); Δμ, roundfunc)

function pvalue_welch_t(m, X̄, U, n, Ȳ, V; Δμ=0.0, roundfunc=round)
    t = (X̄ - Ȳ - Δμ) / √(U/m + V/n)
    df = degree_of_freedom(m, U, n, V)
    2ccdf(TDist(roundfunc(df)), abs(t))
end

function degree_of_freedom(m, U, n, V)
    (U/m + V/n)^2 / ((U/m)^2/(m-1) + (V/n)^2/(n-1))
end

degree_of_freedom(X, Y) = degree_of_freedom(length(X), var(X), length(Y), var(Y))

X = [80, 87, 78, 72]
Y = [119, 78, 95, 124, 85, 92]
@show X Y
@show degree_of_freedom(X, Y);
@show pvalue_welch_t(X, Y);
@show pvalue_welch_t(X, Y; roundfunc=round);

# %%
X = [8, 21, 22, 30]
Y = [5, 42, 43, 83, 83, 119]
@show X Y
@show degree_of_freedom(X, Y);
@show pvalue_welch_t(X, Y);
@show pvalue_welch_t(X, Y; roundfunc=round);

# %%
X = sort(round.(Int, 20*[-0.8392321454441508, -0.13632577223881945, -0.1907178263261891, 0.25603978994522353] .+ 25))
Y = sort(round.(Int, 20*[-0.9932984206053254, 2.913878236092646, 2.8891275856381764, 4.716892383695555, 0.8317174086782941, 0.8766881245471467] .+ 25))
@show X Y
@show degree_of_freedom(X, Y);
@show pvalue_welch_t(X, Y);
@show pvalue_welch_t(X, Y; roundfunc=round);

# %%
X = [-0.8392321454441508, -0.13632577223881945, -0.1907178263261891, 0.25603978994522353]
Y = [-0.9932984206053254, 2.913878236092646, 2.8891275856381764, 4.716892383695555, 0.8317174086782941, 0.8766881245471467]
@show X Y
@show degree_of_freedom(X, Y);
@show pvalue_welch_t(X, Y);
@show pvalue_welch_t(X, Y; roundfunc=round);

# %%
X = [0.8115162043087973, -0.8404363548455143, 1.6727606354453446, 0.533847620210447]
Y = [1.6339095622225246, 2.065488181473473, 2.0345943286959765, 1.7579073472598687, 3.3451457705490997, 1.8267931745140344]
@show X Y
@show degree_of_freedom(X, Y);
@show pvalue_welch_t(X, Y);
@show pvalue_welch_t(X, Y; roundfunc=round);

# %%
# using HypothesisTests

for _ in 1:10^6
    X = round.(Int, rand(Normal(80, 10), 4))
    Y = round.(Int, rand(Normal(80, 20), 6))
    abs(rem(degree_of_freedom(X, Y), 1) - 0.5) ≥ 0.05 && continue
    (pvalue_welch_t(X, Y) - 0.05) * (pvalue_welch_t(X, Y; roundfunc=round) - 0.05) < 0 && break
end

@show X Y
@show degree_of_freedom(X, Y);
@show pvalue_welch_t(X, Y);
@show pvalue_welch_t(X, Y; roundfunc=round);
# UnequalVarianceTTest(X, Y)

# %%
