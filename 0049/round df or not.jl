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
#using StatsPlots
#default(fmt=:png)

function pvalue_welch_t(m, X̄, U, n, Ȳ, V; Δμ=0.0, roundfunc=identity, dffunc=degree_of_freedom)
    t = tvalue_welch(m, X̄, U, n, Ȳ, V; Δμ)
    dfhat = dffunc(m, U, n, V)
    2ccdf(TDist(roundfunc(dfhat)), abs(t))
end

pvalue_welch_t(X, Y; Δμ=0.0, roundfunc=identity, dffunc=degree_of_freedom) =
    pvalue_welch_t(length(X), mean(X), var(X), length(Y), mean(Y), var(Y); Δμ, roundfunc, dffunc)

tvalue_welch(m, X̄, U, n, Ȳ, V; Δμ=0.0) = (X̄ - Ȳ - Δμ) / √(U/m + V/n)
tvalue_welch(X, Y; Δμ=0.0) =
    tvalue_welch(length(X), mean(X), var(X), length(Y), mean(Y), var(Y); Δμ)

degree_of_freedom(m, U, n, V) = (U/m + V/n)^2 / ((U/m)^2/(m-1) + (V/n)^2/(n-1))
degree_of_freedom(X, Y) = degree_of_freedom(length(X), var(X), length(Y), var(Y))

degree_of_freedom_not_Satterthwaite(m, U, n, V) =
    (U/m + V/n)^2 / ((U/m)^2/(m+1) + (V/n)^2/(n+1)) - 2
degree_of_freedom_not_Satterthwaite(X, Y) =
    degree_of_freedom_not_Satterthwaite(length(X), var(X), length(Y), var(Y))

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

# %% [markdown]
# ## R

# %%
using RCall
X = [80, 87, 78, 72]
Y = [119, 78, 95, 124, 85, 92]
@rput X Y
R"""t.test(X, Y)"""

# %%
@show tvalue_welch(X, Y);
@show degree_of_freedom(X, Y);
@show pvalue_welch_t(X, Y);

# %% [markdown]
# ## STATA

# %% [markdown]
# https://www.statology.org/welchs-t-test-stata/
#
# <img src="IMG_5204.png">

# %%
m, X̄, U = 12, 21.00, 2.730301^2
n, Ȳ, V = 12, 22.75, 3.250874^2
@show (; N=m, Mean=X̄, Std=√U, SE=√(U/m));
@show (; N=n, Mean=Ȳ, Std=√V, SE=√(V/n));
@show tvalue_welch(m, X̄, U, n, Ȳ, V);
println()
@show degree_of_freedom(m, U, n, V);
@show pvalue_welch_t(m, X̄, U, n, Ȳ, V);
println()
@show degree_of_freedom_not_Satterthwaite(m, U, n, V);
@show pvalue_welch_t(m, X̄, U, n, Ȳ, V; dffunc=degree_of_freedom_not_Satterthwaite);

# %%
using SymPy
@syms U V m n
expr = m + n - degree_of_freedom_not_Satterthwaite(m, U, n, V)
factor(expr)

# %%
using SymPy
@syms U V m n
expr = m + n - 2 - degree_of_freedom(m, U, n, V)
factor(expr)

# %%
using SymPy
@syms U V m n
expr = m + n - (U/m + V/n)^2 / ((U/m)^2/m + (V/n)^2/n)
factor(expr)

# %% [markdown]
# https://stats.oarc.ucla.edu/stata/output/t-test/
#
# <img src="IMG_5209.jpeg" width=60%>

# %%
m, X̄, U =  91, 50.12088, 10.30516^2
n, Ȳ, V = 109, 54.99083, 8.133715^2
@show (; N=m, Mean=X̄, Std=√U, SE=√(U/m));
@show (; N=n, Mean=Ȳ, Std=√V, SE=√(V/n));
@show tvalue_welch(m, X̄, U, n, Ȳ, V);
@show degree_of_freedom(m, U, n, V);
@show pvalue_welch_t(m, X̄, U, n, Ȳ, V);

# %% [markdown]
# ## SAS

# %% [markdown]
# https://support.sas.com/documentation/onlinedoc/stat/132/ttest.pdf
#
# <img src="IMG_5210.jpeg" width=60%>
# <img src="IMG_5211.jpeg" width=50%>

# %%
m, X̄, U = 7, 76.8571, 2.5448^2
n, Ȳ, V = 7, 82.7143, 3.1472^2
@show (; N=m, Mean=X̄, Std=√U, SE=√(U/m));
@show (; N=n, Mean=Ȳ, Std=√V, SE=√(V/n));
@show tvalue_welch(m, X̄, U, n, Ȳ, V);
println()
@show degree_of_freedom(m, U, n, V);
@show pvalue_welch_t(m, X̄, U, n, Ȳ, V);

# %% [markdown]
# ## SPSS

# %% [markdown]
# https://www.stats-guild.com/analytics/15678
#
# <img src="IMG_5212.jpeg" width=100%>

# %%
m, X̄, U = 115, 3054.17, 751.684^2
n, Ȳ, V =  74, 2773.16, 660.340^2
@show (; N=m, Mean=X̄, Std=√U, SE=√(U/m));
@show (; N=n, Mean=Ȳ, Std=√V, SE=√(V/n));
@show tvalue_welch(m, X̄, U, n, Ȳ, V);
@show degree_of_freedom(m, U, n, V);
@show pvalue_welch_t(m, X̄, U, n, Ȳ, V);

# %%
