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

# %% [markdown]
# * https://twitter.com/eggplantmed/status/1629147688053923840
#   * https://www.pref.chiba.lg.jp/kyouiku/anzen/hokenn/covid-19.html
#     * https://www.pref.chiba.lg.jp/kyouiku/anzen/hokenn/documents/mokushoku-minaoshi.pdf

# %%
using Distributions
using Random
using StatsBase: ecdf
using StatsPlots
default(fmt=:png)

# %%
safediv(x, y) = x==0 ? zero(x/y) : x/y

function stat_chisq(a, b, c, d)
    safediv((a+b+c+d)*(a*d - b*c)^2, (a+b)*(c+d)*(a+c)*(b+d))
end

function pvalue_chisq(a, b, c, d)
    χ² = stat_chisq(a, b, c, d)
    ccdf(Chisq(1), χ²)
end

# %%
function sim_chisq(p, q, r, N; L=10^5)
    mult = Multinomial(N, [p*r, (1-p)*r, q*(1-r), (1-q)*(1-r)])
    pval = Vector{Float64}(undef, L)
    tmpA = [Vector{Int}(undef, 4) for _ in 1:Threads.nthreads()]
    Threads.@threads for i in 1:L
        A = rand!(mult, tmpA[Threads.threadid()])
        pval[i] = pvalue_chisq(A...)
    end
    _ecdf_pval = ecdf(pval)
    ecdf_pval(x) = _ecdf_pval(x)
    ecdf_pval
end

# %%
function plot_sim(p, q, r, N; L=10^5)
    ecdf_pval = sim_chisq(p, q, r, N; L)
    plot(ecdf_pval, 0, 1; label="power")
    plot!(identity; label="", ls=:dot)
    plot!(xtick=0:0.05:1, ytick=0:0.05:1, xrotation=90)
    plot!(size=(400, 400))
end

# %%
# test
plot_sim(0.012, 0.012, 1/5, 2500)

# %%
# test
plot_sim(0.009, 0.014, 1/5, 2500)

# %%
# test
plot_sim(0.009, 0.014, 1/5, 23000)

# %%
plot_sim(5/551, 27/1950, 551/(551+1950), 551+1950)

# %%
plot_sim(5/551, 27/1950, 551/(551+1950), round(Int, 9.5*(551+1950)))

# %%
