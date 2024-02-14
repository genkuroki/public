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
using BenchmarkTools
using Distributions

function pvalue_one_sided_naive(k, n, p=1/2)
    pval = 0.0
    for i in k:n
        pval += binomial(n, i) * p^i * (1-p)^(n-i)
    end
    pval
end

pvalue_one_sided_ccdf(k, n, p=1/2) = ccdf(Binomial(n, p), k-1)

@show pvalue_one_sided_naive(21, 30, 0.5)
@show pvalue_one_sided_ccdf(21, 30, 0.5);

# %%
@btime pvalue_one_sided_naive(21, 30, 0.5)
@btime pvalue_one_sided_ccdf(21, 30, 0.5);

# %%
@btime pvalue_one_sided_naive(16, 30, 0.5)
@btime pvalue_one_sided_ccdf(16, 30, 0.5);

# %%
pvalue_one_sided_naive(60, 100, 0.5)

# %%
using SpecialFunctions

function pvalue_one_sided_sf1(k, n, p=1/2)
    pval = 0.0
    for i in k:n
        pval += exp(logabsbinomial(n, i)[1]) * p^i * (1-p)^(n-i)
    end
    pval
end

@show pvalue_one_sided_sf1(60, 100, 0.5)
@show pvalue_one_sided_ccdf(60, 100, 0.5);

# %%
@btime pvalue_one_sided_sf1(60, 100, 0.5)
@btime pvalue_one_sided_ccdf(60, 100, 0.5);

# %%
pvalue_one_sided_sf1(501000, 10^6, 0.5)

# %%
using SpecialFunctions

function pvalue_one_sided_sf2(k, n, p=1/2)
    pval = 0.0
    for i in k:n
        pval += exp(logabsbinomial(n, i)[1] + i*log(p) + (n-i)*log(1-p))
    end
    pval
end

@show pvalue_one_sided_sf1(60, 100, 0.5)
@show pvalue_one_sided_sf2(60, 100, 0.5)
@show pvalue_one_sided_ccdf(60, 100, 0.5);

# %%
@btime pvalue_one_sided_sf1(60, 100, 0.5)
@btime pvalue_one_sided_sf2(60, 100, 0.5)
@btime pvalue_one_sided_ccdf(60, 100, 0.5);

# %%
@show pvalue_one_sided_sf2(501000, 10^6, 0.5)
@show pvalue_one_sided_ccdf(501000, 10^6, 0.5);

# %%
@btime pvalue_one_sided_sf2(501000, 10^6, 0.5)
@btime pvalue_one_sided_ccdf(501000, 10^6, 0.5);

# %%
ccdf(Binomial(30, 0.5), 20)

# %%
@which ccdf(Binomial(30, 0.5), 20)

# %%
Distributions.binomccdf(30, 0.5, 20)

# %%
@which Distributions.binomccdf(30, 0.5, 20)

# %%
Distributions.betacdf(20+1, 30-20, 0.5)

# %%
@which Distributions.betacdf(20+1, 30-20, 0.5)

# %%
beta_inc(20+1, 30-20, 0.5)

# %%
@which beta_inc(20+1, 30-20, 0.5)

# %%
SpecialFunctions._beta_inc(float(20+1), float(30-20), 0.5)

# %%
@which SpecialFunctions._beta_inc(float(20+1), float(30-20), 0.5)

# %% [markdown]
# https://github.com/JuliaMath/SpecialFunctions.jl/blob/master/src/beta_inc.jl#L738

# %%
