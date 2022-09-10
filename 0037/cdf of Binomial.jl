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
#     display_name: Julia 1.8.0
#     language: julia
#     name: julia-1.8
# ---

# %%
using BenchmarkTools
using Distributions
using StatsFuns: logaddexp

# %%
?logaddexp

# %%
function mycdf_naive(bin::Binomial{T}, k) where T
    n, p = Distributions.params(bin)
    odds = p / (1 - p)
    s = c = (1 - p)^n
    for i in 1:k
        c = c * (n-i+1) / i * odds
        s = s + c
    end
    s
end

function mycdf(bin::Binomial{T}, k) where T
    n, p = Distributions.params(bin)
    logodds = log(p) - log(1-p)
    logc = n*log(1-p)
    s = exp(logc) 
    for i in 1:k
        logc = logc + log(n-i+1) - log(i) + logodds
        s += exp(logc)
    end
    s
end

# %% [markdown]
# ## n = 100

# %%
n, p = 10^2, 0.5
bin = Binomial(n, p)

# %%
@btime cdf($bin, $(n÷2))

# %%
@btime mycdf_naive($bin, $(n÷2))

# %%
@btime mycdf($bin, $(n÷2))

# %%
@btime sum(pdf($bin, i) for i in $(0:n÷2))

# %% [markdown]
# ## n = 1000

# %%
n, p = 10^3, 0.5
bin = Binomial(n, p)

# %%
@btime cdf($bin, $(n÷2))

# %%
@btime mycdf_naive($bin, $(n÷2))

# %%
@btime mycdf($bin, $(n÷2))

# %%
@btime sum(pdf($bin, i) for i in $(0:n÷2))

# %% [markdown]
# ## n = 10000

# %%
n, p = 10^4, 0.5
bin = Binomial(n, p)

# %%
@btime cdf($bin, $(n÷2))

# %%
@btime mycdf_naive($bin, $(n÷2))

# %%
@btime mycdf($bin, $(n÷2))

# %%
@btime sum(pdf($bin, i) for i in $(0:n÷2))

# %% [markdown]
# ## n = 100000

# %%
n, p = 10^5, 0.5
bin = Binomial(n, p)

# %%
@btime cdf($bin, $(n÷2))

# %%
@btime mycdf_naive($bin, $(n÷2))

# %%
@btime mycdf($bin, $(n÷2))

# %%
@btime sum(pdf($bin, i) for i in $(0:n÷2))

# %% [markdown]
# ## n = 1000000

# %%
n, p = 10^6, 0.5
bin = Binomial(n, p)

# %%
@btime cdf($bin, $(n÷2))

# %%
@btime mycdf_naive($bin, $(n÷2))

# %%
@btime mycdf($bin, $(n÷2))

# %%
@btime sum(pdf($bin, i) for i in $(0:n÷2))

# %%
