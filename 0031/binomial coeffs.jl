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
binomial(big(2000), 300)

# %%
using SpecialFunctions
logbinomial(n, k) = logfactorial(n) - logfactorial(k) - logfactorial(n-k)
a = logbinomial(2000, 300)

# %%
exp(a)

# %%
logprob(n, p, k) = logbinomial(n, k) + k*log(p) + (n-k)*log(1-p)
b = logprob(2000, 0.3, 200)

# %%
exp(b)

# %%
using StatsFuns
@show c = logsumexp(logprob(2000, 0.3, k) for k in 0:560)
exp(c)

# %%
using Distributions
cdf(Binomial(2000, 0.3), 560)

# %%
