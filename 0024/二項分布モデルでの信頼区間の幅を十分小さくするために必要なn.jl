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
#     display_name: Julia 1.6.4
#     language: julia
#     name: julia-1.6
# ---

# %% [markdown]
# https://twitter.com/snap_tck/status/1463455803709423617 の確認

# %%
using Distributions
using Roots
using Plots
using StatsFuns

# %%
?logistic

# %%
ccdf(Binomial(1608, 0.001), 1 - 1)

# %%
ccdf(Binomial(1609, 0.001), 1 - 1)

# %% tags=[]
ccdf(Binomial(2000, 0.001), 1 - 1)

# %%
n = 16000
p = 0.001
bin = Binomial(n, p)
cdf(bin, 0.5n*p), ccdf(bin, 1.5n*p - 1)

# %%
normal = Normal(mean(bin), std(bin))
cdf(normal, n*(p - 0.0005)), ccdf(normal, n*(p + 0.0005))

# %%
2*√(0.001(1 - 0.001))/√1600

# %%
2*√(0.001(1 - 0.001))/√4000

# %%
2*√(0.001(1 - 0.001))/√16000

# %%
x ⪅ y = x < y || x ≈ y

function pval(n, p, k)
    bin = Binomial(n, p)
    p0 = pdf(bin, k)
    min(1, sum(pdf(bin, j) for j in support(bin) if pdf(bin, j) ⪅ p0))
end

function ci(n, k, α = 0.05)
    CI = logistic.(find_zeros(t -> pval(n, logistic(t), k) - α, -50, 50))
end

# %%
CI = ci(16000, 16)

# %%
plot(p -> pval(16000, p, 16), 0, 0.003; label="p-value function for n = 16000, k = 16")
plot!(; xtick=0:0.0005:0.003)

# %%
