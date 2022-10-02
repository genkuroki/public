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
#     display_name: Julia 1.8.2
#     language: julia
#     name: julia-1.8
# ---

# %%
using Distributions

# %%
hg = Hypergeometric(4, 136-4, 13)
[(k, pdf(hg, k)) for k in 0:4]

# %%
hg = Hypergeometric(4, 136-4, 13)
[(k, ccdf(hg, k-1)) for k in 0:4]

# %%
hg = Hypergeometric(4, 136-4, 14)
[(k, pdf(hg, k)) for k in 0:4]

# %%
hg = Hypergeometric(4, 136-4, 14)
[(k, ccdf(hg, k-1)) for k in 0:4]

# %%
hg = Hypergeometric(7, 136-7, 13)
[(k, pdf(hg, k)) for k in 0:7]

# %%
hg = Hypergeometric(7, 136-7, 13)
[(k, ccdf(hg, k-1)) for k in 0:7]

# %%
hg = Hypergeometric(7, 136-7, 14)
[(k, pdf(hg, k)) for k in 0:7]

# %%
hg = Hypergeometric(7, 136-7, 14)
[(k, ccdf(hg, k-1)) for k in 0:7]

# %%
hg = Hypergeometric(7, 136-7, 19)
[(k, pdf(hg, k)) for k in 0:7]

# %%
hg = Hypergeometric(7, 136-7, 19)
[(k, ccdf(hg, k-1)) for k in 0:7]

# %%
hg = Hypergeometric(7, 136-7, 18)
[(k, ccdf(hg, k-1)) for k in 0:7]

# %%
