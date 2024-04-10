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
#     display_name: Julia 1.10.2
#     language: julia
#     name: julia-1.10
# ---

# %%
a = b = 1e-200
c = d = 1e200

@show a * b * c * d
@show a * (b * c) * d
@show a * b * (c * d);

# %%
@show exp(log(a) + log(b) + log(c) + log(d))
@show exp(log(a) + (log(b) + log(c)) + log(d))
@show exp(log(a) + log(b) + (log(c) + log(d)));

# %%
