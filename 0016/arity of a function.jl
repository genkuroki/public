# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:hydrogen
#     text_representation:
#       extension: .jl
#       format_name: hydrogen
#       format_version: '1.3'
#       jupytext_version: 1.11.2
#   kernelspec:
#     display_name: Julia 1.6.2
#     language: julia
#     name: julia-1.6
# ---

# %%
arity(f) = extrema(m -> m.nargs, methods(f)) .- 1
findarity(f, k) = (ms = collect(methods(f)); ms[findall(m -> m.nargs - 1 == k, ms)])

# %%
arity(+)

# %%
findarity(+, 5)

# %%
findarity(+, 4)

# %%
@which +(1, 2)

# %%
@which +(1, 2, 3)

# %%
@which +(1, 2, 3, 4)

# %%
]add CairoMakie@0.6.3

# %%
