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
using SymPy

M = Sym[
     6  5 -2 -3
    -3 -1  3  3
     2  1 -2 -3
    -1  1  5  5
]
P, J = M.jordan_form()
display(P)
display(J)
P*J/P

# %%
Q = Rational.(Int.(P))
K = Int.(J)
display(Q)
display(K)
Q*K/Q

# %%
