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
ENV["COLUMNS"] = 100
using SymPy

a(n) = 8(n-1) + 15
b(n) = 2^n - 1

# %%
@syms n
@show a(4) a(7)
@show b(1) (b(n+1) - (2b(n)+1)).simplify()
;

# %%
as = [a(n) for n in 1:45]
as'

# %%
bs = [b(k) for k in 1:20]
bs'

# %%
cs = setdiff(as, bs)
cs'

# %%
sum(cs)

# %%
