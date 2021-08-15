# -*- coding: utf-8 -*-
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
using SymPy

function diffeq_3()
    @vars z t
    P, ∇Pk, Uk, dUdtk = symbols("P ∇Pk Uk dUdtk", cls = sympy.Function)
    P = cos(t)*cos(z)
    ∇Pk = diff.(P, z)
    dUdtk = diff.(Uk(z, t), t)
    diffeq = Eq(dUdtk, ∇Pk)
    display(diffeq)
    solution = pdsolve(diffeq)
end

diffeq_3()

# %%
using SymPy

function diffeq_3()
    @vars z t
    P, ∇Pk, Uk, dUdtk = symbols("P ∇Pk Uk dUdtk", cls = sympy.Function)
    P = cos(t)*cos(z)
    ∇Pk = diff.(P, z)
    dUdtk = diff.(Uk(z, t), t)
    diffeq = Eq(dUdtk, ∇Pk)
    solution = pdsolve(diffeq)
    [diffeq, solution]
end

diffeq_3()

# %%
@show diffeq_3();

# %%
