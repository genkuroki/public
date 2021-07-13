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
#     display_name: Julia 1.8.0-DEV
#     language: julia
#     name: julia-1.8
# ---

# %%
using GaloisFields, LinearAlgebra
const F5 = @GaloisField 5
A = F5[1 2; 0 3]

# %%
det(A)

# %%
B = F5[1 2; 3 4]

# %%
det(B)

# %%
F5(3) |> dump

# %%
F5 |> dump

# %%
Base.abs(x::GaloisFields.PrimeField) = x.n # type-piracy!
det(B)

# %%
F3 = @GaloisField 3
[det(F3.(rand(0:2, 3, 3))) for _ in 1:20]

# %%
Base.abs(x::GaloisFields.ExtensionField) = x.coeffs # type-piracy!
Base.isless(x::GaloisFields.PrimeField, y::GaloisFields.PrimeField) = isless(x.n, y.n) # type-piracy!
F9 = @GaloisField! 3 a^2 + 1
C = F9[1+a 2+a; 2+a 2+a]

# %%
det(C)

# %%
