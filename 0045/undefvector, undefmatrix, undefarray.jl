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
#     display_name: Julia 1.9.3
#     language: julia
#     name: julia-1.9
# ---

# %%
undefvector(T::Type, m) = Vector{T}(undef, m)
undefvector(m) = undefvector(Float64, m)
undefmatrix(T::Type, m, n) = Matrix{T}(undef, m, n)
undefmatrix(m, n) = undefmatrix(Float64, m, n)
undefarray(T::Type, n...) = Array{T}(undef, n...)
undefarray(n...) = undefarray(Float64, n...)

# %%
undefvector(4)

# %%
undefmatrix(4, 5)

# %%
undefmatrix(Int, 4, 5)

# %%
undefarray(4, 3, 2)

# %%
?undefarray

# %%
