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
#     display_name: Julia 1.10.0
#     language: julia
#     name: julia-1.10
# ---

# %%
using BenchmarkTools

# %%
L = 3
N = (4, 3, 2)
T = zeros(N...)
A = [collect(1:N[l]) for l in 1:L]

# %%
keys(T)

# %%
collect(keys(T))

# %%
for i in keys(T)
    T[i] = prod(k -> A[k][i[k]], eachindex(A))
end
T

# %%
@btime begin
    for i in keys(TT)
        TT[i] = prod(k -> AA[k][i[k]], eachindex(AA))
    end
    TT
end setup = begin
    TT = copy(T)
    AA = deepcopy(A)
end

# %%
