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
#     display_name: Julia 1.9.0-beta3
#     language: julia
#     name: julia-1.9
# ---

# %% [markdown]
# * https://twitter.com/YutaTasaki/status/1617668862165778432
# * https://gist.github.com/tagdtm/3d923b706430e80973463d76a0b3605e

# %%
nx = 2
ny = 2
a = zeros(2,2)
b = zeros(2,2)
c = zeros(2,2)

a = [1 2; 3 4]
b = [1.1 1.8; 1.4 1.7]

for i in 1:nx
    for j in 1:ny
      c[i,j] = a[i,j] * a[floor.(Int, b[i,j]), ceil.(Int, b[i,j])]
    end
end

#=
I would like to find simpler way of this without for loop.
The expected form is something like
c = a .* a[floor.(Int64,b),ceil.(Int64,b)]
But it generates multi dimension array, which in this case 4d array.
=#

# %%
c

# %%
c1 = zeros(2,2)
@. c1 = a * getindex((a,), floor(Int, b), ceil(Int, b))

# %%
c1 == c

# %%
m, n = 5, 5
A = round.(10rand(m, n))
B = round.(1 .+ 3rand(m, n); digits=1)
C = similar(A)

for k in keys(A)
    C[k] = A[k] * A[floor.(Int, B[k]), ceil.(Int, B[k])]
end
C

# %%
C1 = similar(A)
@. C1 = A * getindex((A,), floor(Int, B), ceil(Int, B))

# %%
C1 == C

# %%
