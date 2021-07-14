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
#     display_name: Julia 1.6.1
#     language: julia
#     name: julia-1.6
# ---

# %%
using SparseArrays

A = zeros( Integer, 2400, 2400 )
for i in 1:2400
    if rand(1:4) == 1
        for j in i+1:2400
            if rand(1:15) == 1
                A[i,j] = 1
                A[j,i] = 1
            end
        end
    end
end

sum( A)
A = sparse( A )
B = A^2

println( "some stats A: size ", size(A), " min value ", minimum( A ), " max value ", maximum( A ) )
println( "some stats B: size ", size(B), " min value ", minimum( B ), " max value ", maximum( B ) )

# %%
A

# %%
B

# %%
a = spzeros(Integer, 5, 5)
a .= 1
a

# %%
b = a^2

# %%
