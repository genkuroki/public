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
#     display_name: Julia 1.10.5
#     language: julia
#     name: julia-1.10
# ---

# %%
module O

using OffsetArrays

mutable struct POA
    n::Int
    a::OffsetVector{Float64}
    A::OffsetArray{Float64, 2}
    x::OffsetVector{Float64}
    POA() = new()
end

function test_oa!(p::POA)
    range = 0:p.n-1
    a, A, x = p.a, p.A, p.x
    for i in range
        for j in range
            x[i] += A[j,i]*a[i]
        end
    end
end

p = POA()
p.n = 1000
p.a = OffsetVector{Float64}(rand(Float64,p.n), 0:p.n-1)
p.A = OffsetArray{Float64}(rand(Float64,p.n,p.n), 0:p.n-1, 0:p.n-1)
p.x = OffsetVector{Float64}(zeros(Float64, p.n), 0:p.n-1)

end

@time O.test_oa!(O.p)
@time O.test_oa!(O.p)

# %%
@code_warntype O.test_oa!(O.p)

# %%
module Q

using OffsetArrays

mutable struct POA{
        Ta<:AbstractVector,
        TA<:AbstractMatrix,
        Tx<:AbstractVector
    }
    a::Ta
    A::TA
    x::Tx
end

function test_oa!(p::POA)
    (; a, A, x) = p
    for i in axes(A, 2)
        for j in axes(A, 1)
            x[i] += A[j,i]*a[i]
        end
    end    
end

m = 10
n = 20
a = OffsetVector{Float64}(collect(1:n), 0:n-1)
A = OffsetMatrix{Float64}((1:m) * (1:n)', 0:m-1, 0:n-1)
x = OffsetVector{Float64}(zeros(n), 0:n-1)
p = POA(a, A, x)

end

@time Q.test_oa!(Q.p)
@time Q.test_oa!(Q.p)

# %%
@code_warntype Q.test_oa!(Q.p)

# %%
q = Q.POA(Q.p.a.parent, Q.p.A.parent, zero(Q.p.x.parent))
@time Q.test_oa!(q)
@time Q.test_oa!(q)
q.x == Q.p.x.parent

# %%
