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
struct F{p, T<:Integer} <: Integer
    a::T
    F{p, T}(a::Integer) where {p, T<:Integer} = new{p, T}(mod(T(a), p))
end
F{p}(a::Integer) where p = F{p, typeof(a)}(mod(a, p))
F{p, S}(x::F{p, T}) where {S<:Integer, p, T<:Integer} = F{p, S}(x.a)

Base.promote_rule(::Type{F{p, T}}, ::Type{S}) where {p, T<:Integer, S<:Integer} =
    F{p, promote_type(T, S)}

Base.zero(::Type{F{p, T}}) where {p, T<:Integer} = F{p}(mod(zero(T), p))
Base.one(::Type{F{p, T}}) where {p, T<:Integer} = F{p}(mod(one(T), p))
for op in (:-, :+)
    @eval Base.$op(x::F{p}) where p = F{p}(mod($op(x.a), p))
end
for op in (:-, :+, :*)
    @eval Base.$op(x::F{p}, y::F{p}) where p = F{p}(mod($op(x.a, y.a), p))
end
Base.inv(x::F{p}) where p = F{p}(invmod(x.a, p))
Base.:/(x::F{p}, y::F{p}) where p = x * inv(y)
Base.:\(x::F{p}, y::F{p}) where p = inv(x) * y
Base.:<(x::F{p}, y::F{p}) where p = x.a < y.a

Base.show(io::IO, x::F{p}) where p = print(io, "F", p, '(', x.a, ')')

F7 = F{7, BigInt}
x, y = F7(10), F7(-2)
@show(x, y, zero(x), one(x), +x, -x, x + y, x - y, x * y, x / y, x \ y)
println()
A = F7[1 2; 3 4]
@show(A, 4A, A/4)
println()
using LinearAlgebra
L, U = lu(A, NoPivot())
@show(L * U, det(A), inv(A), inv(A) * A, A * inv(A));

# %%
