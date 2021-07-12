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
Base.:(==)(x::F{p}, y::F{p}) where p = x.a == y.a
Base.:<(x::F{p}, y::F{p}) where p = x.a < y.a

Base.show(io::IO, x::F{p}) where p = print(io, "F", p, '(', x.a, ')')

F7 = F{7, BigInt}
x, y = F7(10), F7(-2)
@show(x, y, zero(x), one(x), +x, -x, x + y, x - y, x * y, x / y, x \ y, x^3, x^-5, x == F7(3), x == 3)
println()
A = F7[1 2; 3 4]
@show(A, 4A, A/4)
println()
using LinearAlgebra
L, U = lu(A)
@show(L, U, L * U, det(A), inv(A), inv(A) * A, A * inv(A));

# %%
"""See https://docs.julialang.org/en/v1/manual/interfaces/"""
Base.iterate(Fp::Type{F{p, T}}) where {p, T<:Integer} = (zero(Fp), zero(T))
function Base.iterate(Fp::Type{F{p, T}}, state) where {p, T<:Integer}
    nextstate = state + 1
    nextstate < p ? (Fp(nextstate), nextstate) : nothing
end
Base.IteratorSize(Fp::Type{F{p, T}}) where {p, T<:Integer} = Base.HasLength()
Base.length(Fp::Type{F{p, T}}) where {p, T<:Integer} = p
Base.eltype(Fp::Type{F{p, T}}) where {p, T<:Integer} = Fp

squares(Fp) = Fp[x^2 for x in Fp]
squareroots(k, Fp) = Fp[x for x in Fp if x^2 == k]
@show(collect(F7), squares(F7), squareroots.(0:6, Ref(F7)));

# %%
using GaloisFields, LinearAlgebra
GF7 = @GaloisField 7
B = GF7[1 2; 3 4]
det(B)

# %%
inv(B)

# %%
lu(B)

# %%
using AbstractAlgebra

@show GF7 = GF(7)
@show P, x = PolynomialRing(GF7, "x")
@show C = GF7[1 2; 3 4]
squares(Fp) = [x^2 for x in Fp]
squareroots(k, Fp) = [x for x in Fp if x^2 == k]
println()
@show(det(C), inv(C), lu(C), charpoly(P, C))
println()
@show(collect(GF7), squares(GF7), squareroots.(0:6, Ref(GF7)));

# %%
@which iterate(GF7)

# %%
@code_warntype iterate(F7, big(1))

# %%
[x for x in F7 if x^2 == 4]

# %%
[x for x in F7 if x^2 == 5]

# %%
@which collect(x for x in F7 if x^2 == 5)

# %%
Base.@default_eltype(x for x in F7 if x^2 == 5)

# %%
Base.@default_eltype(x for x in F7 if x^2 == 5)

# %%
(@macroexpand Base.@default_eltype(x for x in F7 if x^2 == 5)) |> Base.remove_linenums!

# %%
iter = (x for x in F7 if x^2 == 5)

# %%
Base._iterator_upper_bound(iter)

# %%
Core.Compiler.return_type(Base._iterator_upper_bound, Tuple{typeof(iter)})

# %%
@which Base.inferencebarrier(nothing)

# %%
iter2 = (x for x in F7 if x^2 == 4)

# %%
Base._iterator_upper_bound(iter2)

# %%
iter3 = (x for x in GF7 if x^2 == 4)
collect(iter3)

# %%
iter4 = (x for x in GF7 if x^2 == 5)
collect(iter4)

# %%
Core.Compiler.return_type(Base._iterator_upper_bound, Tuple{typeof(iter4)})

# %%
Base._iterator_upper_bound(iter4)

# %%
@code_warntype Base._iterator_upper_bound(iter4)

# %%
@code_warntype Base._iterator_upper_bound(iter2)

# %%
@code_warntype iterate(iter4)

# %%
@code_warntype iterate(iter)

# %%
Core._apply_iterate(Base.iterate, Base.iterate, Core.tuple(iter3.iter), ())

# %%
Core._apply_iterate(Base.iterate, Base.iterate, Core.tuple(iter2.iter), ())

# %%
Core.tuple(iter.iter)

# %%
Core.tuple(iter4.iter)

# %%
