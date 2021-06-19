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
#     display_name: Julia 1.8.0-DEV
#     language: julia
#     name: julia-1.8
# ---

# %%
# privateとpublicの区別は`Base.propertynames`の定義で行える。
@doc propertynames

# %%
module O

"""
抽象平面ベクトルの型
"""
abstract type AbstractPlanarVector{T} end

"""タプル版和集合(重複を削除)"""
tupleunion(x) = Tuple(unique(x))
tupleunion(x, y) = Tuple(unique((x..., y...)))
tupleunion(x, y, z...) = tupleunion(tupleunion(x, y), z...)

"""
抽象平面ベクトルのpublic propertiesは x と y

タブ補完でもx, yしか表示されなくなる。
"""
function Base.propertynames(p::AbstractPlanarVector, private=false)
    public_properties = (:x, :y)
    private && return tupleunion(public_properties, fieldnames(p))
    public_properties
end

"""
抽象平面ベクトルは (x, y) の形式で表示される
"""
Base.show(io::IO, p::AbstractPlanarVector) = show(io, (p.x, p.y))

"""
AbstractPlanarVector{T}型の抽象平面ベクトルの成分の型は T
"""
Base.eltype(x::AbstractPlanarVector{T}) where T = T

"""
平面ベクトルの型

x, yをfieldsに持ち、それらがそのままpublic propertiesになる。
"""
struct PlanarVector{T} <: AbstractPlanarVector{T} x::T; y::T end

"""
標準基底の型

x, yをfieldsとして持たず、型のパラメータ i で何番目の基底ベクトルであるかを指定。
"""
struct CanonBasis{T, i} <: AbstractPlanarVector{T} end
function CanonBasis{T}(i) where T
    @assert i == 1 || i == 2
    CanonBasis{T, i}()
end
CanonBasis(i) = CanonBasis{Int}(i)

"""
標準基底のpublic propertiesの取得メソッド

x, yをfieldsとして持たないので定義してやる必要がある。

* i == 1 ⇒ (x, y) = (1, 0)
* i == 2 ⇒ (x, y) = (0, 1)
"""
function Base.getproperty(::CanonBasis{T, i}, name::Symbol) where {T, i}
    name === :x && return ifelse(i == 1, one(T), zero(T))
    name === :y && return ifelse(i == 2, one(T), zero(T))
    error("type CanonBasis has no public property $name")
end

########## ベクトルの演算達

# 抽象平面ベクトル p, q について、 +p, -p, p + q, p - q を定義
for op in (:+, :-)
    @eval Base.$op(p::AbstractPlanarVector) = PlanarVector($op(p.x), $op(p.y))
    @eval Base.$op(p::AbstractPlanarVector, q::AbstractPlanarVector) =
        PlanarVector($op(p.x, q.x), $op(p.y, q.y))
end

# 抽象平面ベクトル p とスカラー a について、 a * p, p * a, a \ p, p / a を定義
Base.:*(a, p::AbstractPlanarVector) = PlanarVector(a * p.x, a * p.y)
Base.:*(p::AbstractPlanarVector, a) = PlanarVector(p.x * a, p.y * a)
Base.:\(a, p::AbstractPlanarVector) = PlanarVector(a \ p.x, a \ p.y)
Base.:/(p::AbstractPlanarVector, a) = PlanarVector(p.x / a, p.y / a)

# LinearAlgebra の dot の定義の準備
using LinearAlgebra

##### 抽象ベクトルの内積を定義

# まず、一般的な定義式で内積を定義
LinearAlgebra.dot(p::AbstractPlanarVector, q::AbstractPlanarVector) =
    conj(p.x) * q.x + conj(p.y) * q.y

# 次に、標準基底の場合に特殊化した定義を行う。多重ディスパッチを本質的に使っている。
LinearAlgebra.dot(p::CanonBasis{T, i}, q::AbstractPlanarVector{U}) where {T, i, U} =
    (P = promote_type(T, U); P(i == 1 ? q.x : q.y))
LinearAlgebra.dot(p::AbstractPlanarVector{U}, q::CanonBasis{T, i}) where {U, T, i} =
    (P = promote_type(T, U); P(conj(i == 1 ? p.x : p.y)))

# Juliaの多重ディスパッチでは上の2つの場合のintersectionの場合の定義もしておく必要がある。
LinearAlgebra.dot(p::CanonBasis{T, i}, q::CanonBasis{U, j}) where {T, i, U, j} =
    (P = promote_type(T, U); i == j ? one(P) : zero(P))

end

# %%
p = O.PlanarVector(2+im, 3+im)

# %%
q = O.PlanarVector(-5.0, 10.0)

# %%
p + q

# %%
eltype(p + q)

# %%
-3p + 4q

# %%
O.dot(p, q)

# %%
e1 = O.CanonBasis{Float64}(1)

# %%
e2 = O.CanonBasis(2)

# %%
O.dot(e1, e2)

# %%
O.dot(O.CanonBasis{ComplexF64}(1), e1)

# %%
e1 - e2

# %%
3*e1 - e2/2

# %%
p.x, p.y

# %%
e1.x, e1.y

# %%
getx(p::O.AbstractPlanarVector) = p.x
gety(p::O.AbstractPlanarVector) = p.y

# %%
getx(p), gety(p)

# %%
@code_typed getx(p)

# %%
@code_typed gety(q)

# %%
@code_typed getx(e2)

# %%
@code_typed gety(e2)

# %%
@code_typed O.dot(p, q)

# %%
@code_typed O.dot(e1, q)

# %%
@code_typed O.dot(p, e2)

# %%
@code_llvm debuginfo=:none O.dot(p, e2)

# %%
@code_typed O.dot(e2, e2)

# %%
@code_typed O.dot(e1, e2)

# %%
@code_llvm debuginfo=:none O.dot(e1, e2)

# %%
using BenchmarkTools
n = 10^6
v = [O.PlanarVector(randn(2)...) for _ in 1:n]
a = O.PlanarVector(1, 0)
e = O.CanonBasis(1)

@show sum(Base.Fix1(O.dot, a), v) == sum(Base.Fix1(O.dot, e), v)
@btime sum($(Base.Fix1(O.dot, a)), $v)
@btime sum($(Base.Fix1(O.dot, e)), $v)

# %%
@code_native debuginfo=:none O.dot(a, q)

# %%
@code_native debuginfo=:none O.dot(e, q)

# %%
