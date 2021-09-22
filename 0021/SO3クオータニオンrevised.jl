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
using SymPy, LinearAlgebra

# %%
struct Quaternion{T}
  w::T
  x::T
  y::T
  z::T
end

# %%
Base.:*(p::Quaternion, q::Quaternion) = (
  Quaternion(p.w*q.w -p.x*q.x -p.y*q.y -p.z*q.z
    ,p.w*q.x + p.x*q.w + p.y*q.z -p.z*q.y 
    ,p.w*q.y + p.y*q.w + p.z*q.x -p.x*q.z 
    ,p.w*q.z + p.z*q.w + p.x*q.y -p.y*q.x 
    )
)

Base.:*(α, q::Quaternion) =  Quaternion( α*q.w, α*q.x, α*q.y, α*q.z)
Base.:+(p::Quaternion, q::Quaternion) = Quaternion(p.w+q.w, p.x+q.x, p.y+q.y, p.z+q.z)
Base.:-(p::Quaternion, q::Quaternion) = Quaternion(p.w-q.w, p.x-q.x, p.y-q.y, p.z-q.z)
Base.:(==)(p::Quaternion, q::Quaternion) = (p.w, p.x, p.y, p.z) == (q.w, q.x, q.y, q.z)

Base.show(io::IO, q::Quaternion) = show(io::IO, [q.w, q.x, q.y, q.z])
Base.show(io::IO, ::MIME"text/plain", q::Quaternion) = show(io::IO, MIME("text/plain"), [q.w, q.x, q.y, q.z])
Base.show(io::IO, ::MIME"text/latex", q::Quaternion{Sym}) = show(io, MIME("text/latex"), [q.w, q.x, q.y, q.z])

SymPy.expand(q::Quaternion) = Quaternion(expand(q.w), expand(q.x), expand(q.y), expand(q.z))

# %%
@syms w1::real w2::real x1::real x2::real y1::real y2::real z1::real z2::real 
@syms w3::real x3::real y3::real z3::real
@syms w::real x::real y::real z::real
@syms a::real b::real c::real
@syms a1::real b1::real c1::real
@syms a2::real b2::real c2::real
@syms α::real, β::real
@syms θ::real

# %%
q = Quaternion(w,x,y,z)
q1 = Quaternion(w1,x1,y1,z1)
q2 = Quaternion(w2,x2,y2,z2)
q1 * q2

# %% [markdown]
# ## 線形性の確認

# %%
Q1 = expand(q*(α*q1 + β*q2))

# %%
Q2 = expand(α*(q*q1) + β*(q*q2))

# %%
Q1 == Q2

# %%
Base.adjoint(q::Quaternion) = Quaternion(q.w, -q.x, -q.y , -q.z)
Base.abs2(p::Quaternion) = p.w*p.w + p.x*p.x + p.y*p.y + p.z*p.z

# %%
q1*adjoint(q1)

# %%
abs2(q1)

# %%
abs2(q1*q2)

# %%
expand(subs(expand(abs2(q1*q2)), w1^2 => 1-x1^2-y1^2-z1^2))

# %%
subs(abs2(q1*q2), w1^2 => 1-x1^2-y1^2-z1^2)

# %%
expand(subs(expand(abs2(q1*q2)), w1^2 +x1^2+y1^2+z1^2=> 1))

# %%
subs(expand(abs2(q1*q2)), w1^2 => 1-x1^2-y1^2-z1^2)

# %%
v = Quaternion(zero(Sym),a,b,c)

# %%
q = Quaternion(w,x,y,z)

# %%
expand(q*(v)*adjoint(q))

# %%
(q*(v)*adjoint(q))

# %%
expand(q*Quaternion(0,1,0,0)*adjoint(q))
expand(q*Quaternion(0,0,1,0)*adjoint(q))
expand(q*Quaternion(0,0,0,1)*adjoint(q))

# %%
Base.collect(q::Quaternion) = [q.w, q.x, q.y, q.z]

# %%
M = [collect(expand(q*Quaternion(0,1,0,0)*adjoint(q))) collect(expand(q*Quaternion(0,0,1,0)*adjoint(q))) collect(expand(q*Quaternion(0,0,0,1)*adjoint(q)))]

# %%
MtM = expand.(M' * M)

# %%
expand.(subs.(MtM, w^2 => 1-x^2-y^2-z^2))

# %%
detM = det(M[2:end,:])

# %%
expand.(subs(detM, w^2 => 1-x^2-y^2-z^2))

# %%
M2 = subs.(subs.(M,x,0),y,0)

# %%
M3 = subs.(M2, w, cos(θ))

# %%
M4 = subs.(M3, z, sin(θ))

# %%
simplify.(M4)

# %%
M3 = subs.(M2, w, cos(θ/2))

# %%
M4 = subs.(M3, z, sin(θ/2))

# %%
simplify.(M4)

# %%
