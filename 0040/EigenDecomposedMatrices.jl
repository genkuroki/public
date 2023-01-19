# -*- coding: utf-8 -*-
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
#     display_name: Julia 1.9.0-beta2
#     language: julia
#     name: julia-1.9
# ---

# %% [markdown]
# ## The Original EigenDecompression.EigenDecompose
#
# See
#
# * https://mobile.twitter.com/realize_ss/status/1615160291108745216
# * https://qiita.com/lelele/items/8408410a94f5c6b8f76e

# %%
using LinearAlgebra
M = rand(100,100)#対角化したい行列
E, P = eigen(M)

# %%
exp(eigen(M))

# %%
module EigenDecompression

export EigenDecompose, eigDecomp
using LinearAlgebra
import Base.*, Base./

#対角化された行列型
struct EigenDecompose{T<:Number} <: AbstractMatrix{T}
    P::AbstractMatrix{T}
    D::Diagonal{T}
    invP::AbstractMatrix{T}
end

#普通のMatrixを対角化する
function eigDecomp(mat::AbstractMatrix)
    E, P = eigen(mat)
    EigenDecompose(P, Diagonal(E), inv(P))
end

#EigenDecompose型に対する関数
Base.exp(eig::EigenDecompose) = EigenDecompose(eig.P, exp(eig.D), eig.invP)
*(eig::EigenDecompose, vec::AbstractVector) = eig.P * eig.D * eig.invP * vec
*(eig::EigenDecompose, sc::Number) = EigenDecompose(eig.P, eig.D*sc, eig.invP)
/(eig::EigenDecompose, sc::Number) = EigenDecompose(eig.P, eig.D/sc, eig.invP)

#普通のMatrixに戻す
Base.Array(eig::EigenDecompose) = eig.P * eig.D * eig.invP

end

# %%
using .EigenDecompression

M = rand(100, 100)
eM = eigDecomp(M)
for i in 1:100
    v = rand(100)
    rnd = rand()
    @assert exp(M*rnd)*v ≈ exp(eM*rnd)*v
end

# %%
using BenchmarkTools

M = rand(100, 100);

#普通な方
function bench1(M)
    for i in 1:100
        v = rand(100)
        exp(M*rand())*v
    end
end

#今回実装した方
function bench2(M)
    eM = eigDecomp(M)
    for i in 1:100
        v = rand(100)
        exp(eM*rand())*v
    end
end

# %%
@benchmark bench1(M)

# %%
@benchmark bench2(M)

# %% [markdown]
# ## EigenDecomposedMatrices.EigenDecomposed

# %%
module EigenDecomposedMatrices

export EigenDecomposed

using LinearAlgebra

struct EigenDecomposed{
        T,
        TA<:AbstractMatrix,
        TE<:AbstractVector,
        TP<:AbstractMatrix,
        TinvP<:AbstractMatrix
    } <: AbstractMatrix{T}
    A::TA
    E::TE
    P::TP
    invP::TinvP
end

function EigenDecomposed(A::AbstractMatrix, E::AbstractVector, P::AbstractMatrix, invP::AbstractMatrix)
    EigenDecomposed{eltype(A), typeof(A), typeof(E), typeof(P), typeof(invP)}(A, E, P, invP)
end

function EigenDecomposed(M::AbstractMatrix)
    E, P = eigen(M)
    A = oftype(P, M)
    invP = ishermitian(M) ? P' : inv(P)
    EigenDecomposed(A, E, P, invP)
end

Base.size(ed::EigenDecomposed) = size(ed.A)
Base.getindex(ed::EigenDecomposed, I...) = getindex(ed.A, I...)
Base.convert(::Type{Array}, ed::EigenDecomposed) = ed.A

function exp_old(ed::EigenDecomposed)
    (; A, E, P, invP) = ed
    expE = exp.(E)
    expA = P * Diagonal(expE) * invP
    EigenDecomposed(expA, expE, P, invP)
end
Base.exp(ed::EigenDecomposed) = exp_eigendecomposed(ed)
Base.:*(c::Number, ed::EigenDecomposed) = EigenDecomposed(c*ed.A, c*ed.E, ed.P, ed.invP)
Base.:*(ed::EigenDecomposed, c::Number) = EigenDecomposed(ed.A*c, ed.E*c, ed.P, ed.invP)
Base.:\(c::Number, ed::EigenDecomposed) = EigenDecomposed(c\ed.A, c\ed.E, ed.P, ed.invP)
Base.:/(ed::EigenDecomposed, c::Number) = EigenDecomposed(ed.A/c, ed.E/c, ed.P, ed.invP)
Base.:*(ed::EigenDecomposed, v::AbstractVector) = ed.A * v

"""
    exp_eigendecomposed!(Y, ed::EigenDecomposed, expE=similar(ed.E), tmpY=similar(Y))

returns the exponential of `ed` and stores the result in `Y`, overwriting the existing value of `Y`. 
It does not overwrite `ed` and uses `expE` and `tmpY` as workspaces.
"""
function exp_eigendecomposed!(Y, ed::EigenDecomposed, expE=similar(ed.E), tmpY=similar(Y))
    (; A, E, P, invP) = ed
    @. expE = exp.(E)
    mul!(tmpY, P, Diagonal(expE))
    mul!(Y, tmpY, invP)
end
exp_eigendecomposed(ed::EigenDecomposed) = exp_eigendecomposed!(similar(ed.P), ed)
LinearAlgebra.lmul!(c::Number, ed::EigenDecomposed) = (lmul!(c, ed.A); lmul!(c, ed.E))
LinearAlgebra.rmul!(ed::EigenDecomposed, c::Number) = (rmul!(ed.A, c); rmul!(ed.E, c))
LinearAlgebra.ldiv!(c::Number, ed::EigenDecomposed) = (ldiv!(c, ed.A); ldiv!(c, ed.E))
LinearAlgebra.rdiv!(ed::EigenDecomposed, c::Number) = (rdiv!(ed.A, c); rdiv!(ed.E, c))

for T in (AbstractVector, AbstractMatrix)
    @eval function LinearAlgebra.mul!(y::$T, ed::EigenDecomposed, x::$T, alpha::Number, beta::Number)
        mul!(y, ed.A, x, alpha, beta)
    end
end

end

# %%
?EigenDecomposedMatrices.exp_eigendecomposed!

# %%
methods(EigenDecomposedMatrices.EigenDecomposed)

# %%
methodswith(EigenDecomposedMatrices.EigenDecomposed)

# %%
methods(EigenDecomposedMatrices.exp_eigendecomposed!)

# %%
n = 2^8
M = 5I + randn(n, n)
v = randn(n)
c = randn()

edM = EigenDecomposedMatrices.EigenDecomposed(M)

Y = similar(edM.A)
expE = similar(edM.E)
tmpY = similar(Y)

y = similar(v)
tmpy = oftype(edM.E, y)
alpha = randn()
beta = randn();

# %%
edM

# %%
dump(edM)

# %%
M ≈ Matrix(edM)

# %%
M ≈ edM

# %%
c*M ≈ c*edM ≈ edM*c

# %%
c\M ≈ c\edM ≈ edM/c

# %%
(
    exp(M)
    ≈ exp(edM)
    ≈ EigenDecomposedMatrices.exp_old(edM)
    ≈ EigenDecomposedMatrices.exp_eigendecomposed!(Y, edM)
    ≈ EigenDecomposedMatrices.exp_eigendecomposed!(Y, edM, expE, tmpY)
)

# %%
@show typeof(y)
mul!(y, M, v, alpha, beta) ≈ mul!(y, edM, v, alpha, beta)

# %% tags=[]
@show typeof(tmpy)
(
    exp(M) * v
    ≈ exp(edM) * v
    ≈ EigenDecomposedMatrices.exp_old(edM) * v
    ≈ mul!(tmpy, EigenDecomposedMatrices.exp_eigendecomposed!(Y, edM), v)
    ≈ mul!(tmpy, EigenDecomposedMatrices.exp_eigendecomposed!(Y, edM, expE, tmpY), v)
)

# %%
@btime edM = EigenDecomposedMatrices.EigenDecomposed(M);

# %%
@btime exp($M) * $v
@btime exp($edM) * $v
@btime EigenDecomposedMatrices.exp_old($edM) * $v
@btime mul!($tmpy, EigenDecomposedMatrices.exp_eigendecomposed!($Y, $edM), $v)
@btime mul!($tmpy, EigenDecomposedMatrices.exp_eigendecomposed!($Y, $edM, $expE, $tmpY), $v);

# %%
n2 = 2^8
M2 = Symmetric(5I + randn(n2, n2))
v2 = randn(n)
c2 = randn()

edM2 = EigenDecomposedMatrices.EigenDecomposed(M2)

Y2 = similar(edM2.A)
expE2 = similar(edM2.E)
tmpY2 = similar(Y2)

y2 = similar(v2)
tmpy2 = oftype(edM2.E, y2)
alpha2 = randn()
beta2 = randn();

# %%
edM2

# %%
dump(edM2)

# %%
M2 ≈ Matrix(edM2)

# %%
M2 ≈ edM2

# %%
c2*M2 ≈ c2*edM2 ≈ edM2*c2

# %%
c2\M2 ≈ c2\edM2 ≈ edM2/c2

# %%
(
    exp(M2)
    ≈ exp(edM2)
    ≈ EigenDecomposedMatrices.exp_old(edM2)
    ≈ EigenDecomposedMatrices.exp_eigendecomposed!(Y2, edM2)
    ≈ EigenDecomposedMatrices.exp_eigendecomposed!(Y2, edM2, expE2, tmpY2)
)

# %%
@show typeof(y2)
mul!(y2, M2, v2, alpha2, beta2) ≈ mul!(y2, edM2, v2, alpha2, beta2)

# %%
@show typeof(tmpy2)
(
    exp(M2) * v2
    ≈ exp(edM2) * v2
    ≈ EigenDecomposedMatrices.exp_old(edM2) * v2
    ≈ mul!(tmpy2, EigenDecomposedMatrices.exp_eigendecomposed!(Y2, edM2), v2)
    ≈ mul!(tmpy2, EigenDecomposedMatrices.exp_eigendecomposed!(Y2, edM2, expE2, tmpY2), v2)
)

# %%
@btime edM2 = EigenDecomposedMatrices.EigenDecomposed(M2);

# %%
@btime exp($M2) * $v2
@btime exp($edM2) * $v2
@btime EigenDecomposedMatrices.exp_old($edM2) * $v2
@btime mul!($tmpy2, EigenDecomposedMatrices.exp_eigendecomposed!($Y2, $edM2), $v2)
@btime mul!($tmpy2, EigenDecomposedMatrices.exp_eigendecomposed!($Y2, $edM2, $expE2, $tmpY2), $v2);

# %%
