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
using LinearAlgebra
using BenchmarkTools

# %%
module EigenDecomposedMatrices

export EigenDecomposed

using LinearAlgebra
using Memoization

struct EigenDecomposed{
        T,
        TE<:AbstractVector{T},
        TP<:AbstractMatrix{T},
        TinvP<:AbstractMatrix{T}
    } <: AbstractMatrix{T}
    E::TE
    P::TP
    invP::TinvP
end

function EigenDecomposed(A::AbstractMatrix)
    E, P = eigen(A)
    invP = ishermitian(A) ? P' : inv(P)
    EigenDecomposed(E, P, invP)
end

LinearAlgebra.eigvals(ed::EigenDecomposed) = ed.E
LinearAlgebra.eigvecs(ed::EigenDecomposed) = ed.P
inveigvecs(ed::EigenDecomposed) = ed.invP
@memoize Base.parent(ed::EigenDecomposed) = eigvecs(ed) * Diagonal(eigvals(ed)) * inveigvecs(ed)
Base.convert(::Type{Array}, ed::EigenDecomposed) = convert(Array, parent(ed))
for op in (:eltype, :size)
    @eval Base.$op(ed::EigenDecomposed) = $op(eigvecs(ed))
end
Base.getindex(ed::EigenDecomposed, I...) = getindex(parent(ed), I...)

Base.:*(c::Number, ed::EigenDecomposed) = EigenDecomposed(c*eigvals(ed), eigvecs(ed), inveigvecs(ed))
Base.:*(ed::EigenDecomposed, c::Number) = EigenDecomposed(eigvals(ed)*c, eigvecs(ed), inveigvecs(ed))
Base.:\(c::Number, ed::EigenDecomposed) = EigenDecomposed(c\eigvals(ed), eigvecs(ed), inveigvecs(ed))
Base.:/(ed::EigenDecomposed, c::Number) = EigenDecomposed(eigvals(ed)/c, eigvecs(ed), inveigvecs(ed))
for T in (AbstractVector, AbstractMatrix)
    @eval function Base.:*(ed::EigenDecomposed, v::$T)
        E, P, invP = eigvals(ed), eigvecs(ed), inveigvecs(ed)
        P * (Diagonal(E) * (invP * v))
    end
end

function exp_old(ed::EigenDecomposed)
    E, P, invP = eigvals(ed), eigvecs(ed), inveigvecs(ed)
    expE = exp.(E)
    expA = P * Diagonal(expE) * invP 
    EigenDecomposed(expE, P, invP)
end

LinearAlgebra.lmul!(c::Number, ed::EigenDecomposed) = lmul!(c, eigvals(ed))
LinearAlgebra.rmul!(ed::EigenDecomposed, c::Number) = rmul!(eigvals(ed), c)
LinearAlgebra.ldiv!(c::Number, ed::EigenDecomposed) = ldiv!(c, eigvals(ed))
LinearAlgebra.rdiv!(ed::EigenDecomposed, c::Number) = rdiv!(eigvals(ed), c)

for op in (:exp, :log, :sin, :cos)
    opE = Symbol(op, "E")
    op_eigendecomposed = Symbol(op, "_eigendecomposed")
    op_eigendecomposed! = Symbol(op_eigendecomposed, "!")
    op_eigendecomposed!_doc =
        """
        $op_eigendecomposed!(Y, ed::EigenDecomposed, $opE=similar(ed.E), tmpY=similar(Y))

        returns the `$op` of `ed` and stores the result in `Y`, overwriting the existing value of `Y`. 
        It does not overwrite `ed` and uses `$opE` and `tmpY` as workspaces.
        """
    @eval begin
        @doc $op_eigendecomposed!_doc
        function $op_eigendecomposed!(Y, ed::EigenDecomposed, $opE=similar(ed.E), tmpY=similar(Y))
            E, P, invP = eigvals(ed), eigvecs(ed), inveigvecs(ed)
            @. $opE = $op.(E)
            mul!(tmpY, P, Diagonal($opE))
            mul!(Y, tmpY, invP)
        end
        $op_eigendecomposed(ed::EigenDecomposed) = $op_eigendecomposed!(similar(eigvecs(ed)), ed)
        Base.$op(ed::EigenDecomposed) = $op_eigendecomposed(ed)
    end
end

end

# %%
?EigenDecomposedMatrices.exp_eigendecomposed!

# %%
?EigenDecomposedMatrices.log_eigendecomposed!

# %%
methods(EigenDecomposedMatrices.EigenDecomposed)

# %%
methods(EigenDecomposedMatrices.EigenDecomposed{Float64, Vector{Float64}, Matrix{Float64}, Matrix{Float64}})

# %%
methodswith(EigenDecomposedMatrices.EigenDecomposed)

# %%
methods(EigenDecomposedMatrices.exp_eigendecomposed!)

# %%
methods(EigenDecomposedMatrices.log_eigendecomposed!)

# %%
A = [
    2 -1 0
    -1 2 -1
    0 -1 2
]

edA = EigenDecomposedMatrices.EigenDecomposed(A)

# %%
log(edA)

# %%
log(A)

# %%
log(edA) ≈ log(A)

# %%
n = 2^8
M = 5I + randn(n, n)
v = randn(n)
c = randn()

edM = EigenDecomposedMatrices.EigenDecomposed(M)

Y = similar(eigvecs(edM))
expE = similar(eigvals(edM))
tmpY = similar(Y)

y = similar(eigvals(edM))
alpha = randn()
beta = randn();

# %%
edM

# %%
dump(edM)

# %%
M ≈ parent(edM) == Matrix(edM)

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

# %% tags=[]
@show typeof(y)
(
    exp(M) * v
    ≈ exp(edM) * v
    ≈ EigenDecomposedMatrices.exp_old(edM) * v
    ≈ mul!(y, EigenDecomposedMatrices.exp_eigendecomposed!(Y, edM), v)
    ≈ mul!(y, EigenDecomposedMatrices.exp_eigendecomposed!(Y, edM, expE, tmpY), v)
)

# %%
@btime edM = EigenDecomposedMatrices.EigenDecomposed(M);

# %%
@btime exp($M) * $v
@btime exp($edM) * $v
@btime EigenDecomposedMatrices.exp_old($edM) * $v
@btime mul!($y, EigenDecomposedMatrices.exp_eigendecomposed!($Y, $edM), $v)
@btime mul!($y, EigenDecomposedMatrices.exp_eigendecomposed!($Y, $edM, $expE, $tmpY), $v);

# %%
n2 = 2^8
M2 = Symmetric(5I + randn(n2, n2))
v2 = randn(n)
c2 = randn()

edM2 = EigenDecomposedMatrices.EigenDecomposed(M2)

Y2 = similar(eigvecs(edM2))
expE2 = similar(eigvals(edM2))
tmpY2 = similar(Y2)

y2 = similar(eigvals(edM2))
alpha2 = randn()
beta2 = randn();

# %%
edM2

# %%
dump(edM2)

# %%
M2 ≈ parent(edM2) == Matrix(edM2)

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
(
    exp(M2) * v2
    ≈ exp(edM2) * v2
    ≈ EigenDecomposedMatrices.exp_old(edM2) * v2
    ≈ mul!(y2, EigenDecomposedMatrices.exp_eigendecomposed!(Y2, edM2), v2)
    ≈ mul!(y2, EigenDecomposedMatrices.exp_eigendecomposed!(Y2, edM2, expE2, tmpY2), v2)
)

# %%
@btime edM2 = EigenDecomposedMatrices.EigenDecomposed(M2);

# %%
@btime exp($M2) * $v2
@btime exp($edM2) * $v2
@btime EigenDecomposedMatrices.exp_old($edM2) * $v2
@btime mul!($y2, EigenDecomposedMatrices.exp_eigendecomposed!($Y2, $edM2), $v2)
@btime mul!($y2, EigenDecomposedMatrices.exp_eigendecomposed!($Y2, $edM2, $expE2, $tmpY2), $v2);

# %%
