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
#     display_name: Julia 1.8.0
#     language: julia
#     name: julia-1.8
# ---

# %%
using LinearAlgebra
using Random
Random.seed!(4648373)
using TensorOperations

# %%
A = [
    1 3
    3 4
]

B = [
    5 6
    7 8
]

c = @tensor A[i, j] * B[j, i]

# %%
tr(A*B)

# %%
i = rand(5:13, 4)
j = rand(3:6, 4)
i[1] = 100
i, j

# %%
G = [randn(i[mod1(k, 4)], j[k], i[mod1(k+1, 4)]) for k in 1:4]
size.(G)

# %%
@macroexpand @tensor A[j1, j2, j3, j4] := G[1][i1, j1, i2] * G[2][i2, j2, i3] * G[3][i3, j3, i4] * G[4][i4, j4, i1]

# %%
@tensor A[j1, j2, j3, j4] := G[1][i1, j1, i2] * G[2][i2, j2, i3] * G[3][i3, j3, i4] * G[4][i4, j4, i1]
A

# %%
function refA_Gs(A, G, N)
    Gs = :($G[1][i1, j1, i2])
    for k in 2:N
        ik = Symbol(:i, k)
        jk = Symbol(:j, k)
        ikp1 = Symbol(:i, mod1(k+1, N))
        Gs = :($Gs * $G[$k][$ik, $jk, $ikp1])
    end
    refA = Expr(:ref, :A)
    for k in 1:N
        push!(refA.args, Symbol(:j, k))
    end
    refA, Gs
end

# %%
refA_Gs(:A, :G, 4)

# %%
macro multitrace(G, N)
    refA, Gs = refA_Gs(:A, :G, N)
    quote
        let G = $(esc(G))
            @tensor $refA := $Gs
        end
    end
end

# %%
@macroexpand @multitrace G 4

# %%
B = @multitrace G 4
A == B

# %%
macro multitrace(A, G, N)
    refA, Gs = refA_Gs(:A, :G, N)
    quote
        let A = $(esc(A)), G = $(esc(G))
            @tensor $refA = $Gs
        end
    end
end

# %%
@macroexpand @multitrace C G 4

# %%
C = similar(zeros(), j...)
@multitrace C G 4
A == B == C

# %% [markdown]
# * https://twitter.com/LirimyDh/status/1559364996865216512
# * https://twitter.com/physics303/status/1559490168477732865

# %%
using StaticArrays

# ３次元配列を行列の配列に変換する
"Just specify that [an element of] the argument is an array of matrices"
struct MatrixArray end

"""
    R = matrixarray(S::AbstractArray{T,3})

`R[k][i, j] = S[i, k, j]`
"""
matrixarray(S) = [SMatrix{size(S[:, i, :])...}(S[:, i, :]) for i in axes(S, 2)]
#matrixarray(S) = [S[:, i, :] for i in axes(S, 2)] # without StaticArrays


# 行列積の最終段を省き、トレースを直接求める
reconst(G) = reconst(MatrixArray(), matrixarray.(G))
reconst(::MatrixArray, G) = _reconst(G...)

# Gs を真ん中で分けないと遅くなる
function _reconst(Gs...)
    h = length(Gs) ÷ 2
    _reconst(reduce(eachprod, Gs[begin:h]), reduce(eachprod, Gs[h+1:end]))
end

#_reconst(G1, G2, Gs...) = _reconst(eachprod(G1, G2), Gs...) # 遅い
_reconst(G1, G2) = trprod.(G1, expanddim(G2, G1))

"""
    trprod(A, B)

Returns `tr(A * B)`
"""
trprod(A, B) = dot(vec(A'), vec(B))

"""
    C = eachprod(A, B)

`C[i, j, k] = A[i, j] * B[k]`

`A, B, C :: Array{Matrix}`
"""
eachprod(A, B) = A .* expanddim(B, A)


"""
    Bx = expanddim(B, A)

`Bx = reshape(B, (1, 1, 1, m, n))` where `ndims(A) == 3`, `size(B) == (m, n)`
"""
expanddim(B, A) = reshape(B, (ntuple(_ -> 1, ndims(A))..., size(B)...))

# %%
D = reconst(G)
A == B == C ≈ D

# %%
@macroexpand reconst4(G) = let X = G; @multitrace X 4 end

# %%
using BenchmarkTools

reconst4(G) = @multitrace G 4

@show size.(G)
@show Threads.nthreads()

D = @btime reconst($G)
E = @btime reconst4($G)
F = @btime @tensor A[j1, j2, j3, j4] := $G[1][i1, j1, i2] * $G[2][i2, j2, i3] * $G[3][i3, j3, i4] * $G[4][i4, j4, i1]
A == B == C == E == F ≈ D

# %%
H = [G[mod1(k-1, 4)] for k in 1:4]

@show size.(H)
@show Threads.nthreads()

X = @btime reconst($H)
Y = @btime reconst4($H)
X ≈ Y

# %%
