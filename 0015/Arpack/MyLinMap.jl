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
using LinearAlgebra
using LoopVectorization
using LinearMaps
using Arpack: Arpack
using ArnoldiMethod

function prepareDiag(L)
    diag = zeros(2^L)
    for state = 1:2^L
        for i = 1:L
            j = i==L ? 1 : i+1
            @inbounds diag[state] -= (((state >> (i-1))&1) == ((state >> (j-1))&1)) ? 1 : -1
        end
    end
    diag
end
    
function Hfunc!(C, B, diag, L)
    N = length(diag)
    @tturbo for state = 1:N
        C[state] = diag[state] * B[state]
    end
    for i = 1:L
        @tturbo for state = 1:N
            newstate = (state&(~(2^L))) ⊻ (1<<(i-1))
            c = newstate == 0
            newstate = !c*newstate + c*N # remove if statement
            C[newstate] -= B[state]
        end
    end
end
prepareHfunc!(diag, L) = (C, B) -> Hfunc!(C, B, diag, L)

struct MyLinMap{T, F, I<:Integer}
    f!::F; N::I; issymmetric::Bool; isposdef::Bool
end
MyLinMap(f!, N; issymmetric=false, isposdef=false) =
    MyLinMap{Float64, typeof(f!), typeof(N)}(f!, N, issymmetric, isposdef)
LinearAlgebra.mul!(y, A::MyLinMap, x) = A.f!(y, x)
Base.size(A::MyLinMap) = (A.N, A.N)
Base.size(A::MyLinMap, i::Integer) = size(A::MyLinMap)[i]
Base.eltype(A::MyLinMap{T}) where T = T
LinearAlgebra.issymmetric(A::MyLinMap) = A.issymmetric
LinearAlgebra.isposdef(A::MyLinMap) = A.isposdef

L = 20
diag_ = prepareDiag(L)
H! = prepareHfunc!(diag_, L)

H_LinearMap = LinearMap(H!, 2^L, ismutating=true, issymmetric=true, isposdef=false)
H_MyLinMap = MyLinMap(H!, 2^L, issymmetric = true, isposdef = false)
sol_exact = -2sum(abs(sin((n-1/2) * pi/L)) for n in 1:L)
@show sol_exact;

# %%
@time e1, v1 = Arpack.eigs(H_LinearMap, nev=1, which=:SR)
@time e1, v1 = Arpack.eigs(H_LinearMap, nev=1, which=:SR)
e1[1]

# %%
@time e2, v2 = Arpack.eigs(H_MyLinMap, nev=1, which=:SR)
@time e2, v2 = Arpack.eigs(H_MyLinMap, nev=1, which=:SR)
e2[1]

# %%
@time d1, h1 = partialschur(H_LinearMap, nev=1, which=SR())
@time d1, h1 = partialschur(H_LinearMap, nev=1, which=SR())
d1.eigenvalues[1]

# %%
@time d2, h2 = partialschur(H_MyLinMap, nev=1, which=SR())
@time d2, h2 = partialschur(H_MyLinMap, nev=1, which=SR())
d2.eigenvalues[1]

# %%
