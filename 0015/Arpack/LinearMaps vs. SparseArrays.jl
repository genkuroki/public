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
prepareHfunc!(L) = (C, B) -> Hfunc!(C, B, prepareDiag(L), L)

function f_LinearMap(L)
    H_LinearMap = LinearMap(prepareHfunc!(L), 2^L, ismutating=true, issymmetric=true, isposdef=false)
    d_LM, h_LM = partialschur(H_LinearMap, nev=1, which=SR())
    d_LM.eigenvalues[1]
end

f_exact(L) = -2sum(abs(sin((n-1/2) * pi/L)) for n in 1:L)

# %%
# https://juliaphysics.github.io/PhysicsTutorials.jl/tutorials/general/quantum_ising/quantum_ising.html

using LinearAlgebra
using SparseArrays
using ArnoldiMethod
⊗(x,y) = kron(x,y)

function TransverseFieldIsing_sparse(;N,h)
    id = [1 0; 0 1] |> sparse
    σˣ = [0 1; 1 0] |> sparse
    σᶻ = [1 0; 0 -1] |> sparse
    
    first_term_ops = fill(id, N)
    first_term_ops[1] = σᶻ
    first_term_ops[2] = σᶻ
    
    second_term_ops = fill(id, N)
    second_term_ops[1] = σˣ
    
    H = spzeros(Int, 2^N, 2^N) # note the spzeros instead of zeros here
    #for i in 1:N-1
    for i in 1:N # periodic boundary condition
        H -= foldl(⊗, first_term_ops)
        first_term_ops = circshift(first_term_ops,1)
    end
    
    for i in 1:N
        H -= h*foldl(⊗, second_term_ops)
        second_term_ops = circshift(second_term_ops,1)
    end
    H
end

function f_sparse(L)
    H = TransverseFieldIsing_sparse(N=L, h=1)
    d, h = partialschur(H; nev=1, which=SR())
    d.eigenvalues[1]
end

# %%
f_exact(20)

# %%
@time f_LinearMap(20)
@time f_LinearMap(20)

# %%
@time f_sparse(20)
@time f_sparse(20)

# %%
@time H = TransverseFieldIsing_sparse(N=20, h=1)
@time d, h = partialschur(H; nev=1, which=SR())

# %%
@time f_LinearMap(10)
@time f_LinearMap(10)

# %%
@time f_sparse(10)
@time f_sparse(10)

# %%
L = 10
@time H_LM = LinearMap(prepareHfunc!(L), 2^L, ismutating=true, issymmetric=true, isposdef=false)
@time d_LM, h_LM = partialschur(H_LM, nev=1, which=SR())
@time d_LM, h_LM = partialschur(H_LM, nev=1, which=SR())

# %%
@time H_sparse = TransverseFieldIsing_sparse(N=10, h=1)
@time d, h = partialschur(H_sparse; nev=1, which=SR())
@time d, h = partialschur(H_sparse; nev=1, which=SR())

# %%
using BenchmarkTools

L = 10
@show L
println("LinearMaps:")
H_LM = @btime LinearMap(prepareHfunc!(L), 2^L, ismutating=true, issymmetric=true, isposdef=false)
@btime d_LM, h_LM = partialschur($H_LM, nev=1, which=$(SR()))
println("SparseArrays:")
H_sparse = @btime TransverseFieldIsing_sparse(N=L, h=1)
@btime d, h = partialschur($H_sparse; nev=1, which=$(SR()));

# %%
using BenchmarkTools

L = 20
@show L
println("LinearMaps:")
H_LM = @btime LinearMap(prepareHfunc!(L), 2^L, ismutating=true, issymmetric=true, isposdef=false)
@btime d_LM, h_LM = partialschur($H_LM, nev=1, which=$(SR()))
println("SparseArrays:")
H_sparse = @btime TransverseFieldIsing_sparse(N=L, h=1)
@btime d, h = partialschur($H_sparse; nev=1, which=$(SR()));

# %%
