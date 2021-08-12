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

# %% [markdown]
# https://twitter.com/yujitach/status/1424030835771023363

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

f_LinearMap(20), f_exact(20)

# %%
@time f_LinearMap(20)
@time f_LinearMap(20)
@time f_LinearMap(20)

# %% [markdown]
# https://juliaphysics.github.io/PhysicsTutorials.jl/tutorials/general/quantum_ising/quantum_ising.html

# %%
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
    
    H = spzeros(Int, 2^N, 2^N)
    for i in 1:N
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

f_exact(L) = -2sum(abs(sin((n-1/2) * pi/L)) for n in 1:L)

f_sparse(20), f_exact(20)

# %%
L = 20
@time H = TransverseFieldIsing_sparse(N=L, h=1)
@time d, h = partialschur(H; nev=1, which=SR())

# %%
L = 20
@time H = TransverseFieldIsing_sparse(N=L, h=1)
@time d, h = partialschur(H; nev=1, which=SR())

# %%
L = 20
@time H = TransverseFieldIsing_sparse(N=L, h=1)
@time d, h = partialschur(H; nev=1, which=SR())

# %%
TransverseFieldIsing_sparse(N=20, h=1)

# %%
function TransverseFieldIsing_sparse_revised(;N, h)
    id = [1 0; 0 1] |> sparse
    σˣ = [0 1; 1 0] |> sparse
    σᶻ = [1 0; 0 -1] |> sparse
    
    first_term_ops = fill(id, N)
    first_term_ops[1] = σᶻ
    first_term_ops[2] = σᶻ
    
    second_term_ops = fill(id, N)
    second_term_ops[1] = σˣ
    
    tmp = [spzeros(Int, 2^k, 2^k) for k in 1:N]
    H = spzeros(Int, 2^N, 2^N)
    
    for i in 1:N
        tmp[1] .= first_term_ops[mod1(i, N)]
        for k in 2:N
            kron!(tmp[k], tmp[k-1], first_term_ops[mod1(i+k-1, N)])
        end
        H .-= tmp[N]
    end
    
    for i in 1:N
        tmp[1] .= second_term_ops[mod1(i, N)]
        for k in 2:N
            kron!(tmp[k], tmp[k-1], second_term_ops[mod1(i+k-1, N)])
        end
        H .-= h .* tmp[N]
    end
    H
end

function f_sparse_revised(L)
    H = TransverseFieldIsing_sparse_revised(N=L, h=1)
    d, h = partialschur(H; nev=1, which=SR())
    d.eigenvalues[1]
end

TransverseFieldIsing_sparse_revised(N=10, h=1) == TransverseFieldIsing_sparse(N=10, h=1)

# %%
L = 20
@time H = TransverseFieldIsing_sparse_revised(N=L, h=1)
@time d, h = partialschur(H; nev=1, which=SR())

# %%
L = 20
@time H = TransverseFieldIsing_sparse_revised(N=L, h=1)
@time d, h = partialschur(H; nev=1, which=SR())

# %%
L = 20
@time H = TransverseFieldIsing_sparse_revised(N=L, h=1)
@time d, h = partialschur(H; nev=1, which=SR())

# %%
TransverseFieldIsing_sparse_revised(N=20, h=1)

# %%
bit(a, k) = (a >> k) & 1

function sigmaz2(N)
    d = zeros(Int, 2^N)
    for k in 1:N, a in 0:2^N-1
        d[a+1] -= ifelse(bit(a, N-k) == bit(a, N - mod1(k+1, N)), 1, -1)
    end
    d
end

function sigmax(N, k)
    v = zeros(Int, 2^N - 2^(N-k))
    for i in 1:2:2^k-1
        v[(2^(N-k)*(i-1) + 1):2^(N-k)*i] .= -1
    end
    v
end

function TransverseFieldIsing_sparse_revised2(;N, h)
    d = sigmaz2(N)    
    v = sigmax.(N, 1:N)    
    H = spdiagm(
        (-2^(N-k) => v[k] for k in 1:N)...,
        0 => d,
        ( 2^(N-k) => v[k] for k in 1:N)...)
end

function f_sparse_revised2(L)
    H = TransverseFieldIsing_sparse_revised2(N=L, h=1)
    d, h = partialschur(H; nev=1, which=SR())
    d.eigenvalues[1]
end

TransverseFieldIsing_sparse_revised2(N=10, h=1) == TransverseFieldIsing_sparse(N=10, h=1)

# %%
L = 20
@time H = TransverseFieldIsing_sparse_revised2(N=L, h=1)
@time d, h = partialschur(H; nev=1, which=SR())

# %%
L = 20
@time H = TransverseFieldIsing_sparse_revised2(N=L, h=1)
@time d, h = partialschur(H; nev=1, which=SR())

# %%
L = 20
@time H = TransverseFieldIsing_sparse_revised2(N=L, h=1)
@time d, h = partialschur(H; nev=1, which=SR())

# %%
function TransverseFieldIsing_sparse_revised2_time(;N, h)
    @time d = sigmaz2(N)    
    @time v = sigmax.(N, 1:N)    
    @time H = spdiagm(
        (-2^(N-k) => v[k] for k in 1:N)...,
        0 => d,
        ( 2^(N-k) => v[k] for k in 1:N)...)
end

TransverseFieldIsing_sparse_revised2_time(N=20, h=1); println()
TransverseFieldIsing_sparse_revised2_time(N=20, h=1); println()
TransverseFieldIsing_sparse_revised2_time(N=20, h=1)

# %%
bit(a, k) = (a >> k) & 1

function sigmaz2_sparse(N)
    d = zeros(Int, 2^N)
    for k in 1:N, a in 0:2^N-1
        d[a+1] -= ifelse(bit(a, N-k) == bit(a, N - mod1(k+1, N)), 1, -1)
    end
    sparse(d) # should be sparse
end

function sigmax_sparse(N, k)
    v = zeros(Int, 2^N - 2^(N-k))
    for i in 1:2:2^k-1
        v[(2^(N-k)*(i-1) + 1):2^(N-k)*i] .= -1
    end
    sparse(v) # should be sparse
end

function TransverseFieldIsing_sparse_revised3(;N, h)
    d = sigmaz2_sparse(N)    
    v = sigmax_sparse.(N, 1:N)    
    H = spdiagm(
        (-2^(N-k) => v[k] for k in 1:N)...,
        0 => d,
        ( 2^(N-k) => v[k] for k in 1:N)...)
end

function f_sparse_revised3(L)
    H = TransverseFieldIsing_sparse_revised3(N=L, h=1)
    d, h = partialschur(H; nev=1, which=SR())
    d.eigenvalues[1]
end

TransverseFieldIsing_sparse_revised3(;N=10, h=1) == TransverseFieldIsing_sparse(;N=10, h=1)

# %%
L = 20
@time H = TransverseFieldIsing_sparse_revised3(N=L, h=1)
@time d, h = partialschur(H; nev=1, which=SR())

# %%
L = 20
@time H = TransverseFieldIsing_sparse_revised3(N=L, h=1)
@time d, h = partialschur(H; nev=1, which=SR())

# %%
L = 20
@time H = TransverseFieldIsing_sparse_revised3(N=L, h=1)
@time d, h = partialschur(H; nev=1, which=SR())

# %%
function TransverseFieldIsing_sparse_revised3_time(;N, h)
    @time d = sigmaz2_sparse(N)    
    @time v = sigmax_sparse.(N, 1:N)    
    @time H = spdiagm(
        (-2^(N-k) => v[k] for k in 1:N)...,
        0 => d,
        ( 2^(N-k) => v[k] for k in 1:N)...)
end

TransverseFieldIsing_sparse_revised3_time(N=20, h=1); println()
TransverseFieldIsing_sparse_revised3_time(N=20, h=1); println()
TransverseFieldIsing_sparse_revised3_time(N=20, h=1)

# %%
