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

# %% [markdown]
# https://github.com/genkuroki/public/blob/main/0036/TensorOperations.jl.ipynb の続き

# %%
using TensorOperations
using BenchmarkTools

# %%
function refA_Gs(A, G, N)
    Gs = :($G[1][i1, j1, i2])
    for k in 2:N
        ik = Symbol(:i, k)
        jk = Symbol(:j, k)
        ikp1 = Symbol(:i, mod1(k+1, N))
        Gs = :($Gs * $G[$k][$ik, $jk, $ikp1])
    end
    refA = Expr(:ref, :A, Symbol.(:j, 1:N)...)
    refA, Gs
end

macro multitrace(G, N)
    refA, Gs = refA_Gs(:A, :G, N)
    quote
        let G = $(esc(G))
            @tensor $refA := $Gs
        end
    end
end

macro multitrace(A, G, N)
    refA, Gs = refA_Gs(:A, :G, N)
    quote
        let A = $(esc(A)), G = $(esc(G))
            @tensor $refA = $Gs
        end
    end
end

# %%
@generated function multr(G, ::Val{N}) where N
    refA, Gs = refA_Gs(:A, :G, N)
    quote
        @tensor $refA := $Gs
    end
end

multr(G) = multr(G, Val(length(G)))

@generated function multr!(A, G, ::Val{N}) where N
    refA, Gs = refA_Gs(:A, :G, N)
    quote
        @tensor $refA = $Gs
    end
end

multr!(A, G) = multr!(A, G, Val(length(G)))

# %%
i = [100, 8, 6, 5]
j = [6, 7, 5, 4]
H = [randn(i[mod1(k, 4)], j[k], i[mod1(k+1, 4)]) for k in 1:4]
size.(H)

# %%
A = @multitrace H 4
typeof(A), size(A)

# %%
B = multr(H)
A == B

# %%
C = similar(zeros(), j...)
multr!(C, H)
A == B == C

# %%
@show size.(H)

A = @btime @multitrace $H 4
B = @btime multr($H)
C = @btime multr!($C, $H)
A == B == C

# %%
K = [H[mod1(k-1, length(H))] for k in 1:length(H)]
@show size.(K)
F = similar(zeros(), size.(K, 2)...)

D = @btime @multitrace $K 4
E = @btime multr($K)
F = @btime multr!($F, $K)
D == E == F

# %%
permutedims(A, (4, 1, 2, 3)) ≈ D

# %%
