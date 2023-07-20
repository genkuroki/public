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
#     display_name: Julia 1.9.2
#     language: julia
#     name: julia-1.9
# ---

# %% [markdown]
# https://theconversation.com/momentum-isnt-magic-vindicating-the-hot-hand-with-the-mathematics-of-streaks-74786

# %%
using Distributions
using Random

function prob(X, k)
    n = length(X)
    den = 0
    num = 0
    for i in k+1:n
        if all(X[i-j] for j in 1:k)
            den += 1
            num += X[i]
        end
    end
    num/den
end

function sim(n, k; L=10^7)
    nths = Threads.nthreads()
    Xtmp = [trues(n) for _ in 1:nths]
    P = zeros(nths)
    N = zeros(Int, nths)
    Threads.@threads for _ in 1:L
        tid = Threads.threadid()
        X = rand!(Xtmp[tid])
        p = prob(X, k)
        if !isnan(p)
            N[tid] += 1
            P[tid] += p
        end
    end
    sum(N)/L, sum(P)/sum(N)
end

function sim1(n, k; L=10^7)
    nths = Threads.nthreads()
    Xtmp = [trues(n) for _ in 1:nths]
    M = zeros(Int, nths)
    N = zeros(Int, nths)
    Threads.@threads for _ in 1:L
        tid = Threads.threadid()
        X = rand!(Xtmp[tid])
        i = rand(k+1:n)
        if all(X[i-j] for j in 1:k)
            N[tid] += 1
            M[tid] += X[i]
        end
    end
    sum(N)/L, sum(M)/sum(N)
end

function isadmissible(X, k)
    c = 0
    for x in X[begin:end-1]
        if x
            c += 1
            c ≥ k && return true
        else
            c = 0
        end
    end
    false
end

function sim2(n, k; L=10^7)
    nths = Threads.nthreads()
    Xtmp = [trues(n) for _ in 1:nths]
    M = zeros(Int, nths)
    N = zeros(Int, nths)
    Threads.@threads for _ in 1:L
        tid = Threads.threadid()
        X = rand!(Xtmp[tid])
        isadmissible(X, k) || continue
        N[tid] += 1
        c = 0
        while true
            c += 1
            i = rand(k+1:n)
            if all(X[i-j] for j in 1:k)
                M[tid] += X[i]
                break
            end
        end
    end
    sum(N)/L, sum(M)/sum(N)
end

function alladmissibles(n, k)
    @assert n ≤ 20
    A = Tuple{NTuple{n, Bool}, Vector{Bool}}[]
    for X in Iterators.product(fill((false, true), n)...)
        if isadmissible(X, k)
            Y = Bool[]
            for i in 1+k:n
                if all(X[i-j] for j in 1:k)
                    push!(Y, X[i])
                end
            end
            push!(A, (X, Y))
        end
    end
    A
end

# %%
alladmissibles(5, 3)

# %%
alladmissibles(5, 3) .|> (t -> t[2]) .|> mean |> mean

# %%
@show sim(5, 3) sim1(5, 3) sim2(5, 3);

# %%
alladmissibles(10, 3)

# %%
alladmissibles(10, 3) .|> (t -> t[2]) .|> mean |> mean

# %%
@show sim(10, 3) sim1(10, 3) sim2(10, 3);

# %%
@show sim(30, 3) sim1(30, 3) sim2(30, 3);

# %%
@show sim(100, 3) sim1(100, 3) sim2(100, 3);

# %%
@show sim(300, 3) sim1(300, 3) sim2(300, 3);

# %%
@show sim(1000, 3) sim1(1000, 3) sim2(1000, 3);

# %%
