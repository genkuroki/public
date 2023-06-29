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
#     display_name: Julia 1.9.1
#     language: julia
#     name: julia-1.9
# ---

# %% [markdown]
# * https://twitter.com/dannchu/status/1672880182427004929
# * https://twitter.com/aoki_taichi/status/1673834965694566401

# %%
n=10^8 #1億回
k=0
@time for i=1:n
    a1=randn(2);a2=randn(2);a3=randn(2) #3点ランダムにとる。
    x1=a2-a1;x2=a3-a1;t1=(x1[1]*x2[1]+x1[2]*x2[2])/(sqrt(x1[1]^2+x1[2]^2)*sqrt(x2[1]^2+x2[2]^2))# 1つ目の角のcos
    x1=a1-a2;x2=a3-a2;t2=(x1[1]*x2[1]+x1[2]*x2[2])/(sqrt(x1[1]^2+x1[2]^2)*sqrt(x2[1]^2+x2[2]^2))# 2つ目の角のcos
    x1=a1-a3;x2=a2-a3;t3=(x1[1]*x2[1]+x1[2]*x2[2])/(sqrt(x1[1]^2+x1[2]^2)*sqrt(x2[1]^2+x2[2]^2))# 3つ目の角のcos
    if t1>0.0 && t2>0.0 && t3>0.0 # 3つの角が鋭角という条件
        k+=1 # 条件を満たしたらカウントする。
    end
end
println(k*100/n,"%")

# %%
function f(n = 10^8)
    k=0
    for i=1:n
        a1=randn(2);a2=randn(2);a3=randn(2) #3点ランダムにとる。
        x1=a2-a1;x2=a3-a1;t1=(x1[1]*x2[1]+x1[2]*x2[2])/(sqrt(x1[1]^2+x1[2]^2)*sqrt(x2[1]^2+x2[2]^2))# 1つ目の角のcos
        x1=a1-a2;x2=a3-a2;t2=(x1[1]*x2[1]+x1[2]*x2[2])/(sqrt(x1[1]^2+x1[2]^2)*sqrt(x2[1]^2+x2[2]^2))# 2つ目の角のcos
        x1=a1-a3;x2=a2-a3;t3=(x1[1]*x2[1]+x1[2]*x2[2])/(sqrt(x1[1]^2+x1[2]^2)*sqrt(x2[1]^2+x2[2]^2))# 3つ目の角のcos
        if t1>0.0 && t2>0.0 && t3>0.0 # 3つの角が鋭角という条件
            k+=1 # 条件を満たしたらカウントする。
        end
    end
    println(k*100/n,"%")
end
@time f()

# %%
using Random
using StaticArrays

function f1(n = 10^8)
    a1, a2, a3 = (MVector{2, Float64}(undef) for _ in 1:3)
    k = 0
    for i in 1:n
        randn!(a1); randn!(a2); randn!(a3) #3点ランダムにとる。
        x1=a2-a1; x2=a3-a1; t1=(x1[1]*x2[1]+x1[2]*x2[2])/(sqrt(x1[1]^2+x1[2]^2)*sqrt(x2[1]^2+x2[2]^2)) # 1つ目の角のcos
        x1=a1-a2; x2=a3-a2; t2=(x1[1]*x2[1]+x1[2]*x2[2])/(sqrt(x1[1]^2+x1[2]^2)*sqrt(x2[1]^2+x2[2]^2)) # 2つ目の角のcos
        x1=a1-a3; x2=a2-a3; t3=(x1[1]*x2[1]+x1[2]*x2[2])/(sqrt(x1[1]^2+x1[2]^2)*sqrt(x2[1]^2+x2[2]^2)) # 3つ目の角のcos
        if t1>0.0 && t2>0.0 && t3>0.0 # 3つの角が鋭角という条件
            k+=1 # 条件を満たしたらカウントする。
        end
    end
    100k/n
end
@time f1()
@time f1()
@time f1()

# %%
using LinearAlgebra
using Random
using StaticArrays

function f2(n = 10^8)
    a, b, c = (MVector{2, Float64}(undef) for _ in 1:3)
    k = 0
    for i in 1:n
        randn!(a); randn!(b); randn!(c)
        dot(b - a, c - a) ≤ 0 && continue
        dot(c - b, a - b) ≤ 0 && continue
        dot(a - c, b - c) ≤ 0 && continue
        k += 1
    end
    100k/n
end
@time f2()
@time f2()
@time f2()

# %%
using LinearAlgebra
using Random
using StaticArrays

function f3(n = 10^8)
    a, b, c = (MVector{2, Float64}(undef) for _ in 1:3)
    k = 0
    for i in 1:n
        randn!(a); randn!(b); randn!(c)
        k += (dot(b - a, c - a) > 0) & (dot(c - b, a - b) > 0) & (dot(a - c, b - c) > 0)
    end
    100k/n
end
@time f3()
@time f3()
@time f3()

# %%
using LinearAlgebra
using Random
using StaticArrays

@show ENV["JULIA_NUM_THREADS"]
@show Threads.nthreads();

function f2threads(n = 10^8)
    A, B, C = ([MVector{2, Float64}(undef) for _ in 1:Threads.nthreads()] for _ in 1:3)
    K = zeros(Int, Threads.nthreads())
    Threads.@threads for i in 1:n
        tid = Threads.threadid()
        a, b, c = A[tid], B[tid], C[tid]
        randn!(a); randn!(b); randn!(c)
        dot(b - a, c - a) ≤ 0 && continue
        dot(c - b, a - b) ≤ 0 && continue
        dot(a - c, b - c) ≤ 0 && continue
        K[tid] += 1
    end
    100sum(K)/n
end
@time f2threads()
@time f2threads()
@time f2threads()

# %%
using LinearAlgebra
using Random
using StaticArrays

@show ENV["JULIA_NUM_THREADS"]
@show Threads.nthreads();

function f3threads(n = 10^8)
    A, B, C = ([MVector{2, Float64}(undef) for _ in 1:Threads.nthreads()] for _ in 1:3)
    K = zeros(Int, Threads.nthreads())
    Threads.@threads for i in 1:n
        tid = Threads.threadid()
        a, b, c = A[tid], B[tid], C[tid]
        randn!(a); randn!(b); randn!(c)
        K[tid] += (dot(b - a, c - a) > 0) & (dot(c - b, a - b) > 0) & (dot(a - c, b - c) > 0)
    end
    100sum(K)/n
end
@time f3threads()
@time f3threads()
@time f3threads()

# %%
