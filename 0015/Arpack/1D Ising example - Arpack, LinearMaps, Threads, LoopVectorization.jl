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
VERSION

# %%
]st

# %%
Threads.nthreads()

# %%
"""
    Original

* https://gist.github.com/yujitach/c30d7a174bbc3d3d3e40a3c0f9f9d47f
* TAB を "    " で置換
"""
module Original

using LinearAlgebra,LinearMaps
import Arpack

const L=20
    
diag_ = zeros(Float64,2^L)

function prepareDiag(diag)
    for state = 1 : 2^L
        for i = 1 : L
            j = i==L ? 1 : i+1
            diag[state] -= (((state >> (i-1))&1) == ((state >> (j-1))&1)) ? 1 : -1
        end
    end
end
    
function Hfunc!(C,B,diag)
    for state = 1 : 2^L
        C[state] = diag[state] * B[state]
    end
    for state = 1 : 2^L
        for i = 1 : L
            newstate = (state&(~(2^L))) ⊻ (1<<(i-1))
            if newstate==0
                newstate = 2^L
            end
            C[newstate] -= B[state]
        end
    end
end


println("preparing...")
prepareDiag(diag_)

println("computing the lowest eigenvalue...")
H=LinearMap((C,B)->Hfunc!(C,B,diag_),2^L,ismutating=true,issymmetric=true,isposdef=false)
@time e,v = Arpack.eigs(H,nev=1,which=:SR)
@time e,v = Arpack.eigs(H,nev=1,which=:SR)

println("obtained:")
println(e[1])

println("theoretical:")
println(-2sum([ abs(sin((n-1/2) * pi/L)) for n in 1 : L]))

end;

# %%
"""
    Rev0

* This revision is almost equivalent to the original.
* Stop using constants.
* Always pass global variables to functions as arguments.
* Swap the order of the for loop.
* Add @inbounds macro.
* Revise sum([f(x) for x in X]) to sum(f(x) for x in X).
"""
module Rev0

using LinearAlgebra, LinearMaps
import Arpack
    
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
    for state = 1:2^L
        @inbounds C[state] = diag[state] * B[state]
    end
    for i = 1:L
        for state = 1:2^L
            newstate = (state&(~(2^L))) ⊻ (1<<(i-1))
            if newstate == 0
                newstate = 2^L
            end
            @inbounds C[newstate] -= B[state]
        end
    end
end
prepareHfunc!(diag, L) = (C, B) -> Hfunc!(C, B, diag, L)

L = 20

println("preparing...")
diag_ = prepareDiag(L)

println("computing the lowest eigenvalue...")
H = LinearMap(prepareHfunc!(diag_, L), 2^L, ismutating=true, issymmetric=true, isposdef=false)
@time e, v = Arpack.eigs(H, nev=1, which=:SR)
@time e, v = Arpack.eigs(H, nev=1, which=:SR)

println("obtained:")
println(e[1])

println("theoretical:")
println(-2sum(abs(sin((n-1/2) * pi/L)) for n in 1:L))

end;

# %%
"""
    Rev1

* Use Threads.@threads macro.
"""
module Rev1

using LinearAlgebra, LinearMaps
import Arpack
    
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
    Threads.@threads for state = 1:2^L
        @inbounds C[state] = diag[state] * B[state]
    end
    for i = 1:L
        Threads.@threads for state = 1:2^L
            newstate = (state&(~(2^L))) ⊻ (1<<(i-1))
            if newstate == 0
                newstate = 2^L
            end
            @inbounds C[newstate] -= B[state]
        end
    end
end
prepareHfunc!(diag, L) = (C, B) -> Hfunc!(C, B, diag, L)

L = 20

println("preparing...")
diag_ = prepareDiag(L)

println("computing the lowest eigenvalue...")
H = LinearMap(prepareHfunc!(diag_, L), 2^L, ismutating=true, issymmetric=true, isposdef=false)
@time e, v = Arpack.eigs(H, nev=1, which=:SR)
@time e, v = Arpack.eigs(H, nev=1, which=:SR)

println("obtained:")
println(e[1])

println("theoretical:")
println(-2sum(abs(sin((n-1/2) * pi/L)) for n in 1:L))

end;

# %%
"""
    Rev2

* Use LoopVectorization.@tturbo macro.
"""
module Rev2

using LinearAlgebra, LinearMaps
import Arpack
using LoopVectorization

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

L = 20

println("preparing...")
diag_ = prepareDiag(L)

println("computing the lowest eigenvalue...")
H = LinearMap(prepareHfunc!(diag_, L), 2^L, ismutating=true, issymmetric=true, isposdef=false)
@time e, v = Arpack.eigs(H, nev=1, which=:SR)
@time e, v = Arpack.eigs(H, nev=1, which=:SR)

println("obtained:")
println(e[1])

println("theoretical:")
println(-2sum(abs(sin((n-1/2) * pi/L)) for n in 1:L))

end;

# %%
using BenchmarkTools
using LinearAlgebra
using Arpack: Arpack

H = Original.H
H0 = Rev0.H
H1 = Rev1.H
H2 = Rev2.H

B, C = similar(Original.diag_), similar(Original.diag_)
println("Hamiltonian bench")
print("  Original:                ")
@btime mul!($C, $H, $B)
print("  Rev0 (almost original):  ")
@btime mul!($C, $H0, $B)
print("  Rev1 (Threads.@threads): ")
@btime mul!($C, $H1, $B)
print("  Rev2 (LoopVectorization):")
@btime mul!($C, $H2, $B);

# %%
println("Arpack.eigs bench")
print("  Original:                ")
@btime e, v = Arpack.eigs($H, nev=1, which=:SR)
print("  Rev0 (almost original):  ")
@btime e, v = Arpack.eigs($H0, nev=1, which=:SR)
print("  Rev1 (Threads.@threads): ")
@btime e, v = Arpack.eigs($H1, nev=1, which=:SR)
print("  Rev2 (LoopVectorization):")
@btime e, v = Arpack.eigs($H2, nev=1, which=:SR);

# %%
@show a = 2
@show foo = x -> a*x
@show foo(3)
@show a = 10
@show foo(3);

# %%
@show a = 2
@show makebar(a) = x -> a*x
@show bar = makebar(a)
@show bar(3)
@show a = 10
@show bar(3);

# %%
