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
# https://github.com/haampie/ArnoldiMethod.jl

# %%
]st

# %%
using ArnoldiMethod, LinearAlgebra, SparseArrays

# %%
A = spdiagm(
    -1 => fill(-1.0, 99),
    0 => fill(2.0, 100), 
    1 => fill(-1.0, 99)
)

# %%
decomp, history = partialschur(A, nev=10, tol=1e-6, which=SR());

# %%
history

# %%
decomp

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
diag_ = prepareDiag(L)
H = LinearMap(prepareHfunc!(diag_, L), 2^L, ismutating=true, issymmetric=true, isposdef=false)
sol_exact = -2sum(abs(sin((n-1/2) * pi/L)) for n in 1:L)

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
diag_ = prepareDiag(L)
H = LinearMap(prepareHfunc!(diag_, L), 2^L, ismutating=true, issymmetric=true, isposdef=false)
sol_exact = -2sum(abs(sin((n-1/2) * pi/L)) for n in 1:L)

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
diag_ = prepareDiag(L)
H = LinearMap(prepareHfunc!(diag_, L), 2^L, ismutating=true, issymmetric=true, isposdef=false)
sol_exact = -2sum(abs(sin((n-1/2) * pi/L)) for n in 1:L)

end;

# %%
Rev0.sol_exact

# %%
@time decomp0, history0 = partialschur(Rev0.H, nev=1, which=SR())
@time decomp0, history0 = partialschur(Rev0.H, nev=1, which=SR())
decomp0.eigenvalues[1]

# %%
@time decomp1, history1 = partialschur(Rev1.H, nev=1, which=SR())
@time decomp1, history1 = partialschur(Rev1.H, nev=1, which=SR())
decomp1.eigenvalues[1]

# %%
@time decomp2, history2 = partialschur(Rev2.H, nev=1, which=SR())
@time decomp2, history2 = partialschur(Rev2.H, nev=1, which=SR())
decomp2.eigenvalues[1]

# %%
