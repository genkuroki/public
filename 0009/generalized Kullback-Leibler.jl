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
using StatsFuns
gkl_gen(Q, P) = sum(xlogy(p, p/q) for (p, q) in zip(P, Q)) + sum(Q) - sum(P)

# %%
n = 10^7
P = rand(n)
Q = abs.(randn(n))
gkl_gen(Q, P), gkl_gen(P, Q)

# %%
gkl_gen(Q, Q)

# %%
gkl_gen(reverse(Q), Q)

# %%
using BenchmarkTools
@benchmark gkl_gen(Q, P)

# %%
using LoopVectorization

unsafe_xlogy(x, y) = x*log(y)
safe_xlogy(x, y) = iszero(x) ? x : x*log(y)

function gkl_lv(Q, P)
    s = sum(Q) - sum(P)
    @tturbo for i in eachindex(P, Q)
        s += unsafe_xlogy(P[i], P[i]/Q[i])
    end
    s
end

gkl_lv(Q, P)

# %%
@benchmark gkl_lv(Q, P)

# %%
Threads.nthreads()

# %%
