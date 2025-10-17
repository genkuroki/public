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
#     display_name: Julia
#     language: julia
#     name: julia
# ---

# %%
@show VERSION

using Random, Statistics
mynthreads() = Threads.nthreads() + Threads.nthreads(:interactive)

n, L = 10, 30
nth = mynthreads()
Xtmp = [rand(n) for _ in 1:nth]
M = Vector{Tuple{Int64, Int64}}(undef, L)

for i in 1:L
#Threads.@threads for i in 1:L
    tid = Threads.threadid()
    X = rand!(Xtmp[tid])
    M[i] = (tid, round(Int, 100mean(X)))
end
@show M;

#for i in 1:L
Threads.@threads for i in 1:L
    tid = Threads.threadid()
    X = rand!(Xtmp[tid])
    M[i] = (tid, round(Int, 100mean(X)))
end
@show M;

# %%
@show VERSION

using Random, Statistics
mynthreads() = Threads.nthreads() + Threads.nthreads(:interactive)

n, L = 10, 30
nth = mynthreads()
Xtmp = [rand(n) for _ in 1:nth]
M = Vector{Tuple{Int64, Int64}}(undef, L)

for i in 1:L
#Threads.@threads for i in 1:L
    tid = Threads.threadid()
    X = rand!(Xtmp[tid])
    M[i] = (tid, round(Int, 100mean(X)))
end
@show M;

#for i in 1:L
Threads.@threads for i in 1:L
    tid = Threads.threadid()
    X = rand!(Xtmp[tid])
    M[i] = (tid, round(Int, 100mean(X)))
end
@show M;

# %%
