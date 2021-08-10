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
#     display_name: Julia 1.8.0-DEV -t 6
#     language: julia
#     name: julia-1.8-t4
# ---

# %%
a = Threads.Atomic{Int}(2)
b = a[]

# %%
using Printf
@show Threads.nthreads()
@show a = Threads.Atomic{Int}(99)
Threads.@threads for i in 1:10
    a[] = i
    @printf "(id: %2d)  i = %2d,  a[] = %2d\n" Threads.threadid() i a[]
    b = a[]
    @printf "(id: %2d)  i = %2d,  b = %2d\n" Threads.threadid() i b
end

# %%
