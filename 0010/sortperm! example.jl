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
# https://twitter.com/physics303/status/1414522243216859143?s=21

# %%
f(A) = sortperm(sort.(collect(eachrow(A))); rev=true)

function g!(A, ix, tmp)
    for i in axes(A, 1) tmp[i] .= @view A[i, :] end
    sort!.(tmp)
    sortperm!(ix, tmp; rev=true)
end

function maketmp(A)
    ix = Vector{Int}(undef, size(A, 1))
    tmp = [Vector{eltype(A)}(undef, size(A, 2)) for _ in 1:size(A, 1)]
    ix, tmp
end

B = [
    2 4 4 4 5 9
    1 2 2 3 5 9
    1 2 3 3 3 9
    1 2 2 5 6 9
    1 3 4 4 4 9
]

ix, tmp = maketmp(B)
@show ix
@show tmp
@show f(B) == g!(B, ix, tmp)
@time f(B)
@time g!(B, ix, tmp)

# %%
using BenchmarkTools
@btime f($B)
@btime g!($B, $ix, $tmp)

# %%
