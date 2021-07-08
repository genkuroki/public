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
using Random

f(x) = string(first(keys(x))) * " => " * string(first(values(x)))

function test_namedtuple(n)
    S = String[]
    for _ in 1:n
        symb = Symbol(randstring(rand(3:10)))
        str = randstring(rand(3:10))
        nt = NamedTuple{(symb,)}((str,))
        push!(S, f(nt))
    end
    S
end

function test_dict(n)
    S = String[]
    for _ in 1:n
        symb = Symbol(randstring(rand(3:10)))
        str = randstring(rand(3:10))
        dic = Dict(symb => str)
        push!(S, f(dic))
    end
    S
end

@show test_namedtuple(5)
@show test_dict(5)
println()
print("test_namedtuple(2^8):")
@time test_namedtuple(2^8)
print("test_dict(2^8):      ")
@time test_dict(2^8);

# %%
