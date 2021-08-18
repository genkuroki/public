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
#     display_name: Julia 1.6.2
#     language: julia
#     name: julia-1.6
# ---

# %%
?@generated

# %%
function _staticsum(f, a, b)
    expr = :(f($a))
    for i in a+1:b
        expr = :($expr + f($i))
    end
    expr
end

# %%
_staticsum(:f, -3, 4)

# %%
@generated function staticsum(f, ::Val{a}, ::Val{b}) where {a, b}
    @assert b - a < 100
    expr = :(f($a))
    for i in a+1:b
        expr = :($expr + f($i))
    end
    expr
end

# %%
staticsum(identity, Val(1), Val(10))

# %%
2staticsum(x -> sin(x)/x, Val(1), Val(100)) + 1

# %%
F(n) = 2sum(x -> sin(x)/x, 1:n) + 1
F_static(n) = 2staticsum(x -> sin(x)/x, Val(1), Val(n)) + 1

# %%
using BenchmarkTools
@show F(100) == F_static(100)
@btime F(100)
@btime F_static(100)

# %%
@benchmark F(100)

# %%
@benchmark F_static(100)

# %%
@code_warntype staticsum(x -> sin(x)/x, Val(1), Val(5))

# %%
@code_typed staticsum(x -> sin(x)/x, Val(1), Val(5))

# %%
@code_typed sum(x -> sin(x)/x, 1:5)

# %%
