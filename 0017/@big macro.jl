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
#     display_name: Julia 1.6.2
#     language: julia
#     name: julia-1.6
# ---

# %%
"""
    big!(expr)

replaces `Int` and `Float64` values in `expr` by big ones and returns the result.
"""
big!(expr) = expr
big!(expr::Int) = big(expr)
big!(expr::Float64) = big(expr)
function big!(expr::Expr)
    for i in eachindex(expr.args)
        expr.args[i] = big!(expr.args[i])
    end
    expr
end

"""
    @big(expr)
    @big expr

replaces `Int` and `Float64` values in `expr` by big ones and executes the result.
"""
macro big(expr) big!(expr) end

# %%
@doc big!

# %%
@doc @big

# %%
1.0π

# %%
@big 1.0π

# %%
for n in 20:24
    println("factorial(", n, ") = ", factorial(n))
end

# %%
@big for n in 20:24
    println("factorial(", n, ") = ", factorial(n))
end

# %%
@big @eval begin
    n = 57
    c = Tuple(1/factorial(k) for k in 0:n)
end
evalpoly(1, c)

# %%
@big exp(1)

# %%
typeof(c)

# %%
