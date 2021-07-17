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
macro h(x, p...)
    ex = :($(p[end]))
    for i in length(p)-1:-1:1
        ex = :(muladd($(esc(x)), $ex, $(p[i])))
    end
    ex
end

expr = @macroexpand f(x) = (x2 = x^2; x * @h(x2, 1.0, -0.16666666666666666, 0.008333333333333333))
write("def.jl", string(Base.remove_linenums!(expr)))
println(read("def.jl", String)); println()
include("def.jl")
@show f(π/6);

# %%
expr = @macroexpand g(x) = @time (sleep(1); sin(x))
write("foo.jl", string(Base.remove_linenums!(expr)))
println(read("def.jl", String)); println()
include("foo.jl")
@show g(π/6);

# %%
string_remove_linenums!(expr) = string(Base.remove_linenums!(expr))
expr = @macroexpand g(x) = @time (sleep(1); sin(x))
str = string_remove_linenums!(expr)
@show string_remove_linenums!(Meta.parse(str)) == str;

# %%
macro h(x, p...)
    ex = :($(p[end]))
    for i in length(p)-1:-1:1
        ex = :(muladd($(esc(x)), $ex, $(p[i])))
    end
    ex
end

expr = @macroexpand f(x) = (x2 = x^2; x * @h(x2, 1.0, -0.16666666666666666, 0.008333333333333333))
str = string_remove_linenums!(expr)
@show string_remove_linenums!(Meta.parse(str)) == str;

# %%
using Serialization
expr = Base.remove_linenums!(@macroexpand g(x) = @time (sleep(1); sin(x)))
serialize("def.jls", expr)
ast = deserialize("def.jls")
eval(ast)
g(π/6)

# %%
