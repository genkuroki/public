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
fib_l = 
((:lambda, (:u,), (:u, :u)),
    (:lambda, (:u,),
        (:lambda, (:n, :a, :b),
            (:if, (:(==), :n, 0), :a,
                ((:u, :u), (:-, :n, 1), :b, (:+, :a, :b))))))

l2j(s) = s
function l2j(s::Tuple)
    s[1] === :if     && return Expr(l2j.(s)...)
    s[1] === :lambda && return Expr(:(->), Expr(:tuple, s[2]...), l2j(s[3]))
    Expr(:call, l2j.(s)...)
end
macro l(x) l2j(Core.eval(__module__, x)) end

fib= @l fib_l
@show methods(fib).ms[1]
@show fib(21, 0, 1);

# %%
f_l =
((:lambda, (:u,), (:u, :u)),
    (:lambda, (:u,),
        (:lambda, (:n,),
            (:if, (:(≤), :n, 0), (:vcat,),
                (:vcat, ((:u, :u), (:-, :n, 1)), :n)))))

f = @l f_l
f(10.0)

# %%
