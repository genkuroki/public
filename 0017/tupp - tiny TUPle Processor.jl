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
t2j(s) = s
function t2j(s::Tuple)
    a = (:block, :if, :tuple, :parameters, :kw, :quote, :(=), :(.=), :.)
    s[1] in a && return Expr(s[1], t2j.(s[2:end])...)
    s[1] === :lambda && return Expr(:(->), Expr(:tuple, s[2]...), t2j(s[3]))
    s[1] === :for && return Expr(s[1], Expr(:(=), t2j(s[2]), t2j(s[3])), t2j(s[4]))
    s[1] === :using && return Expr(s[1], Expr.(:., s[2:end])...)
    Expr(:call, t2j.(s)...)
end
macro tupp(x) t2j(Core.eval(__module__, x)) end

# %%
(:for, :k, (:(:), 1, 5),
    (:if, (:iseven, :k),
        (:println, "even: ", :k),
        (:println, "odd:  ", :k))) |> t2j

# %%
@tupp (:for, :k, (:(:), 1, 5),
    (:if, (:iseven, :k),
        (:println, "even: ", :k),
        (:println, "odd:  ", :k)))

# %%
(:block,
    (:(=), :U, (:lambda, (:u,), (:u, :u))),
    (:(=), :F, (:lambda, (:u,), (:lambda, (:n, :a, :b),
                (:if, (:(==), :n, 0), :a,
                    ((:u, :u), (:-, :n, 1), :b, (:+, :a, :b)))))),
    (:(=), :Fib, (:U, :F)),
    (:for, :k, (:(:), 1, 10),
        (:println, "Fib(", :k, ") = ", (:Fib, :k, 0, 1)))) |> t2j

# %%
@tupp (:block,
    (:(=), :U, (:lambda, (:u,), (:u, :u))),
    (:(=), :F, (:lambda, (:u,), (:lambda, (:n, :a, :b),
                (:if, (:(==), :n, 0), :a,
                    ((:u, :u), (:-, :n, 1), :b, (:+, :a, :b)))))),
    (:(=), :Fib, (:U, :F)),
    (:for, :k, (:(:), 1, 10),
        (:println, "Fib(", :k, ") = ", (:Fib, :k, 0, 1))))

# %%
(:block,
    (:using, :Plots),
    (:(=), :n, 20),
    (:(=), :x, (:range, 0, 2, (:kw, :length, :n))),
    (:(=), :y, (:+, (:., :sinpi, (:tuple, :x)), (:*, 0.2, (:randn, :n)))),
    (:(=), :xs, (:range, 0, 2, (:kw, :length, 200))),
    (:(=), :X, (:.^, :x, (:transpose, (:(:), 0, 3)))),
    (:(=), :b, (:\, :X, :y)),
    (:scatter, :x, :y, (:kw, :label, "sample")),
    (:plot!, :xs, (:., :sinpi, (:tuple, :xs)),
        (:kw, :label, "sinpi(x)"),
        (:kw, :color, (:quote, :black)),
        (:kw, :ls, (:quote, :dash))),
    (:plot!, :xs, (:., :evalpoly, (:tuple, :xs, (:Ref, :b))),
        (:kw, :label, "degree-3 polynomial"),
        (:kw, :color, 2), (:kw, :lw, 2))) |> t2j

# %%
@tupp (:block,
    (:using, :Plots),
    (:(=), :n, 20),
    (:(=), :x, (:range, 0, 2, (:kw, :length, :n))),
    (:(=), :y, (:+, (:., :sinpi, (:tuple, :x)), (:*, 0.2, (:randn, :n)))),
    (:(=), :xs, (:range, 0, 2, (:kw, :length, 200))),
    (:(=), :X, (:.^, :x, (:transpose, (:(:), 0, 3)))),
    (:(=), :b, (:\, :X, :y)),
    (:scatter, :x, :y, (:kw, :label, "sample")),
    (:plot!, :xs, (:., :sinpi, (:tuple, :xs)),
        (:kw, :label, "sinpi(x)"),
        (:kw, :color, (:quote, :black)),
        (:kw, :ls, (:quote, :dash))),
    (:plot!, :xs, (:., :evalpoly, (:tuple, :xs, (:Ref, :b))),
        (:kw, :label, "degree-3 polynomial"),
        (:kw, :color, 2), (:kw, :lw, 2)))

# %%
