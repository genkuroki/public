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

# %% [markdown]
# https://qiita.com/ytaki0801/items/fc785de02eb45cde94fa

# %%
fib21 = 
[[["lambda", ["u"], ["u", "u"]],
        ["lambda", ["u"],
            ["lambda", ["n", "a", "b"],
                ["if", ["=", "n", 0], "a",
                    [["u", "u"], ["-", "n", 1], "b", ["+", "a", "b"]]]]]],
    21, 0, 1]

g = Dict("="=>"==","+"=>"+","-"=>"-")
function ev(s, e)
    if isa(s, String) return merge(g, e)[s]
    elseif isa(s, Int) return s
    elseif s[1] == "if" return ev(s[2], e) ? ev(s[3], e) : ev(s[4], e)
    elseif s[1] == "lambda" return vcat(s, [e])
    else
        f = ev(s[1], e); a = map(x -> ev(x, e), s[2:end])
        if isa(f, String) return eval(Expr(:call, Symbol(f), a...))
        else return ev(f[3], merge(f[4], Dict(k for k = zip(f[2], a))))
        end
    end
end

@show ev(fib21, Dict());

# %%
fib21 = 
[[["lambda", ["u"], ["u", "u"]],
        ["lambda", ["u"],
            ["lambda", ["n", "a", "b"],
                ["if", ["=", "n", 0], "a",
                    [["u", "u"], ["-", "n", 1], "b", ["+", "a", "b"]]]]]],
    21, 0, 1]

g = Dict("=" => "==", "+" => "+", "-" => "-")
function ev(s, e)
    s isa String     && return merge(g, e)[s]
    s isa Int        && return s
    s[1] == "if"     && return ev(s[2], e) ? ev(s[3], e) : ev(s[4], e)
    s[1] == "lambda" && return vcat(s, [e])
    f, a = ev(s[1], e), ev.(s[2:end], Ref(e))
    f isa String     && return eval(Expr(:call, Symbol(f), a...))
    return ev(f[3], merge(f[4], Dict(zip(f[2], a))))
end

@show ev(fib21, Dict());

# %%
fib21 = 
[[["lambda", ["u"], ["u", "u"]],
        ["lambda", ["u"],
            ["lambda", ["n", "a", "b"],
                ["if", ["=", "n", 0], "a",
                    [["u", "u"], ["-", "n", 1], "b", ["+", "a", "b"]]]]]],
    21, 0, 1]

g = Dict("=" => ==, "+" => +, "-" => -)
function ev(s, e)
    s isa String     && return merge(g, e)[s]
    s isa Int        && return s
    s[1] == "if"     && return ev(s[2], e) ? ev(s[3], e) : ev(s[4], e)
    s[1] == "lambda" && return vcat(s, [e])
    f, a = ev(s[1], e), ev.(s[2:end], Ref(e))
    f isa Function   && return f(a...)
    return ev(f[3], merge(f[4], Dict(zip(f[2], a))))
end

@show ev(fib21, Dict());

# %%
fib21 = 
[[[:lambda, [:u], [:u, :u]],
        [:lambda, [:u],
            [:lambda, [:n, :a, :b],
                [:if, [:(=), :n, 0], :a,
                    [[:u, :u], [:-, :n, 1], :b, [:+, :a, :b]]]]]],
    21, 0, 1]

g = Dict(:(=) => ==, :+ => +, :- => -)
function ev(s, e)
    s isa Symbol     && return merge(g, e)[s]
    s isa Int        && return s
    s[1] === :if     && return ev(s[2], e) ? ev(s[3], e) : ev(s[4], e)
    s[1] === :lambda && return vcat(s, [e])
    f, a = ev(s[1], e), ev.(s[2:end], Ref(e))
    f isa Function   && return f(a...)
    return ev(f[3], merge(f[4], Dict(zip(f[2], a))))
end

@show ev(fib21, Dict());

# %%
fib(n, a, b) = (u -> u(u))(u -> (n, a, b) -> n == 0 ? a : u(u)(n-1, b, a+b))(n, a, b)
@show fib(21, 0, 1)
@show [fib(n, 0, 1) for n in 1:10];

# %%
f(u, n, a, b) = n == 0 ? a : u(u, n-1, b, a+b)
@show f(f, 21, 0, 1)
@show [f(f, n, 0, 1) for n in 1:10];

# %%
fib21 = 
[[["lambda", ["u"], ["u", "u"]],
        ["lambda", ["u"],
            ["lambda", ["n", "a", "b"],
                ["if", ["=", "n", 0], "a",
                    [["u", "u"], ["-", "n", 1], "b", ["+", "a", "b"]]]]]],
    21, 0, 1]

g = Dict("=" => :(==), "+" => :+, "-" => :-)
function ev(s)
    s isa String     && return get(g, s, Symbol(s))
    s isa Int        && return s
    s[1] == "if"     && return Expr(:if, ev(s[2]), ev(s[3]), ev(s[4]))
    s[1] == "lambda" && return Expr(:(->), Expr(:tuple, Symbol.(s[2])...), ev(s[3]))
    Expr(:call, ev.(s)...)
end
macro ev(x) ev(Core.eval(__module__, x)) end

@show @ev fib21;

# %%
ev(fib21)

# %%
Fib = @ev [["lambda", ["u"], ["u", "u"]],
    ["lambda", ["u"],
        ["lambda", ["n", "a", "b"],
            ["if", ["=", "n", 0], "a",
                [["u", "u"], ["-", "n", 1], "b", ["+", "a", "b"]]]]]]

@show Fib(21, 0, 1)
@show [Fib(n, 0, 1) for n in 1:10];

# %%
@code_native debuginfo=:none Fib(21, 0, 1)

# %%
f(u, n, a, b) = n == 0 ? a : u(u, n-1, b, a+b)
@code_native debuginfo=:none f(f, 21, 0, 1)

# %%
fib21 = 
[[[:lambda, [:u], [:u, :u]],
        [:lambda, [:u],
            [:lambda, [:n, :a, :b],
                [:if, [:(==), :n, 0], :a,
                    [[:u, :u], [:-, :n, 1], :b, [:+, :a, :b]]]]]],
    21, 0, 1]

g = Dict(:(=) => :(==), :+ => :+, :- => :-)
function ev(s)
    s isa Symbol     && return get(g, s, s)
    s isa Int        && return s
    s[1] === :if     && return Expr(:if, ev(s[2]), ev(s[3]), ev(s[4]))
    s[1] === :lambda && return Expr(:(->), Expr(:tuple, s[2]...), ev(s[3]))
    Expr(:call, ev.(s)...)
end
macro ev(x) ev(Core.eval(__module__, x)) end

@show @ev fib21;

# %%
ev(fib21)

# %%
Fib = @ev [[:lambda, [:u], [:u, :u]],
        [:lambda, [:u],
            [:lambda, [:n, :a, :b],
                [:if, [:(==), :n, 0], :a,
                    [[:u, :u], [:-, :n, 1], :b, [:+, :a, :b]]]]]]

@show Fib(21, 0, 1)
@show [Fib(n, 0, 1) for n in 1:10];

# %%
@code_native debuginfo=:none Fib(21, 0, 1)

# %%
f(u, n, a, b) = n == 0 ? a : u(u, n-1, b, a+b)
@code_native debuginfo=:none f(f, 21, 0, 1)

# %%
