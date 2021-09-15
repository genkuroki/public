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
_foo(expr) = expr
_foo(expr::Symbol) = :foo
function _foo(expr::Expr)
    if expr.head === :macrocall
        Expr(expr.head, expr.args[1], _foo.(expr.args[2:end])...)
    else
        Expr(expr.head, _foo.(expr.args)...)
    end
end

macro foo(expr)
    _foo(expr) |> esc
end

# %%
_bar(expr) = expr
_bar(expr::Symbol) = :bar
function _bar(expr::Expr)
    if expr.head === :macrocall
        Expr(expr.head, expr.args[1], _bar.(expr.args[2:end])...)
    else
        Expr(expr.head, _bar.(expr.args)...)
    end
end

macro bar(expr)
    _bar(expr) |> esc
end

# %%
@macroexpand @foo hoge(x)

# %%
@macroexpand @bar hoge(x)

# %%
@macroexpand @bar @foo hoge(x)

# %%
@macroexpand hoge(@foo @bar x)

# %%
foo(x) = "foo"

# %%
bar(x) = "bar"

# %%
@foo hoge(x)

# %%
@bar hoge(x)

# %%
@bar @foo hoge(x)

# %%
@foo @bar hoge(x)

# %%
@macroexpand1 @foo @bar hoge(x)

# %%
@macroexpand @foo @bar hoge(x)

# %%
