# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:hydrogen
#     text_representation:
#       extension: .jl
#       format_name: hydrogen
#       format_version: '1.3'
#       jupytext_version: 1.10.3
#   kernelspec:
#     display_name: Julia current stable release
#     language: julia
#     name: julia
# ---

# %% [markdown]
# # ColabでJuliaを使うためのノートブック v0.2
#
# このノートブックの内容は再配布や改変や部分的コピーその他すべて自由に行って構いません。
#
# このノートブックは[Google Colabで実行できる](https://colab.research.google.com/github/genkuroki/public/blob/main/0054/Julia%20notebook%20for%20Google%20Colab%20v0.2.ipynb).

# %%
# Google Colabと自分のパソコンの両方で使えるようにするための工夫

import Pkg

"""すでにPkg.add済みのパッケージのリスト (高速化のために用意)"""
_packages_added = [info.name for (uuid, info) in Pkg.dependencies() if info.is_direct_dep]

"""
    @_using A
    @_using A: a₁, ..., aₙ
必要ならPkg.add("A")を実行してからusingを実行してくれる.
ただし, aᵢとして@fooを採用したい場合には"@foo"と書く.
"""
macro _using(x)
    modsymb = x isa Symbol ? x : x.head == :call ? x.args[2] : x.args[1].args[2]
    pkg = string(modsymb)
    if !(pkg in _packages_added)
        println(stderr, "# $(pkg).jl is not added yet, so let's add it.")
        Pkg.add(pkg)
    end
    if x isa Symbol
        Expr(:using, Expr(:., modsymb))
    elseif x.head == :call
        Expr(:using, Expr(:(:), Expr(:., modsymb), Expr(:., Symbol(x.args[3]))))
    else
        expr = Expr(:using, Expr(:(:), Expr(:., modsymb), Expr(:., Symbol(x.args[1].args[3]))))
        for i in 2:length(x.args) push!(expr.args[1].args, Expr(:., Symbol(x.args[i]))) end
        expr
    end
end

(@doc @_using) |> display

@_using Distributions
@_using StatsPlots
@_using QuadGK: quadgk
@_using Roots: find_zero, find_zeros
@_using SymPy: SymPy, sympy, "@syms", oo

# %%
dist = Gamma(2, 1/2)
plot(dist; label="$dist")

# %%
quadgk(x -> pdf(dist, x), 0, Inf)

# %%
f(x) = pdf(dist, x) - pdf(dist, 0.1)
find_zero(f, (0, 0.2))

# %%
find_zeros(f, (0, 10))

# %%
@syms x::real, a::positive

# %%
sympy.integrate(exp(-x^2/a), (x, -oo, oo))

# %%
