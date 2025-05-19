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
# # ColabでJuliaを使うためのノートブック v0.3
#
# このノートブックの内容については再配布・改変・部分的コピーその他すべてを自由に行って構いません。
#
# このノートブックは[Google Colabで実行できる](https://colab.research.google.com/github/genkuroki/public/blob/main/0054/Julia%20notebook%20for%20Google%20Colab%20v0.3.ipynb).

# %%
# Google Colabと自分のパソコンの両方で使えるようにするための工夫

import Pkg

"""すでにPkg.add済みのパッケージのリスト (高速化のために用意)"""
_packages_added = [info.name for (uuid, info) in Pkg.dependencies() if info.is_direct_dep]

"""_packages_added内にないパッケージをPkg.addする"""
add_pkg_if_not_added_yet(pkg) = if !(pkg in _packages_added)
    println(stderr, "# $(pkg).jl is not added yet, so let's add it.")
    Pkg.add(pkg)
end

"""ASTからusing内の`.`を含まないモジュール名を抽出"""
function find_using_pkgs(ast::Expr)
    pkgs = String[]
    function traverse(expr::Expr)
        if expr.head == :using
            for arg in expr.args
                if arg.head == :. && length(arg.args) == 1
                    push!(pkgs, string(arg.args[1]))
                elseif arg.head == :(:) && length(arg.args[1].args) == 1
                    push!(pkgs, string(arg.args[1].args[1]))
                end
            end
        else
            for arg in expr.args arg isa Expr && traverse(arg) end
        end
    end
    traverse(ast)
    pkgs
end

"""必要そうなPkg.addを追加するマクロ"""
macro autoadd(x)
    pkgs = find_using_pkgs(x)
    :(add_pkg_if_not_added_yet.($(pkgs)); $x)
end

@autoadd begin
using Distributions
using StatsPlots
using QuadGK: quadgk
using Roots: find_zero, find_zeros
using SymPy: SymPy, sympy, @syms, oo
end

# %%
(@macroexpand @autoadd using A) |> Base.remove_linenums!

# %%
(@macroexpand @autoadd using A, B, C) |> Base.remove_linenums!

# %%
(@macroexpand @autoadd using A: a1, a2, @a3) |> Base.remove_linenums!

# %%
(@macroexpand @autoadd begin
using A: a1
using A.B
using A.C: c1, c2
#using D
using E, A.F, G
using H: h1, h2
using I
end) |> Base.remove_linenums!

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
