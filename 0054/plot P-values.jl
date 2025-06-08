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
# [Google Colabで実行]()

# %%
# Google Colabと自分のパソコンの両方で使えるようにするための工夫

using Pkg

"""すでにPkg.add済みのパッケージのリスト (高速化のために用意)"""
_packages_added = [sort!(readdir(Sys.STDLIB));
    [info.name for (uuid, info) in Pkg.dependencies() if info.is_direct_dep]]

"""_packages_added内にないパッケージをPkg.addする"""
add_pkg_if_not_added_yet(pkg) = if !(pkg in _packages_added)
    println(stderr, "# $(pkg).jl is not added yet, so let's add it.")
    Pkg.add(pkg)
end

"""expr::Exprからusing内の`.`を含まないモジュール名を抽出"""
function find_using_pkgs(expr::Expr)
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
    traverse(expr)
    pkgs
end

"""必要そうなPkg.addを追加するマクロ"""
macro autoadd(expr)
    pkgs = find_using_pkgs(expr)
    :(add_pkg_if_not_added_yet.($(pkgs)); $expr)
end

# %%
@autoadd using Distributions
@autoadd using Plots
default(fmt=:png, legendfontsize=12, guidefontsize=12, titlefontsize=12)

rd(x) = round(x; sigdigits=3)
safediv(x, y) = x == 0 ? zero(x/y) : x/y

function pvalue_bin_score(k, n, p)
    phat = k/n
    se = sqrt(p*(1-p)/n)
    z = safediv(phat - p, se)
    2ccdf(Normal(), abs(z))
end

function plot_pvalue_bin(; k=11, n=20, xann=0.49, yann=0, kwargs...)
    nullpval = pvalue_bin_score(k, n, 0.5)
    plot(p -> pvalue_bin_score(k, n, p), 0, 1; label="")
    vline!([0.5]; ls=:dot, lw=0.5, c=:black, label="")
    scatter!([0.5], [nullpval]; label="", msc=:auto, c=:black)
    annotate!(xann, nullpval+yann, text("P-value of p=0.5: $(rd(100nullpval))%", :right, :black, 12))
    plot!(xtick=0:0.1:1, ytick=0:0.05:1)
    plot!(xguide="p = c", yguide="P-value")
    title!("data: heads k=$k times in n=$n coin tosses")
    plot!(; kwargs...)
end

# %%
plot_pvalue_bin(; k=501, n=1000, xlim=(0.45, 0.55), xann=0.498, xtick=0:0.01:1)

# %%
plot_pvalue_bin(; k=5010, n=10000, xlim=(0.45, 0.55), xann=0.498, xtick=0:0.01:1)

# %%
plot_pvalue_bin(; k=50100, n=100000, xlim=(0.45, 0.55), xann=0.498, xtick=0:0.01:1)

# %%
plot_pvalue_bin(; k=501000, n=10^6, xlim=(0.45, 0.55), xann=0.498, xtick=0:0.01:1)

# %%
