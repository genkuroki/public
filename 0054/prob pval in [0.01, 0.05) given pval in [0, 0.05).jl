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
# https://x.com/enodon/status/1931248617672847786

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
@autoadd begin
using Distributions
using QuadGK
using Roots
using Plots
default(fmt=:png)
end

# %%
function pvalue_bin_score(k, n, p)
    phat = k/n
    sehat = sqrt(phat * (1 - phat) / n)
    z = (phat - p) / sehat
    2ccdf(Normal(), abs(z))
end

function expectval(f, bin::Binomial)
    sum(f(k) * pdf(bin, k) for k in support(bin))
end

function power_bin_score(n, p0, p1; alphamin=0.0, alphamax=0.05)
    expectval(k -> alphamin ≤ pvalue_bin_score(k, n, p0) < alphamax, Binomial(n, p1))
end

function find_p1(n, p0; power=0.8)
    f(p1) = power_bin_score(n, p0, p1) - power
    find_zero(f, (p0, 1.0))
end

# %%
@show n, p0 = 100, 0.3
@show p1 = find_p1(n, p0)
power_bin_score(n, p0, p0), power_bin_score(n, p0, p1), power_bin_score(n, p0, p1; alphamax=0.005)/power_bin_score(n, p0, p1), power_bin_score(n, p0, p1; alphamin=0.01)/power_bin_score(n, p0, p1)

# %%
@show n, p0 = 1000, 0.3
@show p1 = find_p1(n, p0)
power_bin_score(n, p0, p0), power_bin_score(n, p0, p1), power_bin_score(n, p0, p1; alphamax=0.005)/power_bin_score(n, p0, p1), power_bin_score(n, p0, p1; alphamin=0.01)/power_bin_score(n, p0, p1)

# %%
@show n, p0 = 10000, 0.3
@show p1 = find_p1(n, p0)
power_bin_score(n, p0, p0), power_bin_score(n, p0, p1), power_bin_score(n, p0, p1; alphamax=0.005)/power_bin_score(n, p0, p1), power_bin_score(n, p0, p1; alphamin=0.01)/power_bin_score(n, p0, p1)

# %%
@show n, p0 = 10000, 0.03
@show p1 = find_p1(n, p0)
power_bin_score(n, p0, p0), power_bin_score(n, p0, p1), power_bin_score(n, p0, p1; alphamax=0.005)/power_bin_score(n, p0, p1), power_bin_score(n, p0, p1; alphamin=0.01)/power_bin_score(n, p0, p1)

# %%
@show n, p0 = 10000, 0.003
@show p1 = find_p1(n, p0)
power_bin_score(n, p0, p0), power_bin_score(n, p0, p1), power_bin_score(n, p0, p1; alphamax=0.005)/power_bin_score(n, p0, p1), power_bin_score(n, p0, p1; alphamin=0.01)/power_bin_score(n, p0, p1)

# %%
@show n, p0 = 10000, 0.0003
@show p1 = find_p1(n, p0)
power_bin_score(n, p0, p0), power_bin_score(n, p0, p1), power_bin_score(n, p0, p1; alphamax=0.005)/power_bin_score(n, p0, p1), power_bin_score(n, p0, p1; alphamin=0.01)/power_bin_score(n, p0, p1)

# %%
function pvalue_ztest(xbar, n)
    2ccdf(Normal(0, 1/sqrt(n)), abs(xbar))
end

function power_ztest(n, mu1; alphamin=0.0, alphamax=0.05)
    a, b, c, d = quantile.(Normal(0, 1/sqrt(n)), (alphamin/2, alphamax/2, 1-alphamax/2, 1-alphamin/2))
    normal1 = Normal(mu1, 1/sqrt(n))
    cdf(normal1, b) - cdf(normal1, a) + ccdf(normal1, c) - ccdf(normal1, d)
end

function find_mu1(n; power=0.8)
    f(mu1) = power_ztest(n, mu1) - power
    find_zero(f, (0.0, 10.0))
end

# %%
@show n = 100
@show mu1 = find_mu1(n)
power_ztest(n, 0.0), power_ztest(n, mu1), power_ztest(n, mu1; alphamax=0.005)/power_ztest(n, mu1), power_ztest(n, mu1; alphamin=0.01)/power_ztest(n, mu1)

# %%
@show n = 1000
@show mu1 = find_mu1(n)
power_ztest(n, 0.0), power_ztest(n, mu1), power_ztest(n, mu1; alphamax=0.005)/power_ztest(n, mu1), power_ztest(n, mu1; alphamin=0.01)/power_ztest(n, mu1)

# %%
@show n = 10000
@show mu1 = find_mu1(n)
power_ztest(n, 0.0), power_ztest(n, mu1), power_ztest(n, mu1; alphamax=0.005)/power_ztest(n, mu1), power_ztest(n, mu1; alphamin=0.01)/power_ztest(n, mu1)

# %%
