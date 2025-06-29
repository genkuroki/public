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
using StatsPlots
default(fmt=:png)
end

# %%
function wirehistogram3d(x::AbstractArray{T,1}, y::AbstractArray{T,1}; xlims::Tuple{Number,Number}=(0,0), ylims::Tuple{Number,Number}=(0,0), bins::Int64, kws...) where {T}

    # input check, if not given take extrema
    if all(xlims == (0,0))
        xlims = extrema(x)
    end
    if all(ylims == (0,0))
        ylims = extrema(y)
    end

    # create the x-y-grid to check values
    xgrid = range(xlims[1],xlims[end],length=bins+1)
    ygrid = range(ylims[1],ylims[end],length=bins+1)

    # count all x-y-values within a specific x-y-grid-cell (histcount may be nicer)
    z = [sum(all(((x[i],y[i]).>=(xgrid[itx],ygrid[ity])) .& ((x[i],y[i]).<(xgrid[itx+1],ygrid[ity+1]))) for i in eachindex(x)) for ity = 1:bins, itx = 1:bins]

    # Helpers for nice bars
    xdat = repeateps(xgrid,inner=(2));
    ydat = repeateps(ygrid,inner=(2));

    zdat = zeros(length(xdat),length(ydat));
    zdat[2:end-1,2:end-1] .= repeat(z,inner=(2,2));

    wireframe(xdat, ydat, zdat, xlims=xlims, ylims=ylims; kws...)
end

function histogram3d(x::AbstractArray{T,1}, y::AbstractArray{T,1}; xlims::Tuple{Number,Number}=(0,0), ylims::Tuple{Number,Number}=(0,0), bins::Int64, kws...) where {T}

    # input check, if not given take extrema
    if all(xlims == (0,0))
        xlims = extrema(x)
    end
    if all(ylims == (0,0))
        ylims = extrema(y)
    end

    # create the x-y-grid to check values
    xgrid = range(xlims[1],xlims[end],length=bins+1)
    ygrid = range(ylims[1],ylims[end],length=bins+1)

    # count all x-y-values within a specific x-y-grid-cell (histcount may be nicer)
    z = [sum(all(((x[i],y[i]).>=(xgrid[itx],ygrid[ity])) .& ((x[i],y[i]).<(xgrid[itx+1],ygrid[ity+1]))) for i in eachindex(x)) for ity = 1:bins, itx = 1:bins]

    # Helpers for nice bars
    xdat = repeateps(xgrid,inner=(2));
    ydat = repeateps(ygrid,inner=(2));

    zdat = zeros(length(xdat),length(ydat));
    zdat[2:end-1,2:end-1] .= repeat(z,inner=(2,2));
    
    @show size(xdat), size(ydat), size(zdat)

    surface(xdat, ydat, zdat, xlims=xlims, ylims=ylims; fill_z=[x for y in ydat, x in xdat], kws...)
    #surface(xdat, ydat, zdat, xlims=xlims, ylims=ylims; kws...)
end

repeateps(x;kws...) = repeat(x;kws...) .- [k%2 == 1 ? 1e-15 : 0 for k = 1:2*length(x)];

pyplot()

# %%
n = 1000
X, Y = randn(n), randn(n)
histogram3d(X, Y; bins=20, colorbar=false, c=:rainbow)

# %%
histogram3d(X, Y; xlim=(-3.1, 3.1), ylim=(-3.1, 3.1), bins=20, colorbar=false, camera=(30, 45), size=(500, 500), alpha=0.7, c=:rainbow)

# %%
pyplot()
x, y = range(-1, 1, 21), range(-1, 1, 21)
f(x, y) = (x, y) == (0, 0) ? zero(x*y) : x^2*y/(x^2 + y^2)
surface(x, y, f.(x', y); fill_z = f.(x', y) .+ 0.5, c=:rainbow)

# %%
pyplot()
x, y = range(-1, 1, 21), range(-1, 1, 21)
f(x, y) = (x, y) == (0, 0) ? zero(x*y) : x^2*y/(x^2 + y^2)
surface(x, y, f.(x', y); fill_z = [x for y in y, x in x], colorbar=false)
plot!(xguide="x", yguide="y")

# %%
f.(x', y) .+ 0.5

# %%
