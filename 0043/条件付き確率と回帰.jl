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
#     display_name: Julia 1.9.1
#     language: julia
#     name: julia-1.9
# ---

# %%
using Distributions
using StatsBase
using StatsPlots
default(fmt=:png)

# %%
# X の値が xbin[i] 以上 xbin[i+1] 以下になる部分における
# Y の値が ybin[j] 以上 ybin[j+1] 以下になる部分の割合を
# P[i,j] とする。

function condprob(X, Y, xbin, ybin)
    m, n = length(xbin)-1, length(ybin)-1
    P = Matrix{Float64}(undef, m, n)
    for i in 1:m
        xidx = @.(xbin[i] ≤ X ≤ xbin[i+1])
        countx = count(xidx)
        for j in 1:n
            yidx = @.(ybin[j] ≤ Y[xidx] ≤ ybin[j+1])
            P[i,j] = count(yidx)/countx
        end
    end
    P
end

# %%
# Xの分布

distx = MixtureModel([Normal(), Normal(6, 2)], [1/2, 1/2])
plot(x -> pdf(distx, x), -4, 14; label="distx")

# %%
# モデルを y = 1 + 0.5x + rand(distu) とする。

distu = Normal(0, 2)
β = [1, 0.5]

# %%
# サンプル生成

n = 10000
X = rand(distx, n)
Y = evalpoly.(X, (β,)) + rand(distu, n)
scatter(X, Y; label="sample", msc=:auto, alpha=0.2, ms=2)
plot!(xguide="x", yguide="y")

# %%
# 最小二乗法

A = X .^ (0:1)'
@show β̂ = A \ Y

scatter(X, Y; label="sample", msc=:auto, alpha=0.3, ms=2)
plot!(x -> evalpoly(x, β̂); label="regression line", lw=3)
plot!(xguide="x", yguide="y")

# %%
# メッシュで区切る

xbin = -2:0.5:10
ybin = -4:0.5:10

scatter(X, Y; label="sample", msc=:auto, alpha=0.3, ms=2)
plot!(x -> evalpoly(x, β̂); label="regression line", lw=3)
plot!(xguide="x", yguide="y")
vline!(xbin; label="")
hline!(ybin; label="")
plot!(xlim=extrema(xbin), ylim=extrema(ybin))

# %%
# P[i,j] = (格子長方形の内側の点の個数) / (その講師超法権を含む縦の帯の中の点の個数)

P = condprob(X, Y, xbin, ybin)
heatmap(xbin[1:end-1] .+ 0.25, ybin[1:end-1] .+ 0.25, P'; colorbar=false)
plot!(xlim=extrema(xbin), ylim=extrema(ybin))
plot!(x -> evalpoly(x, β̂); label="regression line", c=:cyan, lw=3)
plot!(xguide="x", yguide="y")

# %%
