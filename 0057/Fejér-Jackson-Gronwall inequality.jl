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
#     display_name: Julia
#     language: julia
#     name: julia
# ---

# %% [markdown]
# # Fejér-Jackson-Gronwall 不等式
#
# * 黒木玄
# * 2026-03-28
#
# Google Colabで実行:
#
# * https://colab.research.google.com/github/genkuroki/public/blob/main/0057/Fejér-Jackson-Gronwall%20inequality.ipynb
#
# 参照資料:
#
# * https://math.stackexchange.com/questions/376273/inequality-sum-limits-1-le-k-le-n-frac-sin-kxk-ge-0-fejer-jackson
#
# 以下のグラフを見れば, これに書かれている帰納法による証明をすぐに理解できる.

# %%
using Plots
default(fmt=:png, legendfontsize=10)
f(n, x) = sum(sin(k*x)/k for k in 1:n)

PP = []
for n in 2:11
    P = plot(x -> f(n, x), 0, pi; label="\$f_{$n}(x)\$")
    scatter!(0:2pi/n:pi, x -> f(n, x); msc=:auto, label="local minima")
    #scatter!(pi/(n+1):2pi/(n+1):pi, x -> f(n, x); msc=:auto, label="local maxima")
    plot!(x -> f(n-1, x), 0, pi; ls=:dot, label="\$f_{$(n-1)}(x)\$")
    plot!(ylim=(-0.05, 1.7))
    push!(PP, P)
end

plot(PP...; size=(1000, 1500), layout=(5, 2))
plot!(plot_title=raw"Fejér-Jackson-Gronwall inequality:  $f_n(x) = \sum_{k=1}^n \frac{\sin(kx)}{k} > 0 \quad (0<x<\pi)$")

# %%
