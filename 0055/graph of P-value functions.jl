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

# %%
using Distributions
using StatsPlots
default(fmt=:png, legend_font_halign=:left)
pgfplotsx()
PGFPlotsX.CUSTOM_PREAMBLE = [raw"\usepackage{luatexja}"]
@eval Plots pgfx_sanitize_string(s::AbstractString) = s

function pvalue_central(k, n, p)
    bin = Binomial(n, p)
    min(1, 2cdf(bin, k), 2ccdf(bin, k-1))
end

# %%
gr()
const mincho = "ipamp"
const gothic = "ipagp"

plot()
plot!(p -> pvalue_central(15, 20, p), 0, 1; label="20回中15回当たり (先行報告)", ls=:dash)
plot!(p -> pvalue_central(19, 30, p), 0, 1; label="30回中19回当たり (本追試)")
plot!(p -> pvalue_central(34, 50, p), 0, 1; label="50回中34回当たり (統合データ)", ls=:dot)
vline!([0.5]; label="帰無値 1/2", c=:black, lw=0.5, ls=:dot)
plot!(xtick=0:0.1:1, ytick=0:0.05:1.5)
plot!(xguide="表が出る確率の仮説値", yguide="P値 = 仮説値とデータの値の相性の良さ")
title!("P値関数達のグラフ")
#plot!(legend=:topleft)
plot!(titlefontsize=16, guidefontsize=12, legendfontsize=11)
plot!(fontfamily=gothic, size=(720, 480), leftmargin=6Plots.mm)

# %%
pgfplotsx()

plot()
plot!(p -> pvalue_central(15, 20, p), 0, 1; label="20回中15回当たり (先行報告)", ls=:dash)
plot!(p -> pvalue_central(19, 30, p), 0, 1; label="30回中19回当たり (本追試)")
plot!(p -> pvalue_central(34, 50, p), 0, 1; label="50回中34回当たり (統合データ)", ls=:dot)
vline!([0.5]; label="帰無値 1/2", c=:black, lw=0.5, ls=:dot)
plot!(xtick=0:0.1:1, ytick=0:0.05:1.5)
plot!(xguide="表が出る確率の仮説値", yguide="P値 = 仮説値とデータの値の相性の良さ")
title!("\\bfseries{P値関数達のグラフ}")
plot!(legend=:topleft)
plot!(titlefontsize=16, guidefontsize=12, legendfontsize=11)
#plot!(fontfamily=gothic, size=(750, 500), leftmargin=6Plots.mm)

# %%
