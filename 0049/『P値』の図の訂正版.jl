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
#     display_name: Julia 1.10.4
#     language: julia
#     name: julia-1.10
# ---

# %%
using Distributions
using StatsPlots
default(fmt=:png)

# %%
plot(Normal(); label="", c=:black)
annotate!((1.6, 0.2, text("ϕ(x)", 16)))
plot!([-0.98, -0.02], fill(pdf(Normal(), -1), 2); label="", arrow=:both, c=:black)
annotate!((-0.5, 0.26, text("σ", 16)))
plot!([-1.98, -0.02], fill(pdf(Normal(), -2), 2); label="", arrow=:both, c=:black)
annotate!((-1, 0.07, text("2σ", 16)))
vline!([0]; label="", c=:black, ls=:dot)
plot!(xtick=(-2:2, ["μ-2σ", "μ-σ", "μ", "μ+σ", "μ+2σ"]), ytick=false, tickfontsize=16)
plot!(yaxis=false, ylim=(0, 0.42))

# %%
plot(Normal(0, 10), -25, 25; label="", c=:black, ls=:dash)
annotate!((10, 0.05, text("n=1", 16)))
plot!(Normal(); label="", c=:black)
annotate!((1.8, 0.2, text("n=100", 16, :left)))
plot!([-1.98, -0.02], fill(pdf(Normal(), -2), 2); label="", arrow=:both, c=:black)
plot!([-12, -5, -1], [0.15, 0.13, 0.06]; label="", arrow=true, c=:black)
annotate!((-12.2, 0.155, text("2σ", 16, :right)))
vline!([0]; label="", c=:black, ls=:dot)
plot!(xtick=([0], ["μ"]), ytick=false, tickfontsize=16)
plot!(yaxis=false, ylim=(0, 0.42))

# %%
plot(TDist(10), -4, 4; label="", c=:black)
annotate!((-1.0, 0.27, text("pdf of T", 16, :right)))
vline!([0]; label="", c=:black, ls=:dot)
plot!(fill(-1.8, 2), [0, pdf(TDist(10), -1.8)]; label="", c=:black)
plot!(TDist(10), -4, -1.8; label="", fill=true, c=:black, fillalpha=0.2)
plot!(fill(1.8, 2), [0, pdf(TDist(10), 1.8)]; label="", c=:black)
plot!(TDist(10), 1.8, 4; label="", fill=true, c=:black, fillalpha=0.2)
plot!(xtick=([-1.8, 0, 1.8], ["-|A|", "0", "|A|"]), ytick=false, tickfontsize=16)
plot!(yaxis=false, ylim=(0, 0.42))
plot!([(1.7, 0.3), (-2.1, 0.13), (-2.1, 0.03)]; label="", arrow=true, c=:black, ls=:dot)
plot!([(2.2, 0.28), (2.1, 0.03)]; label="", arrow=true, c=:black, ls=:dot)
annotate!((1.8, 0.3, text("P-value of Δ=0", 16, :left)))

# %%
