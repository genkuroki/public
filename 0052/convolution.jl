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
#     display_name: Julia 1.11.2
#     language: julia
#     name: julia-1.11
# ---

# %% [markdown]
# https://x.com/dannchu/status/1880296123061137818

# %%
using QuadGK,Plots
default(fmt=:png)

f_X(x) = 0 ≤ x ≤ 1 ? 1 : 0
f_Y(y) = 0 ≤ y ≤ 1 ? 1 : 0

function f_Z(z)
    quadgk(t -> f_X(z-t)*f_Y(t) , -Inf, Inf)[1]
end

plot(f_Z,label="Z=X+Y")

# %% [markdown]
# ## 不連続関数の`quadgk`はうまく行かない.

# %%
plot(range(-1, 3, 1001), f_Z; label="Z=X+Y")

# %%
@show f_Z(1.89) f_Z(1.9);

# %% [markdown]
# ## convolutionの実装
#
# <img src="IMG_7804.jpeg" width=500>

# %%
using Distributions
using QuadGK
using StatsPlots
default(fmt=:png, legendfontsize=12)

# Assume that f is continuous on (a, b) and g is continuous on (c, d). 
function _convolution((f, a, b), (g, c, d))
    function h(z)
        (z ≤ a+c || b+d ≤ z) && return zero(float(z))
        if b+c ≤ a+d
            tmin = z ≤ a+d ? a : z-d
            tmax = z ≤ b+c ? z-c : b
            quadgk(t->f(t)*g(z-t), tmin, tmax)[1]
        else
            tmin = z ≤ b+c ? c : z-b
            tmax = z ≤ a+d ? z-a : d
            quadgk(t->f(z-t)*g(t), tmin, tmax)[1]
        end
    end
    return (h, a+c, b+d)
end

pdf_convolution(distx::ContinuousUnivariateDistribution, disty::ContinuousUnivariateDistribution) =
    _convolution((x->pdf(distx, x), extrema(distx)...), (y->pdf(disty, y), extrema(disty)...))

pdf_convolution(distx::ContinuousUnivariateDistribution, (g, c, d)) =
    _convolution((x->pdf(distx, x), extrema(distx)...), (g, c, d))

pdf_convolution((f, a, b), disty::ContinuousUnivariateDistribution) =
    _convolution((f, a, b), (y->pdf(disty, y), extrema(disty)...))

# %%
d = Uniform()
f2 = pdf_convolution(d, d)
f3 = pdf_convolution(d, f2)
f4 = pdf_convolution(d, f3)
f5 = pdf_convolution(d, f4)
@time f5[1](1.0)

# %%
@time plot(f2[1], f2[2]-1, f2[3]+1; label="X₁+X₂")
plot!(Normal(2mean(Uniform()), √2 * std(Uniform())); label="normal approx.", ls=:dash)
title!("Xᵢ ~ Uniform(0, 1), i.i.d.")

# %%
@time plot(f3[1], f3[2]-1, f3[3]+1; label="X₁+X₂+X₃")
plot!(Normal(3mean(Uniform()), √3 * std(Uniform())); label="normal approx.", ls=:dash)
title!("Xᵢ ~ Uniform(0, 1), i.i.d.")

# %%
@time plot(f4[1], f4[2]-1, f4[3]+1; label="X₁+X₂+X₃+X₄")
plot!(Normal(4mean(Uniform()), √4 * std(Uniform())); label="normal approx.", ls=:dash)
title!("Xᵢ ~ Uniform(0, 1), i.i.d.")

# %%
@time plot(f5[1], f5[2]-1, f5[3]+1; label="X₁+X₂+⋯+X₅")
plot!(Normal(5mean(Uniform()), √5 * std(Uniform())); label="normal approx.", ls=:dash)
title!("Xᵢ ~ Uniform(0, 1), i.i.d.")

# %% [markdown]
# 本当は[`Distributions.convolve`](https://juliastats.org/Distributions.jl/latest/convolution/)を任意の単変量連続分布達について実装するべき.

# %%
distx, disty = Normal(1, 3), Normal(2, 4)
@show distx disty
h = pdf_convolution(distx, disty)
distz = convolve(distx, disty)
plot(h[1], -25, 30; label="X+Y (quadgk)")
plot!(distz, -25, 30; label="X+Y (convolve)", ls=:dash)
title!("X ~ Normal(1, 3), Y ~ Normal(2, 4), i.i.d.")

# %%
distx, disty = Gamma(2, 1), Gamma(3, 1)
@show distx disty
h = pdf_convolution(distx, disty)
distz = convolve(distx, disty)
plot(h[1], -1, 18; label="X+Y (quadgk)")
plot!(distz, -1, 18; label="X+Y (convolve)", ls=:dash)
title!("X ~ Gamma(2, 1), Y ~ Gamma(3, 1)")

# %%
distx, disty = Normal(0, 1), Gamma(3, 1)
@show distx disty
h = pdf_convolution(distx, disty)
distz = convolve(distx, disty)
plot(h[1], -4, 12; label="X+Y (quadgk)")
plot!(distz, -4, 18; label="X+Y (convolve)", ls=:dash)
title!("X ~ Normal(0, 1), Y ~ Gamma(3, 1), i.i.d.")

# %%
distx, disty = Normal(0, 1), Gamma(3, 1)
@show distx disty
h = pdf_convolution(distx, disty)
#distz = convolve(distx, disty)
plot(h[1], -4, 12; label="X+Y (quadgk)")
#plot!(distz, -4, 18; label="X+Y (convolve)", ls=:dash)
title!("X ~ Normal(0, 1), Y ~ Gamma(3, 1), i.i.d.")

# %%
