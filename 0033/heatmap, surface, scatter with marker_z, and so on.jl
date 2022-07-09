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
#     display_name: Julia 1.7.3
#     language: julia
#     name: julia-1.7
# ---

# %%
using Plots
f(x, y) = cos(x*y)
x = range(-4, 4, 321)
y = range(-2.5, 2.5, 201)
heatmap(x, y, f; c=:bwr)

# %%
surface(x, y, f; camera=(30, 75), c=:bwr)

# %%
z = f.(x', y)
heatmap(x, y, z; c=:bwr)

# %%
surface(x, y, z; camera=(30, 75), c=:bwr)

# %%
using Plots

f(x, y) = cos(x*y)
x = range(-4, 4, 81)
y = range(-2.5, 2.5, 51)
xx, yy = vec.((z -> (real(z), imag(z)))(complex.(x', y)))
zz = f.(xx, yy)
xx, yy, zz

# %%
scatter(xx, yy; marker_z=zz, c=:bwr, ms=2, msw=0, label="")

# %%
using Plots
using Distributions

f(x, y) = cos(x*y)
n = 2^13
xx = rand(Normal(0, 4), n)
yy = rand(Normal(0, 4), n)
zz = f.(xx, yy)
xx, yy, zz

scatter(xx, yy; marker_z=zz, c=:bwr, ms=3, msw=0, label="", alpha=0.5)
plot!(xlim=(-4, 4), ylim=(-2.5, 2.5))

# %%
using Plots

f(x, y) = cos(x*y)
x = range(-4, 4, 81)
y = range(-2.5, 2.5, 51)
xx, yy = vec.((z -> (real(z), imag(z)))(complex.(x', y)))
zz = f.(xx, yy)

X = reshape(xx, 51, 81)[1,:]
Y = reshape(yy, 51, 81)[:,1]
Z = reshape(zz, 51, 81)
heatmap(X, Y, Z; c=:bwr)

# %%
using Plots
default(titlefontsize=10)
using Distributions

function plot_nbgam(α, θ, L)
    p = 1/(L*θ)
    @show mean((α + NegativeBinomial(α, p))/L), mean(Gamma(α, θ))
    @show std((α + NegativeBinomial(α, p))/L), std(Gamma(α, θ))
    t = range(max(-1/(2L), α*θ - 4√α*θ), α*θ + 6√α*θ, 1000)
    plot(t, t -> pdf(NegativeBinomial(α, p), round(Int, L*t - α))*L; label="(α+NB(α,p))/L")
    plot!(t, t -> pdf(Gamma(α, θ), t); label="Gamma(α, θ)", ls=:dash)
    title!("α=$α, θ=$θ, L=$L, p=1/(Lθ)")
end

# %%
plot(plot_nbgam.(1, 1, (1.2, 3, 10))...; size=(800, 200), layout=(1, 3))

# %%
plot(plot_nbgam.(2, 2, (1.2, 3, 10)./2)...; size=(800, 200), layout=(1, 3))

# %%
plot(plot_nbgam.(5, 4, (1.2, 3, 10)./4)...; size=(800, 200), layout=(1, 3))

# %%
using Plots
default(titlefontsize=10)
using Distributions

function plot_bbnb(L, α, θ; kwargs...)
    bb = BetaBinomial(L, α, L/θ)
    nb = NegativeBinomial(α, 1/(1+θ))
    μ, σ = mean(nb), std(nb)
    m = range(max(-1, μ-4σ), μ+4σ, 1000)
    plot(m, m -> pdf(bb, round(m)); label="BetaBin(L,α,L/θ)")
    plot!(m, m -> pdf(nb, round(m)); label="NegBin(α,1/(1+θ))", ls=:dash)
    title!("L=$L, α=$α, θ=$θ, p=1/(1+θ)=1/$(1+θ)")
    plot!(; kwargs...)
end

# %%
plot(plot_bbnb.((10, 30), 3, 2)...; size=(800, 250))

# %%
plot(plot_bbnb.((100, 300), 3, 2)...; size=(800, 250))

# %%
