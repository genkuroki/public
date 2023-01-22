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
#     display_name: Julia 1.9.0-beta3
#     language: julia
#     name: julia-1.9
# ---

# %% [markdown]
# * https://twitter.com/yudapearl/status/1616051192685277184
# * https://en.wikipedia.org/wiki/Lord%27s_paradox

# %%
using Distributions
using LinearAlgebra
dot2(x) = LinearAlgebra.dot(x, x)
using StatsPlots
default(fmt=:png)

# %%
a, b = 55, 70
A = B = [
    100 50
    50 100
]
distF = MvNormal([a, a], A)
distM = MvNormal([b, b], B)

N = 4000
m = rand(Binomial(N, 0.5))
n = N - m
F, M = rand(distF, m), rand(distM, n)
F0, F1 = F[1,:], F[2,:]
M0, M1 = M[1,:], M[2,:]

@show size(F, 2)
@show size(M, 2)
println()
@show distF
@show mvnF = fit(MvNormal, F)
@show betaF = F0 .^ (0:1)' \ F1
@show sigmaF = √(dot2(F1 - evalpoly.(F0, Ref(betaF)))/(length(F0) - 2))
println("\n")
@show distM
@show mvnM = fit(MvNormal, M)
@show betaM = M0 .^ (0:1)' \ M1
@show sigmaM = √(dot2(M1 - evalpoly.(M0, Ref(betaM)))/(length(M0) - 2))

xs = [minimum([F0; M0]), maximum([F0; M0])]
plot()
scatter!(F0, F1; label="F", msc=:auto, alpha=0.3, ms=2, c=1)
scatter!([mean(F0)], [mean(F1)]; label="", c=1)
plot!(xs, x -> evalpoly(x, betaF); label="", c=:blue)
scatter!(M0, M1; label="M", msc=:auto, alpha=0.3, ms=2, c=2)
scatter!([mean(M0)], [mean(M1)]; label="", c=2)
plot!(xs, x -> evalpoly(x, betaM); label="", c=:red)
plot!(xs, xs; label="", c=:black, ls=:dash, lw=0.5)
plot!(size=(600, 600))

# %%
@show normal0 = fit_mle(Normal, [F0; M0])
histogram([F0; M0]; norm=true, alpha=0.3, bin=50, label="Initial")
plot!(normal0; label="")

# %%
@show normal1 = fit_mle(Normal, [F1; M1])
histogram([F1; M1]; norm=true, alpha=0.3, bin=50, label="Initial")
plot!(normal1; label="")

# %%
plot(x -> pdf(Normal(70, 10), x)/(pdf(Normal(55, 10), x) + pdf(Normal(70, 10), x)), 10, 110; label="M")
plot!(x -> cdf(Logistic((55+70)/2, 6.666), x); label="", ls=:dash)

# %% [markdown]
# $
# (x-a)^2/200 - (x-b)^2/200
# = 2(b-a)x/200 + (a^2-b^2)/200
# = (b-a)x/100 + (a^2-b^2)/200
# $

# %%
100/(b-a)

# %%
σ² = (A[1,1]*A[2,2] - A[1,2]^2)/A[1,1]

# %%
σ = √σ²

# %%
a, b = 55, 70

dist0 = MixtureModel([Normal(a, 10), Normal(b, 10)], [0.5, 0.5])

function distD(x)
    F = pdf(Normal(55, 10), x)
    M = pdf(Normal(70, 10), x)
    p = M/(F+M)
    Bernoulli(p)
end

function dist1(x, d)
    d == 0 ? Normal(0.5x+0.5a, σ) : Normal(0.5x+0.5b, σ)
end

N = 4000
W0 = rand(dist0, N)
W1 = [(d = rand(distD(x)); y = rand(dist1(x, d)); (d, y)) for x in W0]
F = stack([[W0[i], W1[i][2]] for i in eachindex(W0, W1) if W1[i][1] == 0])
M = stack([[W0[i], W1[i][2]] for i in eachindex(W0, W1) if W1[i][1] == 1])

F0, F1 = F[1,:], F[2,:]
M0, M1 = M[1,:], M[2,:]

@show size(F, 2)
@show size(M, 2)
println()
@show distF
@show mvnF = fit(MvNormal, F)
@show betaF = F0 .^ (0:1)' \ F1
@show sigmaF = √(dot2(F1 - evalpoly.(F0, Ref(betaF)))/(length(F0) - 2))
println("\n")
@show distM
@show mvnM = fit(MvNormal, M)
@show betaM = M0 .^ (0:1)' \ M1
@show sigmaM = √(dot2(M1 - evalpoly.(M0, Ref(betaM)))/(length(M0) - 2))

xs = [minimum([F0; M0]), maximum([F0; M0])]
plot()
scatter!(F0, F1; label="F", msc=:auto, alpha=0.3, ms=2, c=1)
scatter!([mean(F0)], [mean(F1)]; label="", c=1)
plot!(xs, x -> evalpoly(x, betaF); label="", c=:blue)
scatter!(M0, M1; label="M", msc=:auto, alpha=0.3, ms=2, c=2)
scatter!([mean(M0)], [mean(M1)]; label="", c=2)
plot!(xs, x -> evalpoly(x, betaM); label="", c=:red)
plot!(xs, xs; label="", c=:black, ls=:dash, lw=0.5)
plot!(size=(600, 600))

# %%
@show normal0 = fit_mle(Normal, [F0; M0])
histogram([F0; M0]; norm=true, alpha=0.3, bin=50, label="Initial")
plot!(normal0; label="")

# %%
@show normal1 = fit_mle(Normal, [F1; M1])
histogram([F1; M1]; norm=true, alpha=0.3, bin=50, label="Initial")
plot!(normal1; label="")

# %%
