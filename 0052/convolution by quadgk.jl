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

# %%
using Distributions
using QuadGK
using StatsPlots
default(fmt=:png)

# %%
function cf(dist::ContinuousUnivariateDistribution, t;
        rtol=1e-4, maxevals=2*10^3, xmin=minimum(dist), xmax=maximum(dist))
    quadgk(x -> pdf(dist, x)*exp(im*t*x), xmin, xmax; rtol, maxevals)[1]
end

function cf2pdf(cf, x; rtol=1e-3, maxevals=10^3, tmin=-Inf, tmax=Inf)
    #g(t) = real(cf(t)) * cos(t*x) + imag(cf(t)) * sin(t*x)
    g(t) = cf(t) * exp(-im*x*t)
    real(quadgk(g, tmin, tmax; rtol, maxevals)[1]) / (2π)
end

function pdf_nconv(dist, n, x; rtol1=1e-4, maxevals1=2*10^3, rtol2=1e-3, maxevals2=10^3, tmin=-1e3, tmax=1e3)
    f(t) = cf(dist, t; rtol=rtol1, maxevals=maxevals1)^n
    cf2pdf(f, x; rtol=rtol2, maxevals=maxevals2, tmin, tmax)
end

# %%
dist = Normal(2, 1)
f(t) = cf(dist, t)
plot(t -> real(f(t)), -5, 5)
plot!(t -> imag(f(t)))

# %%
n = 1
x = 1.0
@time @show pdf_nconv(dist, n, x; maxevals1=20000)
@time @show pdf_nconv(dist, n, x; maxevals1=20000)
@time @show pdf_nconv(dist, n, x; maxevals1=20000)
pdf(dist, x)

# %%
@time plot(range(-5, 10, step=0.5), x -> pdf_nconv(dist, n, x; maxevals1=30000))
plot!(dist)

# %%
dist = Uniform()
f(t) = cf(dist, t)
plot(t -> real(f(t)), -20, 20)
plot!(t -> imag(f(t)))

# %%
dist = Uniform()
n = 5
@time pdf_nconv(dist, n, 1.0)
@time pdf_nconv(dist, n, 1.0)
@time pdf_nconv(dist, n, 1.0)

# %%
n = 6
@time plot(range(-1, n+1, step=0.05), x -> 0 ≤ x ≤ n ? pdf_nconv(dist, n, x) : 0.0)
plot!(Normal(n*mean(dist), √n*std(dist)), -1, n+1; ls=:dash)

# %%
n = 7
@time plot(range(-1, n+1, step=0.05), x -> 0 ≤ x ≤ n ? pdf_nconv(dist, n, x) : 0.0)
plot!(Normal(n*mean(dist), √n*std(dist)), -1, n+1; ls=:dash)

# %%
n = 10
@time plot(range(-1, n+1, step=0.05), x -> 0 ≤ x ≤ n ? pdf_nconv(dist, n, x) : 0.0)
plot!(Normal(n*mean(dist), √n*std(dist)), -1, n+1; ls=:dash)

# %%
