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
#     display_name: Julia 1.9.3
#     language: julia
#     name: julia-1.9
# ---

# %%
using Distributions
using QuadGK

function skewness_quadgk(dist)
    μ = mean(dist)
    σ = std(dist)
    quadgk(x -> ((x-μ)/σ)^3 * pdf(dist, x), extrema(dist)...)[1]
end

function kurtosis_quadgk(dist)
    μ = mean(dist)
    σ = std(dist)
    quadgk(x -> ((x-μ)/σ)^4 * pdf(dist, x), extrema(dist)...)[1] - 3
end

function meanvarstdskku(dist)
    m = mean(dist)
    s2 = var(dist)
    s = √s2
    sk = skewness(dist)
    ku = kurtosis(dist)
    m, s2, s, sk, ku
end

function meanvarstdskku_quadgk(dist)
    m = mean(dist)
    s2 = quadgk(x -> (x - m)^2 * pdf(dist, x), extrema(dist)...)[1]
    s = √s2
    sk = quadgk(x -> ((x - m)/s)^3 * pdf(dist, x), extrema(dist)...)[1]
    ku = quadgk(x -> ((x - m)/s)^4 * pdf(dist, x), extrema(dist)...)[1] - 3
    m, s2, s, sk, ku
end

function meanvarstdskku(dist::MixtureModel)
    M1, M2, M3, M4 = moment1234_quadgk(dist)
    m = mean(dist)
    s2 = var(dist)
    s = √s2
    s3, s4 = s*s2, s2^2
    sk = 1/s3 * (M3 - 3m*s2 - m^3)
    ku = 1/s4 * (M4 - 4m*s3*sk - 6m^2*s2 - m^4) - 3
    m, s2, s, sk, ku
end

function moment1234(dist)
    m, s2, s, sk, ku = meanvarstdskku(dist)
    s3, s4 = s*s2, s2^2
    m1 = m
    m2 = s2 + m^2
    m3 = s3*sk + 3m*s2 + m^3
    m4 = s4*(ku + 3) + 4m*s3*sk + 6m^2*s2 + m^4
    m1, m2, m3, m4
end

function moment1234_quadgk(dist)
    m1 = quadgk(x -> x * pdf(dist, x), extrema(dist)...)[1]
    m2 = quadgk(x -> x^2 * pdf(dist, x), extrema(dist)...)[1]
    m3 = quadgk(x -> x^3 * pdf(dist, x), extrema(dist)...)[1]
    m4 = quadgk(x -> x^4 * pdf(dist, x), extrema(dist)...)[1]
    m1, m2, m3, m4
end

function moment1234(dist::MixtureModel)
    T = float(eltype(dist))
    M1, M2, M3, M4 = zero(T), zero(T), zero(T), zero(T)
    for (d, q) in zip(components(dist), probs(dist))
        m1, m2, m3, m4 = moment1234(d)
        M1 += q * m1
        M2 += q * m2
        M3 += q * m3
        M4 += q * m4
    end
    M1, M2, M3, M4
end

# %%
T = float(eltype(Bernoulli(0.5)))
zero(T)

# %%
dist = Gamma(2, 3)
@show meanvarstdskku(dist) meanvarstdskku_quadgk(dist)
@show moment1234(dist) moment1234_quadgk(dist);

# %%
mix = MixtureModel([Normal(), Normal(20)], [0.95, 0.05])
@show moment1234(mix)
@show moment1234_quadgk(mix)
@show meanvarstdskku(mix)
@show meanvarstdskku_quadgk(mix);

# %%
