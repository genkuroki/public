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

function meanvarstdskku(dist)
    m = mean(dist)
    s2 = var(dist)
    s = √s2
    sk = skewness(dist)
    ku = kurtosis(dist)
    m, s2, s, sk, ku
end

function meanvarstdskku(dist::MixtureModel)
    M1, M2, M3, M4 = moment1234(dist)
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

# Warning: type piracy
Distributions.skewness(dist::MixtureModel) = meanvarstdskku(dist)[4]
Distributions.kurtosis(dist::MixtureModel) = meanvarstdskku(dist)[5]

# %%
mix = MixtureModel([Normal(), Normal(20)], [0.95, 0.05])
@show skewness(mix) kurtosis(mix);

# %%
mixdisc = MixtureModel([Poisson(), Poisson(10)], [0.9, 0.1])
@show skewness(mixdisc) kurtosis(mixdisc);

# %%
mixmix = MixtureModel([mix, Normal(2, 3)], [0.7, 0.3])
@show skewness(mixmix) kurtosis(mixmix);

# %%
using Distributions
using QuadGK

function skewness_approx(dist::ContinuousUnivariateDistribution)
    μ = mean(dist)
    σ = std(dist)
    quadgk(x -> ((x-μ)/σ)^3 * pdf(dist, x), extrema(dist)...)[1]
end

function skewness_approx(dist::DiscreteUnivariateDistribution)
    μ = mean(dist)
    σ = std(dist)
    xmin = max(round(Int, μ-100σ), minimum(dist))
    xmax = min(round(Int, μ+100σ), maximum(dist))
    sum(x -> ((x-μ)/σ)^3 * pdf(dist, x), xmin:xmax)
end

function kurtosis_approx(dist::ContinuousUnivariateDistribution)
    μ = mean(dist)
    σ = std(dist)
    quadgk(x -> ((x-μ)/σ)^4 * pdf(dist, x), extrema(dist)...)[1] - 3
end

function kurtosis_approx(dist::DiscreteUnivariateDistribution)
    μ = mean(dist)
    σ = std(dist)
    xmin = max(round(Int, μ-100σ), minimum(dist))
    xmax = min(round(Int, μ+100σ), maximum(dist))
    sum(x -> ((x-μ)/σ)^4 * pdf(dist, x), xmin:xmax) - 3
end

function meanvarstdskku_approx(dist)
    m = mean(dist)
    s2 = var(dist)
    s = √s2
    sk = skewness_approx(dist)
    ku = kurtosis_approx(dist)
    m, s2, s, sk, ku
end

function _moment_approx(dist::ContinuousUnivariateDistribution, k)
    quadgk(x -> x^k * pdf(dist, x), extrema(dist)...)[1]
end

function _moment_approx(dist::DiscreteUnivariateDistribution, k)
    μ = mean(dist)
    σ = std(dist)
    xmin = max(round(Int, μ-100σ), minimum(dist))
    xmax = min(round(Int, μ+100σ), maximum(dist))
    sum(x -> x^k * pdf(dist, x), xmin:xmax)
end

function moment1234_approx(dist)
    m1, m2, m3, m4 = (_moment_approx(dist, k) for k in 1:4)
    m1, m2, m3, m4
end

# %%
dist = Gamma(2, 3)
@show meanvarstdskku(dist) meanvarstdskku_approx(dist)
@show moment1234(dist) moment1234_approx(dist);

# %%
dist = Poisson(30)
@show meanvarstdskku(dist) meanvarstdskku_approx(dist)
@show moment1234(dist) moment1234_approx(dist);

# %%
mix = MixtureModel([Normal(), Normal(20)], [0.95, 0.05])
@show moment1234(mix)
@show moment1234_approx(mix)
@show meanvarstdskku(mix)
@show meanvarstdskku_approx(mix);

# %%
mixdisc = MixtureModel([Poisson(1), Poisson(10)], [0.9, 0.1])
@show moment1234(mixdisc)
@show moment1234_approx(mixdisc)
@show meanvarstdskku(mixdisc)
@show meanvarstdskku_approx(mixdisc);

# %%
mixmix = MixtureModel([mix, Normal(2, 3)], [0.7, 0.3])
@show moment1234(mixmix)
@show moment1234_approx(mixmix)
@show meanvarstdskku(mixmix)
@show meanvarstdskku_approx(mixmix);

# %%
