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
using StatsPlots
default(fmt=:png)

dist = Gamma(2, 3)
μ, σ = mean(dist), std(dist)
n = 300

L = 10^6
X̄ = zeros(L)
Threads.@threads for i in 1:L
    X = rand(dist, n)
    X̄[i] = mean(X)
end

P = stephist(X̄; norm=true, label="")
plot!(xlim=(μ-5σ, μ+5σ))
vline!([μ]; label="μ", ls=:dot, c=3)
plot!(ytick=[0])

Q = stephist(X̄; norm=true, label="")
plot!(xlim=(μ-5σ/√n, μ+5σ/√n))
plot!(Normal(μ, σ/√n); label="Normal(μ, σ/√n)", ls=:dash)
vline!([μ]; label="μ", ls=:dot)
plot!(ytick=[0])

plot(P, Q; size=(600, 500), layout=(2, 1))
plot!(plot_title="distribution of sample means for n = $n")

# %%
using Distributions
using StatsPlots
default(fmt=:png)

dist = Poisson(1)
dist = Bernoulli(0.1)
m = 5

μ, σ = mean(dist), std(dist)
n = 1000

L = 10^6
X̄ = zeros(L)
Threads.@threads for i in 1:L
    X = rand(dist, n)
    X̄[i] = mean(X)
end

nbin = 101
binmin, binmax = floor(n*(μ-m*σ/√n))/n, ceil(n*(μ+m*σ/√n))/n
binstep = max(1, round(n*(binmax - binmin)/nbin))/n
bin = binmin-binstep/2:binstep:binmax+binstep/2

P = stephist(X̄; norm=true, bin, label="")
plot!(xlim=(μ-m*σ, μ+m*σ))
vline!([μ]; label="μ", ls=:dot, c=3)
plot!(ytick=[0])

Q = stephist(X̄; norm=true, bin, label="")
plot!(xlim=(μ-m*σ/√n, μ+m*σ/√n))
plot!(Normal(μ, σ/√n); label="Normal(μ, σ/√n)", ls=:dash)
vline!([μ]; label="μ", ls=:dot)
plot!(ytick=[0])

plot(P, Q; size=(600, 500), layout=(2, 1))
plot!(plot_title="distribution of sample means for n = $n")

# %%
using Distributions
using Random
using StatsPlots
default(fmt=:png)

function plot_lln_clt(; dist=Gamma(2, 3), m=5, nbin=200, L=10^6)
    μ, σ = mean(dist), std(dist)
    ns = round.(Int, exp.(range(log(1), log(1000), 200)))
    nmax = maximum(ns)

    X̄ = zeros(L, length(ns))
    Xtmp = [zeros(eltype(dist), nmax) for _ in 1:Threads.nthreads()]
    @time Threads.@threads for i in 1:L
        X = rand!(dist, Xtmp[Threads.threadid()])
        for (j, n) in enumerate(ns)
            X̄[i, j] = mean(@view(X[1:n]))
        end
    end

    js = collect(eachindex(ns))
    js_gif = [fill(js[begin], 60); js; fill(js[end], 60)]
    ns_gif = ns[js_gif]
    @time @gif for (j, n) in zip(js_gif, ns_gif)
        bin = if dist isa DiscreteUnivariateDistribution
            binmin, binmax = floor(n*(μ-m*σ/√n))/n, ceil(n*(μ+m*σ/√n))/n
            binstep = max(1, round(n*(binmax - binmin)/nbin))/n
            binmin-binstep/2:binstep:binmax+binstep/2
        else
            :auto
        end
        
        P = stephist(@view(X̄[:, j]); norm=true, bin, label="")
        plot!(xlim=(μ-m*σ, μ+m*σ))
        vline!([μ]; label="μ", ls=:dot, c=3)
        title!("law of large numbers (n = $n)"; titlefontsize=12)
        plot!(ytick=[0])

        Q = stephist(@view(X̄[:, j]); norm=true, bin, label="")
        plot!(xlim=(μ-m*σ/√n, μ+m*σ/√n))
        plot!(Normal(μ, σ/√n), μ-m*σ/√n, μ+m*σ/√n; label="Normal(μ, σ/√n)", ls=:dash)
        vline!([μ]; label="μ", ls=:dot)
        title!("central limit theorem (n = $n)"; titlefontsize=12)
        plot!(ytick=[0])

        plot(P, Q; size=(600, 500), layout=(2, 1))
        plot!(plot_title="distribution of sample means for n = $n",
            plot_titlefontsize=14)
    end
end

# %%
plot_lln_clt(; dist = Uniform())

# %%
plot_lln_clt(; dist = Gamma(2, 3))

# %%
plot_lln_clt(; dist = MixtureModel([Normal(), Normal(10)], [0.9, 0.1]))

# %%
plot_lln_clt(; dist = MixtureModel([Normal(), Normal(20)], [0.95, 0.05]))

# %%
plot_lln_clt(; dist = Poisson(1))

# %%
plot_lln_clt(; dist = Bernoulli(0.1))

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
using Distributions
using Random
using StatsPlots
default(fmt=:png)

function plot_lln_clt_meanvar2x2(; dist=Gamma(2, 3), m=5, nbin=200, nns=200, L=10^6)
    μ, σ = mean(dist), std(dist)
    sk, ku = skewness(dist), kurtosis(dist)
    ns = round.(Int, exp.(range(log(1), log(1000), nns)))
    nmax = maximum(ns)

    X̄ = zeros(L, length(ns))
    S² = zeros(L, length(ns))
    Xtmp = [zeros(eltype(dist), nmax) for _ in 1:Threads.nthreads()]
    @time Threads.@threads for i in 1:L
        X = rand!(dist, Xtmp[Threads.threadid()])
        for (j, n) in enumerate(ns)
            XX = @view X[1:n]
            X̄[i, j] = mean(XX)
            S²[i, j] = n == 1 ? (XX[1] - μ)^2 : var(XX)
        end
    end

    js = collect(eachindex(ns))
    js_gif = [fill(js[1], 60); js; fill(js[end], 60)]
    ns_gif = ns[js_gif]
    
    @time @gif for (j, n) in zip(js_gif, ns_gif)
        bin = if dist isa DiscreteUnivariateDistribution
            binmin, binmax = floor(n*(μ-m*σ/√n))/n, ceil(n*(μ+m*σ/√n))/n
            binstep = max(1, round(n*(binmax - binmin)/nbin))/n
            binmin-binstep/2:binstep:binmax+binstep/2
        else
            :auto
        end
        
        P1 = stephist(@view(X̄[:, j]); norm=true, bin, label="")
        plot!(xlim=(μ-m*σ, μ+m*σ))
        vline!([μ]; label="μ", ls=:dot, c=3)
        title!("law of large numbers for sample mean (n=$n)"; titlefontsize=10)
        plot!(ytick=[0])

        Q1 = stephist(@view(X̄[:, j]); norm=true, bin, label="")
        plot!(xlim=(μ-m*σ/√n, μ+m*σ/√n))
        plot!(Normal(μ, σ/√n), μ-m*σ/√n, μ+m*σ/√n; label="Normal(μ, σ/√n)", ls=:dash)
        vline!([μ]; label="μ", ls=:dot)
        title!("central limit theorem for sample mean (n=$n)"; titlefontsize=10)
        plot!(ytick=[0])
        
        P2 = stephist(@view(S²[:, j]); norm=true, label="")
        plot!(xlim=(σ^2-m*σ^2*√(ku+2), σ^2+m*σ^2*√(ku+2)))
        vline!([σ^2]; label="σ²", ls=:dot, c=3)
        title!("law of large numbers for unbiased variance (n = $n)"; titlefontsize=10)
        plot!(ytick=[0])

        std_S² = σ^2 * √(ku/n + 2/max(1, n-1))
        #alpha = σ^4/std_S²^2
        #theta = std_S²^2/σ^2
        Q2 = stephist(@view(S²[:, j]); norm=true, label="")
        plot!(xlim=(σ^2-m*std_S², σ^2+m*std_S²))
        plot!(Normal(σ^2, std_S²), σ^2-m*std_S², σ^2+m*std_S²; label="normal approx.", ls=:dash)
        vline!([σ^2]; label="σ²", ls=:dot)
        #plot!(Gamma(alpha, theta), σ^2-m*std_S², σ^2+m*std_S²; label="gamma aprrox.", ls=:dashdot)
        title!("central limit theorem for unbiased variance (n = $n)"; titlefontsize=10)
        plot!(ytick=[0])
        
        plot(P1, P2, Q1, Q2; size=(1000, 500), layout=(2, 2))
    end
end

# %%
plot_lln_clt_meanvar2x2(; dist=Gamma(2, 3))

# %%
plot_lln_clt_meanvar2x2(; dist=MixtureModel([Normal(), Normal(20)], [0.95, 0.05]))

# %%
