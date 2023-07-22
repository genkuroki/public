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
#     display_name: Julia 1.9.2
#     language: julia
#     name: julia-1.9
# ---

# %%
using Distributions
using Random
using StatsPlots
default(fmt=:png)

coverage_probability(A, α) = count(≥(α), A)/length(A)

# %%
function sim(dist, n; L=10^6)
    μ, σ² = mean(dist), var(dist)
    nths = Threads.nthreads()
    Xtmp = [Vector{Float64}(undef, n) for _ in 1:nths]
    X̄ = Vector{Float64}(undef, L)
    S² = similar(X̄)
    T = similar(X̄)
    Chi² = similar(X̄)
    pvalT = similar(X̄)
    pvalChi² = similar(X̄)
    Threads.@threads for i in 1:L
        tid = Threads.threadid()
        X = rand!(dist, Xtmp[tid])
        x̄, s² = mean(X), var(X)
        t = (x̄ - μ)/√(s²/n)
        chi² = (n-1)*s²/σ²
        X̄[i] = x̄
        S²[i] = s²
        T[i] = t
        Chi²[i] = chi²
        pvalT[i] = 2ccdf(TDist(n-1), abs(t))
        pvalChi²[i] = min(1, 2cdf(Chisq(n-1), chi²), 2ccdf(Chisq(n-1), chi²))
    end
    X̄, S², T, Chi², pvalT, pvalChi²
end

# %%
dist = Normal(2, 3)
n = 10

plot(dist; label="Normal(2, 3)", legendfontsize=10)
title!("sample size: $n", titlefontsize=10)
plot!(size=(400, 250)) |> display

X̄, S², T, Chi², pvalT, pvalChi² = @time sim(dist, n)
@show coverage_probability(pvalT, 0.05) coverage_probability(pvalChi², 0.05)

P = histogram(T; normed=true, alpha=0.2, label="sample t-values")
plot!(TDist(n-1); label="TDist(n-1)")
plot!(xlim=(-6, 6))

Q = histogram(Chi²; normed=true, alpha=0.2, label="sample (n-1)s²/σ²")
plot!(Chisq(n-1), label="Chisq(n-1)")

R = plot(α -> coverage_probability(pvalT, α), 0, 1;
    label="population mean")
plot!(α -> 1 - α; label="", ls=:dot)
plot!(xtick=0:0.1:1, ytick=0:0.1:1)
plot!(xguide="α", yguide="coverage probability")

S = plot(α -> coverage_probability(pvalChi², α), 0, 1;
    label="population variance")
plot!(α -> 1 - α; label="", ls=:dot)
plot!(xtick=0:0.1:1, ytick=0:0.1:1)
plot!(xguide="α", yguide="coverage probability")

plot(P, Q, R, S; size=(800, 500))

# %%
dist = Normal(2, 3)
n = 100

plot(dist; label="Normal(2, 3)", legendfontsize=10)
title!("sample size: $n", titlefontsize=10)
plot!(size=(400, 250)) |> display

X̄, S², T, Chi², pvalT, pvalChi² = @time sim(dist, n)
@show coverage_probability(pvalT, 0.05) coverage_probability(pvalChi², 0.05)

P = histogram(T; normed=true, alpha=0.2, label="sample t-values")
plot!(TDist(n-1); label="TDist(n-1)")
plot!(xlim=(-6, 6))

Q = histogram(Chi²; normed=true, alpha=0.2, label="sample (n-1)s²/σ²")
plot!(Chisq(n-1), label="Chisq(n-1)")
plot!(xlim=(max(0, n-8√n), n+8√n))

R = plot(α -> coverage_probability(pvalT, α), 0, 1;
    label="population mean")
plot!(α -> 1 - α; label="", ls=:dot)
plot!(xtick=0:0.1:1, ytick=0:0.1:1)
plot!(xguide="α", yguide="coverage probability")

S = plot(α -> coverage_probability(pvalChi², α), 0, 1;
    label="population variance")
plot!(α -> 1 - α; label="", ls=:dot)
plot!(xtick=0:0.1:1, ytick=0:0.1:1)
plot!(xguide="α", yguide="coverage probability")

plot(P, Q, R, S; size=(800, 500))

# %%
dist = Gamma(2, 3)
n = 100

plot(dist; label="Gamma(2, 3)", legendfontsize=10)
title!("sample size: $n", titlefontsize=10)
plot!(size=(400, 250)) |> display

X̄, S², T, Chi², pvalT, pvalChi² = @time sim(dist, n)
@show coverage_probability(pvalT, 0.05) coverage_probability(pvalChi², 0.05)

P = histogram(T; normed=true, alpha=0.2, label="sample t-values")
plot!(TDist(n-1); label="TDist(n-1)")
plot!(xlim=(-6, 6))

Q = histogram(Chi²; normed=true, alpha=0.2, label="sample (n-1)s²/σ²")
plot!(Chisq(n-1), label="Chisq(n-1)")
plot!(xlim=(max(0, n-8√n), n+8√n))

R = plot(α -> coverage_probability(pvalT, α), 0, 1;
    label="population mean")
plot!(α -> 1 - α; label="", ls=:dot)
plot!(xtick=0:0.1:1, ytick=0:0.1:1)
plot!(xguide="α", yguide="coverage probability")

S = plot(α -> coverage_probability(pvalChi², α), 0, 1;
    label="population variance")
plot!(α -> 1 - α; label="", ls=:dot)
plot!(xtick=0:0.1:1, ytick=0:0.1:1)
plot!(xguide="α", yguide="coverage probability")

plot(P, Q, R, S; size=(800, 500))

# %%
