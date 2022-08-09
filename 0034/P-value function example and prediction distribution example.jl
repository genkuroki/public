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
using Distributions
using Random
using Roots
using StatsPlots
default(fmt = :png, size = (640, 400),
    titlefontsize = 10, plot_titlefont = 10)

# %%
function confdist(m, x̄, sx²)
    SEhat = √(sx²/m)
    x̄ + SEhat * TDist(m - 1)
end

function confdist(x)
    m, x̄, sx² = length(x), mean(x), var(x)
    confdist(m, x̄, sx²)
end

function tvalue(m, x̄, sx², μ = 0.0)
    (x̄ - μ) / √(sx²/m)
end

function tvalue(x, μ = 0.0)
    m, x̄, sx² = length(x), mean(x), var(x)
    tvalue(m, x̄, sx², μ)
end

function pvalue(m, x̄, sx², μ = 0.0)
    t = tvalue(m, x̄, sx², μ)
    ν = m - 1
    2ccdf(TDist(ν), abs(t))
end

function pvalue(x, μ = 0.0)
    m, x̄, sx² = length(x), mean(x), var(x)
    pvalue(m, x̄, sx², μ)
end

function confint(m, x̄, sx²; α = 0.05)
    ν = m - 1
    c = quantile(TDist(ν), 1-α/2)
    SEhat = √(sx²/m)
    [x̄ - c*SEhat, x̄ + c*SEhat]
end

function confint(x; α = 0.05)
    m, x̄, sx² = length(x), mean(x), var(x)
    confint(m, x̄, sx²; α)
end

# %%
function dist_pred(m, x̄, sx²)
    SEhat = √(sx²*(1+1/m))
    x̄ + SEhat * TDist(m - 1)
end

function dist_pred(x)
    m, x̄, sx² = length(x), mean(x), var(x)
    dist_pred(m, x̄, sx²)
end

function tvalue_pred(m, x̄, sx², xnew)
    (x̄ - xnew) / √(sx²*(1+1/m))
end

function tvalue_pred(x, xnew)
    m, x̄, sx² = length(x), mean(x), var(x)
    tvalue_pred(m, x̄, sx², xnew)
end

function pvalue_pred(m, x̄, sx², xnew)
    t = tvalue_pred(m, x̄, sx², xnew)
    ν = m - 1
    2ccdf(TDist(ν), abs(t))
end

function pvalue_pred(x, xnew)
    m, x̄, sx² = length(x), mean(x), var(x)
    pvalue_pred(m, x̄, sx², xnew)
end

function predint(m, x̄, sx²; α = 0.05)
    ν = m - 1
    c = quantile(TDist(ν), 1-α/2)
    SEhat = √(sx²*(1+1/m))
    [x̄ - c*SEhat, x̄ + c*SEhat]
end

function predint(x; α = 0.05)
    m, x̄, sx² = length(x), mean(x), var(x)
    confint_pred(m, x̄, sx²; α)
end

# %%
function logmaxlik(x)
    mle = fit_mle(Normal, x)
    loglikelihood(mle, x)
end

function maxlik(x)
    exp(logmaxlik(x))
end

function conditional_mle(x, μ)
    s = std(x; mean = μ, corrected = false)
    Normal(μ, s)
end

function conditional_logmaxlik(x, μ)
    cmle = conditional_mle(x, μ)
    loglikelihood(cmle, x)
end

function conditional_maxlik(x, μ)
    exp(conditional_logmaxlik(x, μ))
end

function logmaxlikrat(x, μ)
    conditional_logmaxlik(x, μ) - logmaxlik(x)
end

function maxlikrat(x, μ)
    exp(logmaxlikrat(x, μ))
end

function neg2loglikrat(x, μ)
    -2logmaxlikrat(x, μ)
end

function pvalue_loglikrat(x, μ)
    ccdf(Chisq(1), neg2loglikrat(x, μ))
end

function confint_loglikrat(x; α = 0.05)
    c = quantile(Chisq(1), 1 - α)
    ci = confint(x; α = α/10)
    find_zeros(μ -> neg2loglikrat(x, μ) - c, ci...)
end

# %%
logmaxlikrat2pvalue(lmlr) = ccdf(Chisq(1), -2lmlr)
maxlikrat2pvalue(mlr) = ccdf(Chisq(1), -2log(mlr))
pvalue2logmaxlikrat(pval) = -quantile(Chisq(1), 1 - pval)/2
pvalue2maxlikrat(pval) = exp(pvalue2logmaxlikrat(pval))

# %%
function baysian_update(x;
        μ₀ = 0, sμ² = 100^2, γ = 2, θ = 0.5,
        pri = (μ₀ = μ₀, λ₀ = 1/(γ*θ*(1 - 1/γ)*sμ²), γ, θ)
    )
    (; μ₀, λ₀, γ, θ) = pri
    m, x̄, s² = length(x), mean(x), var(x)
    μ₀_new = (λ₀*μ₀ + m*x̄)/(λ₀ + m)
    λ₀_new = λ₀ + m
    γ_new = γ + m/2
    θ_new = θ/(1 + θ/2*(m*λ₀/(λ₀ + m)*(x̄ - μ₀)^2 + (m - 1)*s²))
    post = (μ₀ = μ₀_new, λ₀ = λ₀_new, γ = γ_new, θ = θ_new)
end

function pdf_μλ(μ, λ;
        μ₀ = 0, sμ² = 100^2, γ = 2, θ = 0.5,
        pri = (μ₀ = μ₀, λ₀ = 1/(γ*θ*(1 - 1/γ)*sμ²), γ, θ)
    )
    (; μ₀, λ₀, γ, θ) = pri
    pdf(Normal(μ₀, 1/√(λ*λ₀)), μ) * pdf(Gamma(γ, θ), λ)
end

function dist_μ(;
        μ₀ = 0, sμ² = 100^2, γ = 2, θ = 0.5,
        pri = (μ₀ = μ₀, λ₀ = 1/(γ*θ*(1 - 1/γ)*sμ²), γ, θ)
    )
    (; μ₀, λ₀, γ, θ) = pri
    μ₀ + TDist(2γ)/√(λ₀*γ*θ)
end

function dist_λ(;
        μ₀ = 0, sμ² = 100^2, γ = 2, θ = 0.5,
        pri = (μ₀ = μ₀, λ₀ = 1/(γ*θ*(1 - 1/γ)*sμ²), γ, θ)
    )
    (; μ₀, λ₀, γ, θ) = pri
    Gamma(γ, θ)
end

function dist_x(; 
        μ₀ = 0, sμ² = 100^2, γ = 2, θ = 0.5,
        pri = (μ₀ = μ₀, λ₀ = 1/(γ*θ*(1 - 1/γ)*sμ²), γ, θ)
    )
    (; μ₀, λ₀, γ, θ) = pri
    μ₀ + TDist(2γ)*√((1 + 1/λ₀)/(γ*θ))
end

function rand_μλ(; 
        μ₀ = 0, sμ² = 100^2, γ = 2, θ = 0.5,
        pri = (μ₀ = μ₀, λ₀ = 1/(γ*θ*(1 - 1/γ)*sμ²), γ, θ)
    )
    (; μ₀, λ₀, γ, θ) = pri
    λ = rand(Gamma(γ, θ))
    μ = rand(Normal(μ₀, 1/√(λ*λ₀)))
    (μ, λ)
end

function rand_μλ(L; 
        μ₀ = 0, sμ² = 100^2, γ = 2, θ = 0.5,
        pri = (μ₀ = μ₀, λ₀ = 1/(γ*θ*(1 - 1/γ)*sμ²), γ, θ)
    )
    [rand_μλ(; pri) for _ in 1:L]
end

function pvalue_bayes(x, μ;
        μ₀ = 0, sμ² = 100^2, γ = 2, θ = 0.5,
        pri = (μ₀ = μ₀, λ₀ = 1/(γ*θ*(1 - 1/γ)*sμ²), γ, θ)
    )
    post = baysian_update(x; pri)
    post_μ = dist_μ(; pri=post)
    2ccdf(post_μ, mean(post_μ) + abs(μ - mean(post_μ)))
end

function credint(x; α = 0.05,
        μ₀ = 0, sμ² = 100^2, γ = 2, θ = 0.5,
        pri = (μ₀ = μ₀, λ₀ = 1/(γ*θ*(1 - 1/γ)*sμ²), γ, θ)
    )
    post = baysian_update(x; pri)
    post_μ = dist_μ(; pri=post)
    quantile.(post_μ, [α/2, 1-α/2])
end

# %%
function plot_etc(;
        dist_true = Gamma(10, 1),
        m = 10,
        ms = 3,
        α = 0.05,
        μ₀ = 0, sμ² = 100^2, γ = 2, θ = 0.5,
        pri = (μ₀ = μ₀, λ₀ = 1/(γ*θ*(1 - 1/γ)*sμ²), γ, θ),
        seed = nothing
    )
    isnothing(seed) || Random.seed!(seed)
    x = rand(dist_true, m)
    
    @show mean(dist_true) var(dist_true)
    @show mean(x) var(x)
    flush(stdout)
    
    μ, σ = mean(dist_true), std(dist_true)
    a, b = max(minimum(dist_true), μ-3σ), min(maximum(dist_true), μ+6σ)
    
    plot(dist_true, a, b; label="true dist")
    h = pdf(dist_true, mode(dist_true))
    scatter!(x, fill(-0.05h, m); label="sample", ms, msc=:auto, alpha=0.5, c=:red)
    title!("true distribution and size-$m sample")
    plot!() |> display
    
    plot(μ₀ -> pvalue(x, μ₀), a, b; label="P-value function of t-test")
    scatter!(x, fill(-0.05, m); label="sample", ms, msc=:auto, alpha=0.5, c=:red)
    plot!(xguide="μ₀", ytick=0:0.05:1)
    title!("P-value function of hypothesis μ = μ₀ for size-$m sample")
    plot!() |> display
    
    plot(μ₀ -> pvalue_loglikrat(x, μ₀), a, b;
        label="P-value function of log max. lik. ratio test", c=2)
    scatter!(x, fill(-0.05, m); label="sample", ms, msc=:auto, alpha=0.5, c=:red)
    plot!(xguide="μ₀", ytick=0:0.05:1)
    title!("P-value function of hypothesis μ = μ₀ for size-$m sample")
    plot!() |> display
    
    plot(μ₀ -> pvalue(x, μ₀), a, b; label="P-value function of t-test")
    plot!(μ₀ -> pvalue_loglikrat(x, μ₀), a, b;
        label="P-value function of log max. lik. ratio test", c=2, ls=:dash)
    scatter!(x, fill(-0.05, m); label="sample", ms, msc=:auto, alpha=0.5, c=:red)
    plot!(xguide="μ₀", ytick=0:0.05:1)
    title!("P-value function of hypothesis μ = μ₀ for size-$m sample")
    plot!() |> display
    
    @show ci = confint(x; α)
    flush(stdout)

    plot(μ₀ -> pvalue(x, μ₀), a, b; label="P-value function of t-test")
    plot!(ci, fill(α, 2); label="$(100(1-α))% CI of μ", lw=3, c=:blue)
    scatter!(x, fill(-0.05, m); label="sample", ms, msc=:auto, alpha=0.5, c=:red)
    plot!(xguide="μ₀", ytick=0:0.05:1)
    title!("P-value function of hypothesis μ = μ₀ and \
        confidence interval for size-$m sample")
    plot!() |> display
    
    @show ci = confint_loglikrat(x; α)
    flush(stdout)

    plot(μ₀ -> pvalue_loglikrat(x, μ₀), a, b;
        label="P-value function of log max. lik. ratio test", c=2)
    plot!(ci, fill(α, 2); label="$(100(1-α))% CI of μ", lw=3, c=:magenta)
    scatter!(x, fill(-0.05, m); label="sample", ms, msc=:auto, alpha=0.5, c=:red)
    plot!(xguide="μ₀", ytick=0:0.05:1)
    title!("P-value function of hypothesis μ = μ₀ and \
        confidence interval for size-$m sample")
    plot!() |> display
    
    @show ci = confint_loglikrat(x; α)
    flush(stdout)

    plot(μ₀ -> maxlikrat(x, μ₀), a, b; label="maximum likelihood ratio function", c=2)
    c = quantile(Chisq(1), 1-α)
    plot!(ci, fill(exp(-c/2), 2); label="$(100(1-α))% CI of μ", lw=3, c=:magenta)
    scatter!(x, fill(-0.05, m); label="sample", ms=3, msc=:auto, alpha=0.5, c=:red)
    plot!(xguide="μ₀", ytick=0:0.05:1)
    title!("maximum likelihood ratio function of μ = μ₀ for size-$m sample")
    plot!() |> display
    
    plot(μ₀ -> pvalue2maxlikrat(pvalue(x, μ₀)), a, b;
        label="max. lik. rat. conversion of t-test P-value func.")
    plot!(μ₀ -> maxlikrat(x, μ₀), a, b;
        label="maximum likelihood ratio function", c=2, ls=:dash)
    scatter!(x, fill(-0.05, m);
        label="sample", ms, msc=:auto, alpha=0.5, c=:red)
    plot!(xguide="μ₀", ytick=0:0.05:1)
    title!("maximum likelihood ratio function of μ = μ₀ for size-$m sample")
    plot!() |> display
    
    @show pri
    μλ = rand_μλ(10^5; pri)
    μσ² = [(μ, 1/λ) for (μ, λ) in μλ]
    μ, σ² = first.(μσ²), last.(μσ²)
    @show mean(μ), mean(σ²)
    xlim = quantile.(Ref(μ), (0.001, 0.999))
    ylim = quantile.(Ref(σ²), (0.001, 0.999))
    flush(stdout)

    scatter(μσ²; label="", ms=1, msc=:auto, alpha=0.1)
    plot!(xguide="μ", yguide="σ²")
    plot!(; xlim, ylim, yscale=:log10)
    title!("prior of (μ, σ²)")
    plot!() |> display
    
    @show post = baysian_update(x; pri)
    μλ = rand_μλ(10^4; pri=post)
    μσ² = [(μ, 1/λ) for (μ, λ) in μλ]
    μ, σ² = first.(μσ²), last.(μσ²)
    @show mean(μ), mean(σ²)
    xlim = quantile.(Ref(μ), (0.001, 0.999))
    ylim = quantile.(Ref(σ²), (0.001, 0.999))
    flush(stdout)
    
    scatter(μσ²; label="", ms=1, msc=:auto, alpha=0.3)
    plot!(xguide="μ", yguide="σ²")
    plot!(; xlim, ylim, yscale=:log10)
    title!("posterior of (μ, σ²) for size-$m sample")
    plot!() |> display
    
    @show post_μ = dist_μ(; pri=post)
    @show ci = credint(x; α, pri)
    flush(stdout)
    
    plot(post_μ, a, b; label="posterior")
    c = pdf(post_μ, mode(post_μ))
    plot!(μ₀ -> c*maxlikrat(x, μ₀), a, b;
        label="normalized maximum likelihood ratio function",
        c=2, ls=:dash)
    scatter!(x, fill(-0.05c, 2); label="sample",
        ms, msc=:auto, alpha=0.5, c=:red)
    title!("posterior pdf function of μ for size-$m sample")
    plot!() |> display
end

# %%
plot_etc(; dist_true=Gamma(10, 1), m=10, ms=3, seed=4649373)

# %%
plot_etc(; dist_true=Gamma(10, 1), m=20, ms=3, seed=4649373)

# %%
plot_etc(; dist_true=Gamma(10, 1), m=30, ms=3, seed=4649373)

# %%
plot_etc(; dist_true=Gamma(10, 1), m=40, ms=3, seed=4649373)

# %%
plot_etc(; dist_true=Gamma(10, 1), m=80, ms=2, seed=4649373)

# %%
plot_etc(; dist_true=Gamma(10, 1), m=160, ms=1, seed=4649373)

# %%
plot_etc(; dist_true=Gamma(10, 1), m=320, ms=1, seed=4649373)

# %%
plot_etc(; dist_true=Gamma(10, 1), m=640, ms=1, seed=4649373)

# %%
function plot_ttest(;
        dist_true = Gamma(10, 1),
        m = 10,
        μ_0 = 0, sμ² = 100^2, γ = 1.5, θ = 2,
        pri = (μ₀ = 0, λ₀ = 1/(γ*θ*(1 - 1/γ)*sμ²), γ, θ),
        seed = 4649373
    )
    isnothing(seed) || Random.seed!(seed)
    x = rand(dist_true, m)

    μ, σ = mean(dist_true), std(dist_true)
    a, b = max(minimum(dist_true), μ - 5σ), min(maximum(dist_true), μ + 6σ)

    post = baysian_update(x; pri)
    post_μ = dist_μ(; pri=post)
    pred_xnew = dist_x(; pri=post)
    μ_xnew = mean(pred_xnew)
    
    @show dist_true
    @show m
    @show mean(dist_true) var(dist_true)
    @show mean(x) var(x)
    @show mean(dist_pred(x)) var(dist_pred(x))
    @show mean(pred_xnew) var(pred_xnew)
    
    P1 = plot(μ -> pvalue(x, μ), a, b; label="t-test")
    plot!(μ -> pvalue_bayes(x, μ), a, b; label="Bayesian", ls=:dash)
    scatter!(x, fill(-0.05, length(x));
        label="sample", ms=1.5, msc=:auto, alpha=0.5, c=:red)
    title!("P-value functions")

    P2 = plot(μ -> pdf(confdist(x), μ), a, b; label="t-test")
    plot!(μ -> pdf(post_μ, μ), a, b; label="Bayesian", ls=:dash)
    h = pdf(confdist(x), mode(confdist(x)))
    scatter!(x, fill(-0.05h, length(x));
        label="sample", ms=1.5, msc=:auto, alpha=0.5, c=:red)
    title!("parameter distributions")

    P3 = plot(xnew -> pvalue_pred(x, xnew), a, b; label="t-test")
    plot!(xnew -> 2ccdf(pred_xnew, μ_xnew + abs(xnew - μ_xnew)), a, b;
        label="Bayesian", ls=:dash)
    scatter!(x, fill(-0.05, length(x));
        label="sample", ms=1.5, msc=:auto, alpha=0.5, c=:red)
    title!("prediction P-value functions")

    P4 = plot(xnew -> pdf(dist_pred(x), xnew), a, b; label="t-test")
    plot!(xnew -> pdf(pred_xnew, xnew), a, b; label="Bayesian", ls=:dash)
    plot!(dist_true, a, b; label="true dist", ls=:dot, c=:black)
    h = pdf(dist_pred(x), mode(dist_pred(x)))
    scatter!(x, fill(-0.05h, length(x));
        label="sample", ms=1.5, msc=:auto, alpha=0.5, c=:red)
    title!("prediction distributions")

    plot(P1, P2, P3, P4; size=(800, 500), layout=(2, 2))
    plot!(; plot_title="P-value functions, etc. \
        for size-$m sample of $(dist_true)")
end

# %%
plot_ttest(dist_true=Gamma(10, 1), m=10, seed=4649373)

# %%
plot_ttest(dist_true=Gamma(10, 1), m=20, seed=4649373)

# %%
plot_ttest(dist_true=Gamma(10, 1), m=30, seed=4649373)

# %%
plot_ttest(dist_true=Gamma(10, 1), m=40, seed=4649373)

# %%
plot_ttest(dist_true=Gamma(10, 1), m=80, seed=4649373)

# %%
plot_ttest(dist_true=Gamma(10, 1), m=160, seed=4649373)

# %%
plot_ttest(dist_true=Gamma(10, 1), m=320, seed=4649373)

# %%
plot_ttest(dist_true=Gamma(10, 1), m=640, seed=4649373)

# %%
