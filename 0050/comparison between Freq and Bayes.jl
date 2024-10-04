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
#     display_name: Julia 1.10.5
#     language: julia
#     name: julia-1.10
# ---

# %%
ENV["COLUMNS"] = 200
using Distributions
using Optim: optimize, Brent
using Turing
using StatsPlots
default(fmt=:png, titlefontsize=10)

# %% [markdown]
# ## 二項分布モデルの場合

# %% [markdown]
# ### 無情報事前分布の場合

# %%
function confint_score(k_new, n_new, α=0.05; prior_data=(k_prior=0, n_prior=0))
    (; k_prior, n_prior) = prior_data
    k, n = k_prior + k_new, n_prior + n_new
    p̂ = k/n
    z = quantile(Normal(), 1-α/2)
    a, b, c = 1+z^2/n, p̂+z^2/(2n), p̂^2
    # ap² - 2bp + c = 0 を解く.
    sqrtD = √(b^2 - a*c)
    p_L = (b - sqrtD)/a
    p_U = (b + sqrtD)/a
    [p_L, p_U]
end

sigdigits::Int = 3
r(x) = round(x; sigdigits=sigdigits)

function highest_density_interval(dist, α=0.05)
    f(p) = quantile(dist, min(1, p+(1-α))) - quantile(dist, p)
    o = optimize(f, 0, α, Brent())
    p = o.minimizer
    [quantile(dist, p), quantile(dist, min(1, p+(1-α)))]
end

function credint_hdi(k, n, α=0.05; prior=Beta(1, 1))
    κ, λ = params(prior)
    posterior = Beta(κ+k, λ+n-k)
    m = mode(posterior)
    α ≈ 1 && return [m, m]
    highest_density_interval(posterior, α)
end

function plot_ci(; k_new=7, n_new=24, α=0.05,
        prior=Beta(1, 1), prior_data=(k_prior=0, n_prior=0),
        f=Bool[1,1], s=Bool[0,0], c=Bool[1,1], ls2=:solid, kwargs...)
    κ, λ = params(prior)
    @show k_new n_new α prior
    (; k_prior, n_prior) = prior_data
    @show prior_data
    println()
    k, n = k_prior + k_new, n_prior + n_new
    posterior = Beta(κ + k_new, λ + n_new - k_new)
    pointest_score = k/n |> r
    pointest_bayes = mode(posterior) |> r
    @show pointest_score pointest_bayes
    ci_score = confint_score(k_new, n_new, α; prior_data) .|> r
    ci_bayes = credint_hdi(k_new, n_new, α; prior) .|> r
    @show ci_score ci_bayes
    println()
    
    if s[1] | s[2]
        if s[1] & s[2]
            title = "CIs:  k_new=$k_new,  n_new=$n_new"
        elseif s[1]
            title = "confidence intervals:  k_new=$k_new,  n_new=$n_new"
        else
            title = "credible intervals:  k_new=$k_new,  n_new=$n_new"
        end
        yguide = "α"
    else
        title = "P-value function" * (all(f) ? "s" : "") * ":  k_new=$k_new,  n_new=$n_new"
        yguide = "P-value"
    end
    title *= "\nprior=Beta$(params(prior)),  prior_data=$prior_data"

    plot(; title)
    if f[1]
        if s[1]
            l = true
            for α in [0.001; 0.01:0.01:1]
                plot!(confint_score(k_new, n_new, α; prior_data), fill(α, 2);
                    label=l ? "score method" : "", c=1)
                l = false
            end
        else
            αs = [0.001:0.001:0.009; 0.01:0.01:1]
            scatter!(vcat(confint_score.(k_new, n_new, αs; prior_data)...), repeat(αs; inner=2); 
                label="score method", ms=1.5, msc=:auto, c=1)
            c[1] && plot!(ci_score, fill(1.1α, 2); label="$(100(1-α))% score conf.int.", c=1, lw=2)
        end
    end
    if f[2]
        if s[2]
            l = true
            for α in [0.001; 0.01:0.01:1]
                plot!(credint_hdi.(k_new, n_new, α; prior), fill(α, 2);
                    label=l ? "Bayes HDI" : "", c=2, ls=ls2)
                l = false
            end
        else
            αs = [0.001:0.001:0.009; 0.01:0.01:1]
            scatter!(vcat(credint_hdi.(k_new, n_new, αs; prior)...), repeat(αs; inner=2); 
                label="Bayes HDI", ms=1.5, msc=:auto, c=2)
            c[2] && plot!(ci_bayes, fill(0.9α, 2); label="$(100(1-α))% cred.int. HDI", c=2, lw=2)
        end
    end
    plot!(; xguide="p", yguide)
    plot!(; ytick=0:0.05:1, ylim=(-0.03, 1.05))
    plot!(; kwargs...)
end

plot_ci(; k_new=6, n_new=20)

# %%
plot_ci(; k_new=6, n_new=20, f=Bool[1,0], s=Bool[1,0])

# %%
plot_ci(; k_new=6, n_new=20, f=Bool[0,1], s=Bool[0,1])

# %%
plot_ci(; k_new=6, n_new=20, f=Bool[1,1], s=Bool[1,0], c=Bool[1,0])

# %%
plot_ci(; k_new=24, n_new=80)

# %%
plot_ci(; k_new=96, n_new=320)

# %% [markdown]
# ## 偏りのある事前分布の場合

# %%
prior = Beta(6, 4)
plot(prior; label="prior", title="prior = Beta$(params(prior))", c=2)
plot!(xguide="p", yguide="probability density")
plot!(xtick=0:0.1:1)

# %%
plot_ci(; k_new=6, n_new=20, prior)

# %%
@show κ, λ = params(prior);
@show prior_data = (k_prior=κ-1, n_prior=κ+λ-2);

# %%
plot_ci(; k_new=6, n_new=20, prior, prior_data)

# %%
plot_ci(; k_new=24, n_new=80, prior)

# %%
plot_ci(; k_new=24, n_new=80, prior, prior_data)

# %%
plot_ci(; k_new=96, n_new=320, prior=Beta(3, 2))

# %%
plot_ci(; k_new=96, n_new=320, prior=Beta(3, 2), prior_data=(k_prior=2, n_prior=3))

# %% [markdown]
# ## 正規分布モデルの場合

# %%
r1(x) = round(x; digits=1)

function confintmean(n, x̄, s², α=0.05)
    c = cquantile(Normal(), α/2)
    sehat = √(s²/n)
    [x̄ - c*sehat, x̄ + c*sehat]
end

function confintmean(X, α=0.05)
    confintmean(length(X), mean(X), var(X), α)
end

@model function model_normal(X; prior_μ=Flat(), prior_logσ=Flat())
    μ ~ prior_μ
    logσ ~ prior_logσ
    σ = exp(logσ)
    for i in eachindex(X)
        X[i] ~ Normal(μ, σ)
    end
end

function credint(A, α=0.05; sorted=true, α_m=0.95)
    α == 0 && return [-Inf, Inf]
    if α == 1
        m = mean(credint(A, α_m; sorted))
        return [m, m]
    end
    θ = sorted ? A : sort(A)
    n = length(θ)
    m = floor(Int, α*n)
    f(i) = θ[i+n-m-1] - θ[i]
    i = argmin(f, 1:m)
    [θ[i], θ[i+n-m-1]]
end

function plot_ci_mean(X; α=0.05,
        prior_μ=Flat(), prior_logσ=Flat(), N=10^5, nchains=10)
    chain = sample(model_normal(X; prior_μ, prior_logσ), NUTS(), MCMCThreads(), N, nchains)
    display(chain)
    println()
    #plot(chain; lw=0.2)
    println()
    
    μ = sort(vec(chain[:μ]))
    @show n = length(X)
    @show x̄ = mean(X)
    @show s = std(X)
    @show prior_μ prior_logσ
    @show α
    println()
    
    pointest_ttest = x̄ |> r1
    pointest_bayes = mean(credint(μ, 1.0)) |> r1
    @show pointest_ttest pointest_bayes
    ci_ttest = confintmean(X, α) .|> r1
    ci_bayes = credint(μ, α) .|> r1
    @show ci_ttest ci_bayes
    println()
    
    αs = [0.001:0.001:0.009; 0.01:0.01:0.09; 0.1:0.02:1]
    plot(title="P-value functions")
    scatter!(vcat(confintmean.(n, x̄, s^2, αs)...), repeat(αs; inner=2); label="1-sample t-test", ms=1.5, msc=:auto, c=1)
    plot!(ci_ttest, fill(1.1α, 2); label="$(100(1-α))% conf.int.", c=1, lw=2)
    scatter!(vcat(credint.((μ,), αs)...), repeat(αs; inner=2); label="Bayes HDI", ms=1.5, msc=:auto, c=2)
    plot!(ci_bayes, fill(0.9α, 2); label="$(100(1-α))% cred.int. HDI", c=2, lw=2)
    plot!(xguide="x", yguide="P-value")
    plot!(ytick=0:0.05:1)
end

# %%
X = [
    143.7, 149.9, 151.7, 155.1, 155.4, 156.4, 157.1, 157.2, 157.7, 157.8,
    158.1, 158.3, 158.7, 159.2, 159.3, 159.8, 159.9, 160.0, 160.8, 161.1, 
    161.2, 161.9, 162.1, 162.2, 162.5, 163.1, 163.1, 164.2, 165.3, 165.3, 
    165.4, 165.5, 165.7, 166.3, 166.3, 166.7, 166.8, 166.9, 166.9, 167.0, 
    167.0, 167.1, 167.6, 167.8, 167.8, 167.9, 168.0, 168.8, 169.0, 169.1, 
    169.4, 169.6, 170.5, 170.6, 171.5, 172.2, 172.7, 173.0, 173.3, 173.7, 
    173.8, 173.9, 174.0, 175.0, 175.1, 175.7, 175.7, 175.8, 176.0, 176.2, 
    176.6, 177.1, 177.2, 177.3, 177.4, 177.5, 177.8, 177.8, 178.1, 179.2, 
    179.3, 179.6, 180.0, 180.1, 181.4, 181.6, 181.7, 182.3, 182.5, 183.3, 
    183.4, 183.4, 183.6, 183.7, 185.1, 185.8, 186.3, 186.7, 188.0, 195.5
]

N = 10^4
nchains=10
chain = sample(model_normal(X), NUTS(), MCMCThreads(), N, nchains)
println()

μ = sort(vec(chain[:μ]))

@show length(X) mean(X) std(X) log(std(X))
println()
@show confintmean(X) .|> r1
@show credint(μ) .|> r1
println()

display(chain)
println()

plot(chain; lw=0.2)

# %%
X = [
    143.7, 149.9, 151.7, 155.1, 155.4, 156.4, 157.1, 157.2, 157.7, 157.8,
    158.1, 158.3, 158.7, 159.2, 159.3, 159.8, 159.9, 160.0, 160.8, 161.1, 
    161.2, 161.9, 162.1, 162.2, 162.5, 163.1, 163.1, 164.2, 165.3, 165.3, 
    165.4, 165.5, 165.7, 166.3, 166.3, 166.7, 166.8, 166.9, 166.9, 167.0, 
    167.0, 167.1, 167.6, 167.8, 167.8, 167.9, 168.0, 168.8, 169.0, 169.1, 
    169.4, 169.6, 170.5, 170.6, 171.5, 172.2, 172.7, 173.0, 173.3, 173.7, 
    173.8, 173.9, 174.0, 175.0, 175.1, 175.7, 175.7, 175.8, 176.0, 176.2, 
    176.6, 177.1, 177.2, 177.3, 177.4, 177.5, 177.8, 177.8, 178.1, 179.2, 
    179.3, 179.6, 180.0, 180.1, 181.4, 181.6, 181.7, 182.3, 182.5, 183.3, 
    183.4, 183.4, 183.6, 183.7, 185.1, 185.8, 186.3, 186.7, 188.0, 195.5
]

plot_ci_mean(X)

# %%
plot_ci_mean(X; prior_μ=Normal(200, 10), prior_logσ=Normal(0, 10))

# %%
plot(Normal(200, 10), 160, 240; label="", title="Normal(200, 10)", xtick=160:5:240, c=2)

# %%
plot_ci_mean(X; prior_μ=Normal(200, 5), prior_logσ=Normal(0, 10))

# %%
plot(Normal(200, 5), 160, 240; label="", title="Normal(200, 5)", xtick=160:5:240, c=2)

# %%
