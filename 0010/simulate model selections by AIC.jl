# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:hydrogen
#     text_representation:
#       extension: .jl
#       format_name: hydrogen
#       format_version: '1.3'
#       jupytext_version: 1.11.2
#   kernelspec:
#     display_name: Julia 1.8.0-DEV
#     language: julia
#     name: julia-1.8
# ---

# %% tags=[]
using Distributions, StatsPlots

function aic(model, Y)
    mle = fit_mle(model, Y)
    -2loglikelihood(mle, Y) + 2length(params(mle))
end

function simulate_model_selections(models, truedist, samplesize; niters = 10^4)
    selectedmodel = Vector{Int}(undef, niters)
    Threads.@threads for i in 1:niters
        Y = rand(truedist, samplesize)
        selectedmodel[i] = argmin(aic(model, Y) for model in models)
    end
    nselected = zeros(Int, 3)
    for i in 1:niters
        nselected[selectedmodel[i]] += 1
    end
    [model => nselected[i]/niters for (i, model) in enumerate(models)]
end

models = (Gamma, Weibull, LogNormal)

a, b = 0, 40
plot(Gamma(5, 2), a, b; label="Gamma(5, 2)")
plot!(Weibull(2.4, 11.3), a, b; label="Weibull(2.4, 11.3)", ls=:dash)
plot!(LogNormal(2.2, 0.47), a, b; label="LogNormal(2.2, 0.47)", ls=:dashdot)

# %%
simulate_model_selections(models, Gamma(5, 2), 100; niters = 10^4)

# %%
simulate_model_selections(models, Weibull(2.4, 11.3), 100; niters = 10^4)

# %%
simulate_model_selections(models, LogNormal(2.2, 0.47), 100; niters = 10^4)

# %%
