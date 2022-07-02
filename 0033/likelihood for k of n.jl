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
using Base.Threads
using Distributions
using StatsPlots
default(fmt=:png, titlefontsize=12)

# %%
function sim_binomial(n, k, ps; L=10^5)
    lik = similar(ps)
    @threads for i in eachindex(ps)
        bin = Binomial(n, ps[i])
        lik[i] = mean(rand(bin) == k for _ in 1:L) # likelihood
    end
    lik
end

function plot_sim_binomial(n, k, ps; L=10^5, kwargs...)
    @time lik = sim_binomial(n, k, ps; L)
    plot(ps, lik; label="", kwargs...)
end

# %%
plot_sim_binomial(1000, 967, 0.935:0.0001:0.985; L=10^6)

# %%
plot_sim_binomial(1000, 990, 0.97:0.0001:1; L=10^6)

# %%
plot_sim_binomial(1000, 995, 0.98:0.0001:1; L=10^6)

# %%
plot(p -> pdf(Binomial(1000, p), 995), 0.98, 1; label="")

# %%
plot(p -> pdf(Binomial(1000, p), 996), 0.98, 1; label="")

# %%
plot(p -> pdf(Binomial(1000, p), 997), 0.98, 1; label="")

# %%
plot(p -> pdf(Binomial(1000, p), 998), 0.98, 1; label="")

# %%
plot(p -> pdf(Binomial(1000, p), 999), 0.98, 1; label="")

# %%
plot(p -> pdf(Binomial(1000, p), 1000), 0.98, 1; label="")

# %%
