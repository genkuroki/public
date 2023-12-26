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
#     display_name: Julia 1.9.4
#     language: julia
#     name: julia-1.9
# ---

# %%
using Distributions
using StatsPlots
default(fmt=:png)

# %%
n = 10
p = 5/n
dist = Binomial(n, p)
normal = Normal(mean(dist), std(dist))
bar(dist, 0:min(n, round(Int, 4n*p)); label="Binomial($n, $p)", alpha=0.3)
plot!(x -> pdf(Poisson(n*p), round(Int, x)); label="Poisson($(n*p))", ls=:dash)
plot!(normal; label="normal approx.", lw=1.5)

# %%
n = 40
p = 5/n
dist = Binomial(n, p)
normal = Normal(mean(dist), std(dist))
bar(dist, 0:min(n, round(Int, 4n*p)); label="Binomial($n, $p)", alpha=0.3)
plot!(x -> pdf(Poisson(n*p), round(Int, x)); label="Poisson($(n*p))", ls=:dash)
plot!(normal; label="normal approx.", lw=1.5)

# %%
n = 100
p = 5/n
dist = Binomial(n, p)
normal = Normal(mean(dist), std(dist))
bar(dist, 0:min(n, round(Int, 4n*p)); label="Binomial($n, $p)", alpha=0.3)
plot!(x -> pdf(Poisson(n*p), round(Int, x)); label="Poisson($(n*p))", ls=:dash)
plot!(normal; label="normal approx.", lw=1.5)

# %%
n = 400
p = 5/n
dist = Binomial(n, p)
normal = Normal(mean(dist), std(dist))
bar(dist, 0:min(n, round(Int, 4n*p)); label="Binomial($n, $p)", alpha=0.3)
plot!(x -> pdf(Poisson(n*p), round(Int, x)); label="Poisson($(n*p))", ls=:dash)
plot!(normal; label="normal approx.", lw=1.5)

# %%
n = 400
p = 4/n
dist = Binomial(n, p)
normal = Normal(mean(dist), std(dist))
bar(dist, 0:min(n, round(Int, 4n*p)); label="Binomial($n, $p)", alpha=0.3)
plot!(x -> pdf(Poisson(n*p), round(Int, x)); label="Poisson($(n*p))", ls=:dash)
plot!(normal; label="normal approx.", lw=1.5)

# %%
n = 400
p = 3/n
dist = Binomial(n, p)
normal = Normal(mean(dist), std(dist))
bar(dist, 0:min(n, round(Int, 4n*p)); label="Binomial($n, $p)", alpha=0.3)
plot!(x -> pdf(Poisson(n*p), round(Int, x)); label="Poisson($(n*p))", ls=:dash)
plot!(normal; label="normal approx.", lw=1.5)

# %%
n = 400
p = 2/n
dist = Binomial(n, p)
normal = Normal(mean(dist), std(dist))
bar(dist, 0:min(n, round(Int, 4.5n*p)); label="Binomial($n, $p)", alpha=0.3)
plot!(x -> pdf(Poisson(n*p), round(Int, x)); label="Poisson($(n*p))", ls=:dash)
plot!(normal; label="normal approx.", lw=1.5)

# %%
n = 400
p = 1/n
dist = Binomial(n, p)
normal = Normal(mean(dist), std(dist))
bar(dist, 0:min(n, round(Int, 5n*p)); label="Binomial($n, $p)", alpha=0.3)
plot!(x -> pdf(Poisson(n*p), round(Int, x)); label="Poisson($(n*p))", ls=:dash)
plot!(normal; label="normal approx.", lw=1.5)

# %%
