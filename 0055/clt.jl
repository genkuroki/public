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
#     display_name: Julia
#     language: julia
#     name: julia
# ---

# %%
using Random
using Distributions
using StatsPlots
default(fmt=:png, legend=false, size=(320, 200), titlefontsize=10)

# %%
dist = Uniform()
@show skewness(dist), kurtosis(dist)
plot(x -> pdf(dist, x), -0.2, 1.2; title="Uniform()") |> display

PP = []
for n in 1:6
    X = zeros(n)
    Xbar = [mean(rand!(dist, X)) for _ in 1:10^6]
    P = density(Xbar; title="sample size = $n")
    plot!(Normal(mean(dist), std(dist)/sqrt(n)); ls=:dot)
    push!(PP, P)
end
plot(PP...; size=(900, 400), layout=(2, 3))

# %%
dist = Logistic()
@show skewness(dist), kurtosis(dist)
plot(dist; title="Logistic()") |> display

PP = []
for n in 1:6
    X = zeros(n)
    Xbar = [mean(rand!(dist, X)) for _ in 1:10^6]
    P = density(Xbar; title="sample size = $n")
    plot!(Normal(mean(dist), std(dist)/sqrt(n)); ls=:dot)
    push!(PP, P)
end
plot(PP...; size=(900, 400), layout=(2, 3))

# %%
dist = Laplace()
@show skewness(dist), kurtosis(dist)
plot(dist; title="Laplace()") |> display

PP = []
for n in 1:6
    X = zeros(n)
    Xbar = [mean(rand!(dist, X)) for _ in 1:10^6]
    P = density(Xbar; title="sample size = $n")
    plot!(Normal(mean(dist), std(dist)/sqrt(n)); ls=:dot)
    push!(PP, P)
end
plot(PP...; size=(900, 400), layout=(2, 3))

# %%
dist = TDist(3+eps())
@show skewness(dist), kurtosis(dist)
plot(dist, -6, 6; title="TDist(3)") |> display

PP = []
for n in 1:6
    X = zeros(n)
    Xbar = [mean(rand!(dist, X)) for _ in 1:10^6]
    Xbar = Xbar[abs.(Xbar) .< 8/sqrt(n)]
    P = density(Xbar; title="sample size = $n")
    plot!(Normal(mean(dist), std(dist)/sqrt(n)); ls=:dot)
    push!(PP, P)
end
plot(PP...; size=(900, 400), layout=(2, 3))

# %%
dist = TDist(3.01)
@show skewness(dist), kurtosis(dist)
plot(dist, -6, 6; title="TDist(3.01)") |> display

PP = []
for n in (10, 20, 40, 80, 160, 320)
    X = zeros(n)
    Xbar = [mean(rand!(dist, X)) for _ in 1:10^6]
    Xbar = Xbar[abs.(Xbar) .< 8/sqrt(n)]
    P = density(Xbar; title="sample size = $n")
    plot!(Normal(mean(dist), std(dist)/sqrt(n)); ls=:dot)
    push!(PP, P)
end
plot(PP...; size=(900, 400), layout=(2, 3))

# %%
dist = TDist(4.01)
@show skewness(dist), kurtosis(dist)
plot(dist, -6, 6; title="TDist(4.01)") |> display

PP = []
for n in 1:6
    X = zeros(n)
    Xbar = [mean(rand!(dist, X)) for _ in 1:10^6]
    Xbar = Xbar[abs.(Xbar) .< 8/sqrt(n)]
    P = density(Xbar; title="sample size = $n")
    plot!(Normal(mean(dist), std(dist)/sqrt(n)); ls=:dot)
    push!(PP, P)
end
plot(PP...; size=(900, 400), layout=(2, 3))

# %%
dist = TDist(4.01)
@show skewness(dist), kurtosis(dist)
plot(dist, -6, 6; title="TDist(4.01)") |> display

PP = []
for n in (10, 20, 30, 40, 50, 100)
    X = zeros(n)
    Xbar = [mean(rand!(dist, X)) for _ in 1:10^6]
    Xbar = Xbar[abs.(Xbar) .< 8/sqrt(n)]
    P = density(Xbar; title="sample size = $n")
    plot!(Normal(mean(dist), std(dist)/sqrt(n)); ls=:dot)
    push!(PP, P)
end
plot(PP...; size=(900, 400), layout=(2, 3))

# %%
dist = Exponential()
@show skewness(dist), kurtosis(dist)
plot(x -> pdf(dist, x), -0.2, 6.2; title="Exponential()") |> display

PP = []
for n in 1:6
    X = zeros(n)
    Xbar = [mean(rand!(dist, X)) for _ in 1:10^6]
    P = density(Xbar; title="sample size = $n")
    plot!(Normal(mean(dist), std(dist)/sqrt(n)); ls=:dot)
    push!(PP, P)
end
plot(PP...; size=(900, 400), layout=(2, 3))

# %%
dist = Exponential()
@show skewness(dist), kurtosis(dist)
plot(x -> pdf(dist, x), -0.2, 6.2; title="Exponential()") |> display

PP = []
for n in (10, 20, 30, 40, 50, 100)
    X = zeros(n)
    Xbar = [mean(rand!(dist, X)) for _ in 1:10^6]
    P = density(Xbar; title="sample size = $n")
    plot!(Normal(mean(dist), std(dist)/sqrt(n)); ls=:dot)
    push!(PP, P)
end
plot(PP...; size=(900, 400), layout=(2, 3))

# %%
dist = InverseGamma(3.01)
@show skewness(dist), kurtosis(dist)
plot(range(-0.2, 3.2, 300), x -> x ≤ 0 ? 0.0 : pdf(dist, x); 
    title="InverseGamma(3.01)") |> display

PP = []
for n in (10, 20, 30, 40, 50, 100)
    X = zeros(n)
    Xbar = [mean(rand!(dist, X)) for _ in 1:10^6]
    Xbar = Xbar[Xbar .< mean(dist) + 3/sqrt(n)]
    #@show n, length(Xbar)
    P = density(Xbar; title="sample size = $n")
    plot!(Normal(mean(dist), std(dist)/sqrt(n)); ls=:dot)
    push!(PP, P)
end
plot(PP...; size=(900, 400), layout=(2, 3))

# %%
dist = InverseGamma(3.01)
@show skewness(dist), kurtosis(dist)
plot(range(-0.2, 3.2, 300), x -> x ≤ 0 ? 0.0 : pdf(dist, x); 
    title="InverseGamma(3.01)") |> display

PP = []
@time for n in (50, 100, 200, 400, 800, 1600)
    X = zeros(n)
    Xbar = [mean(rand!(dist, X)) for _ in 1:10^6]
    Xbar = Xbar[Xbar .< mean(dist) + 3/sqrt(n)]
    #@show n, length(Xbar)
    P = density(Xbar; title="sample size = $n")
    plot!(Normal(mean(dist), std(dist)/sqrt(n)); ls=:dot)
    push!(PP, P)
end
plot(PP...; size=(900, 400), layout=(2, 3))

# %%
dist = InverseGamma(4.1)
@show skewness(dist), kurtosis(dist)
plot(range(-0.2, 1.5, 300), x -> x ≤ 0 ? 0.0 : pdf(dist, x); 
    title="InverseGamma(4.1)") |> display

PP = []
@time for n in (10, 20, 40, 80, 160, 320)
    X = zeros(n)
    Xbar = [mean(rand!(dist, X)) for _ in 1:10^6]
    Xbar = Xbar[Xbar .< mean(dist) + 1/sqrt(n)]
    #@show n, length(Xbar)
    P = density(Xbar; title="sample size = $n")
    plot!(Normal(mean(dist), std(dist)/sqrt(n)); ls=:dot)
    push!(PP, P)
end
plot(PP...; size=(900, 400), layout=(2, 3))

# %%
