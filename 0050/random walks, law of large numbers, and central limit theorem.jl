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
#     display_name: Julia 1.10.4
#     language: julia
#     name: julia-1.10
# ---

# %%
using Distributions
n, p = 10^4, 1/6
bin = Binomial(n, p)
@show a = sum(pdf(bin, k) for k in support(bin) if k ≤ 0.99n*p);
@show b = sum(pdf(bin, k) for k in support(bin) if k ≥ 1.01n*p);
a + b

# %%
using Distributions
n, p = 10^4, 1/6
bin = Binomial(n, p)
@show a = sum(pdf(bin, k) for k in support(bin) if k ≤ 0.97n*p);
@show b = sum(pdf(bin, k) for k in support(bin) if k ≥ 1.03n*p);
a + b

# %%
using Distributions
n, p = 10^4, 1/6
bin = Binomial(n, p)
@show sum(pdf(bin, k) for k in support(bin) if k ≤ 1593);
@show sum(pdf(bin, k) for k in support(bin) if k ≥ 1741);

# %%
using Distributions
n, p = 10^4, 1/6
bin1 = Binomial(n, 0.99p)
bin2 = Binomial(n, 1.01p)
@show sum(pdf(bin1, k) for k in support(bin) if k ≤ 1593);
@show sum(pdf(bin2, k) for k in support(bin) if k ≥ 1741);

# %%
using Distributions
n, p = 10^4, 1/6
bin1 = Binomial(n, 0.93p)
bin2 = Binomial(n, 1.07p)
@show sum(pdf(bin1, k) for k in support(bin) if k ≤ 1593);
@show sum(pdf(bin2, k) for k in support(bin) if k ≥ 1741);

# %%
n/6-2std(bin)

# %%
using Distributions
using Random
using StatsPlots
default(fmt=:png, legend=false, tickfontsize=6, plot_titlefontsize=8)

function plot_phat(; n=10^4, p=1/6, seed=nothing,
        m=n-1000, relerr=0.03, ylim=((1-relerr)*p, (1+relerr)*p),
        size=(600, 90), kwargs...)
    
    if seed isa Integer
        Random.seed!(seed)
    end
    
    ber = Bernoulli(p)
    X = rand(ber, n)
    S = cumsum(X)
    phat = S ./ (1:n)

    P1 = plot(1:n, phat)
    hline!([p])
    plot!(; ylim)
    
    P2 = plot(m:n, phat[m:n])
    hline!([p])
    plot!(; ylim)

    plot(P1, P2; size, layout=(1, 2))
    if seed isa Integer
        plot!(plot_title="seed = $seed")
    end
    plot!(; kwargs...)
end

for s in 1:50
    plot_phat(seed=s) |> display
end

# %%
p = 1/6
for s in (13, 15, 17, 19, 20, 21, 27, 32, 35, 37, 41)
    plot_phat(seed=s, relerr=0.1) |> display
end

# %%
for s in 1:50
    plot_phat(seed=s, relerr=0.07) |> display
end

# %%
using Distributions
using Random: rand!

function sim_dice(; n=10^4, relerr=0.05, L=10^7)
    p = 1/6
    dice = Multinomial(n, fill(p, 6))
    nth = Threads.nthreads()
    Xtmp = [rand(dice) for _ in 1:nth]
    c = zeros(Int, nth)
    Threads.@threads :static for i in 1:L
        tid = Threads.threadid()
        X = rand!(dice, Xtmp[tid])
        c[tid] += maximum(k->abs(X[k]/(n*p) - 1), 1:6) ≥ relerr
    end
    sum(c)/L
end

for relerr in 0.01:0.01:0.1
    @eval @show round(sim_dice(; relerr=$relerr); digits=3)
end
println()
for relerr in 0.01:0.01:0.1
    @eval @show round(sim_dice(; relerr=$relerr); sigdigits=3)
end

# %%
using Distributions
using Random
using StatsPlots
default(fmt=:png, legend=false, 
    tickfontsize=6, plot_titlefontsize=8, legendfontsize=12)

function plot_randomwalks(; n=10^4, p=1//6, niters=100, seed=nothing, kwargs...)
    if seed isa Integer
        Random.seed!(seed)
    end
    
    ber = Bernoulli(p)
    σ = √(p*(1-p))
    
    plot()
    for i in 1:niters
        X = rand(ber, n)
        Y = cumsum(X) .- (1:n) .* p
        plot!(0:n, [0; Y]; lw=0.1, a=0.3, label="")
    end
    if seed isa Integer
        plot!(plot_title="seed = $seed")
    end
    plot!(0:n, x-> 2σ√x; label="±2√(np(1-p))", c=:red)
    plot!(0:n, x->-2σ√x; label="", c=:red)
    plot!(xguide="n", yguide="X₁ + ⋯ + Xₙ - np")
    plot!(legend=true)
    title!("i.i.d. Xᵢ ~ Bernoulli(p),  p = $p")
    plot!(; kwargs...)
end

Random.seed!(4649373)
plot_randomwalks()

# %%
using Distributions
using Random
using StatsPlots
default(fmt=:png, legend=false, 
    tickfontsize=6, plot_titlefontsize=8, legendfontsize=12)

function plot_lln(; n=10^4, p=1//6, niters=100, seed=nothing, kwargs...)
    if seed isa Integer
        Random.seed!(seed)
    end
    
    ber = Bernoulli(p)
    σ = √(p*(1-p))
    
    plot()
    for i in 1:niters
        X = rand(ber, n)
        Y = Y = (cumsum(X) .- (1:n) .* p) ./ (1:n)
        plot!(0:n, [0; Y]; lw=0.1, a=0.3, label="")
    end
    if seed isa Integer
        plot!(plot_title="seed = $seed")
    end
    plot!(0:n, x->2σ/√x; label="±2√(p(1-p))/√n", c=:red)
    plot!(0:n, x->-2σ/√x; label="", c=:red)
    plot!(xguide="n", yguide="(X₁ + ⋯ + Xₙ - np) / n")
    plot!(ylim=(-σ/3, σ/3))
    plot!(legend=true)
    title!("i.i.d. Xᵢ ~ Bernoulli(p),  p = $p")
    plot!(; kwargs...)
end

Random.seed!(4649373)
plot_lln()

# %%
using Distributions
using Random
using StatsPlots
default(fmt=:png, legend=false, 
    tickfontsize=6, plot_titlefontsize=8, legendfontsize=12)

function plot_clt(; n=10^4, p=1//6, niters=100, L=10^5, seed=nothing, kwargs...)
    if seed isa Integer
        Random.seed!(seed)
    end
    
    ber = Bernoulli(p)
    σ = √(p*(1-p))
    
    P1 = plot()
    for i in 1:niters
        X = rand(ber, n)
        Y = (cumsum(X) .- (1:n) .* p) ./ .√(1:n)
        plot!(0:n, [0; Y]; lw=0.1, a=0.3, label="")
    end
    if seed isa Integer
        plot!(plot_title="seed = $seed")
    end
    plot!(0:n, x-> 2σ; label="±2√(p(1-p))", c=:red)
    plot!(0:n, x->-2σ; label="", c=:red)
    plot!(xguide="n", yguide="(X₁ + ⋯ + Xₙ - np) / √n")
    plot!(legend=true)
    title!("i.i.d. Xᵢ ~ Bernoulli(p),  p = $p")
    plot!(; kwargs...)
    
    Ylast = (rand(Binomial(n, float(p)), L) .- n*p) ./ √n
    P2 = histogram(Ylast; norm=true, alpha=0.3, label="")
    plot!(Normal(0, σ); label="Normal(0, √(p(1-p)))", lw=2)
    plot!(xguide="(X₁ + ⋯ + Xₙ - np) / √n")
    plot!(legend=true, legendfontsize=10)
    title!("i.i.d. Xᵢ ~ Bernoulli(p),  p = $p,  n = $n")
    plot!(; kwargs...)
    
    plot(P1, P2; size=(600, 700), layout=(2, 1))
end

Random.seed!(4649373)
plot_clt()

# %%
