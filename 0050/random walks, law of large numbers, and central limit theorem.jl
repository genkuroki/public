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

function sim_dice_relerr(; n=10^4, p=1/6, relerr=0.05, L=10^7)
    dice = Multinomial(n, fill(p, 6))
    nth = Threads.nthreads()
    Xtmp = [rand(dice) for _ in 1:nth]
    c = zeros(Int, nth)
    Threads.@threads :static for i in 1:L
        tid = Threads.threadid()
        X = rand!(dice, Xtmp[tid])
        c[tid] += any(k-> abs(X[k]/(n*p) - 1) > relerr, 1:6)
    end
    sum(c)/L
end

n = 10^4
println("(X₁,…,X₆) ~ Multinomial(n=$n, [1/6, 1/6, 1/6, 1/6, 1/6, 1/6])")
for relerr in 0.01:0.01:0.1
    prob = sim_dice_relerr(; n, relerr)
    prob1 = round(prob; digits=3)
    prob2 = round(prob; sigdigits=3)
    print("(probability that |Xₖ/(np) - 1| > $relerr for some k) = ")
    if prob1 == prob2
        println(prob1)
    else
        println(prob1, " (= ", prob2, ")")
    end
end

# %%
using Distributions
using Random: rand!

function sim_dice_abserr(; n=10^4, p=1/6, abserr=0.01, L=10^7)
    dice = Multinomial(n, fill(p, 6))
    nth = Threads.nthreads()
    Xtmp = [rand(dice) for _ in 1:nth]
    c = zeros(Int, nth)
    Threads.@threads :static for i in 1:L
        tid = Threads.threadid()
        X = rand!(dice, Xtmp[tid])
        c[tid] += any(k-> abs(X[k]/n - p) > abserr, 1:6)
    end
    sum(c)/L
end

n = 10000
println("(X₁,…,X₆) ~ Multinomial(n=$n, [1/6, 1/6, 1/6, 1/6, 1/6, 1/6])")
for abserr in 0.001:0.001:0.016
    prob = sim_dice_abserr(; n, abserr)
    prob1 = round(prob; digits=3)
    prob2 = round(prob; sigdigits=3)
    print("(probability that |Xₖ/n - p| > $abserr for some k) = ")
    if prob1 == prob2
        println(prob1)
    else
        println(prob1, " (= ", prob2, ")")
    end
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

function plot_lln(; n=10^4, p=1//6, niters=100, seed=nothing,
        σ=√(p*(1-p)), ylim=(-σ/4, σ/4), kwargs...)
    if seed isa Integer
        Random.seed!(seed)
    end
    
    ber = Bernoulli(p)
    
    plot()
    for i in 1:niters
        X = rand(ber, n)
        Y = Y = (cumsum(X) .- (1:n) .* p) ./ (1:n)
        plot!(0:n, [0; Y]; lw=0.1, a=0.3, label="")
    end
    if seed isa Integer
        plot!(plot_title="seed = $seed")
    end
    plot!(0:n, x->2σ/√x; label="±2√(p(1-p)/n)", c=:red)
    plot!(0:n, x->-2σ/√x; label="", c=:red)
    plot!(xguide="n", yguide="(X₁ + ⋯ + Xₙ)/n - p")
    plot!(legend=true)
    hline!([0.0]; label="", c=:black, lw=0.5)
    title!("i.i.d. Xᵢ ~ Bernoulli(p),  p = $p")
    plot!(; ylim, kwargs...)
end

Random.seed!(4649373)
plot_lln(; xtick=0:1000:10000, ytick=-0.1:0.01:0.1)

# %%
Random.seed!(4649373)
plot_lln(n=1000, p=1/4, ylim=(-0.16, 0.16))

ns = [10:10:50; 100:100:1000]
ks = [2, 6, 12, 15, 20, 36, 62, 83, 109, 130, 153, 170, 201, 223, 249]
ys = @. ks/ns - 0.25
plot!(ns, ys; label="data", c=:blue)
plot!(xtick=0:100:1000)

# %%
Random.seed!(4649373)
plot_lln(n=1000, p=1/2, ylim=(-0.16, 0.16))

ns = [
     10,  20,  30,  40,  50,  60,  70,  80,  90,  100,
    150, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
]
ks = [
      6,  10,  18,  23,  28,  34,  38,  45,  49,   54,
     77,  91, 147, 192, 247, 299, 353, 400, 448,  497,
]
ys = @. ks/ns - 0.5
plot!(ns, ys; label="data", c=:blue)
plot!(xtick=0:100:1000)

# %% [markdown]
# <img src="IMG_5571.jpg" width=610>

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
using Distributions
using DataFrames
using StatsPlots
default(fmt=:png, titlefontsize=10, tickfontsize=6)

df = DataFrame(
    n = [50; 100; 200:200:2000],
    k = [9, 14, 31, 69, 91, 126, 163, 205, 238, 267, 301, 333]
)
@show df

PP = []
p = 1//6
for i in 1:12
    n = df.n[i]
    bin = Binomial(n, p)
    μ, σ = mean(bin), std(bin)
    P = bar(bin; alpha=0.3, lc=:match, label="")
    plot!(xlim=round.(Int, (μ-4σ, μ+4σ)))
    plot!(xtick=round.(Int, μ-4σ:σ:μ+4σ))
    vline!([μ]; label="", c=:black)
    vline!([df.k[i]]; label="data", c=:red)
    title!("Binomial(n=$n, p=$p)")
    push!(PP, P)
end
plot(PP...; size=(1000, 1000), layout=(4, 3))

# %% tags=[]
using Distributions
using DataFrames
using StatsPlots
default(fmt=:png, titlefontsize=10, tickfontsize=6)

PP = []
p = 1//6
for m in 2:7
    n = 10^m
    bin = Binomial(n, p)
    μ, σ = mean(bin), std(bin)
    P = plot()
    if n < 10^4
        bar!(bin; alpha=0.3, lc=:match, label="")
    else
        plot!(Normal(μ, σ), μ-4.2σ, μ+4.2σ; label="")
    end
    plot!(xlim=round.(Int, (μ-4.2σ, μ+4.2σ)))
    plot!(xtick=round.(Int, μ-5σ:σ:μ+5σ))
    vline!([μ]; label="", c=:black)
    title!("Binomial(n=10^$m, p=$p)")
    push!(PP, P)
end
plot(PP...; size=(1000, 600), layout=(3, 2))

# %%
