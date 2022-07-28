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
using Base.Threads
using Distributions
using HypothesisTests
using QuadGK
using Random
using Roots
using StatsBase
using StatsFuns
using StatsPlots
default(fmt=:png, titlefontsize=10, guidefontsize=8, tickfontsize=6)

# %%
function prob_x_le_y(distx, disty, a=0.0)
    H(y) = cdf(distx, y) * pdf(disty, y-a)
    quadgk(H, extrema(disty + a)...)[1]
end

function tieshift(distx, disty; probtie=0.5)
    #s = max(std(distx), std(disty))
    #m = median(distx) - median(disty)
    #find_zero(a -> prob_x_le_y(distx, disty, a) - probtie,
    #    (amin, amax), Bisection())
    find_zero(a -> prob_x_le_y(distx, disty, a) - probtie, 0.0)
end

@show tieshift(Normal(0, 1), Normal(2, 2))
@show tieshift(Normal(0, 1), Laplace(2, 2))
@show tieshift(Normal(0, 1), Uniform(0, 1));

# %%
distx, disty = Gamma(6, 1), Gamma(2, 3)
@show median(distx), median(disty)
@show median(distx) - median(disty)
@show tieshift(distx, disty);

# %%
distx, disty = Uniform(), Chisq(1)
@show median(distx), median(disty)
@show median(distx) - median(disty)
@show tieshift(distx, disty);

# %%
function sim(TestFunc = MannWhitneyUTest;
        distx = Normal(0, 1), disty = Normal(0, 4), m = 100, n = 50,
        L = 10^6)
    pval = Vector{Float64}(undef, L)
    tmpX = [Vector{Float64}(undef, m) for _ in 1:nthreads()]
    tmpY = [Vector{Float64}(undef, n) for _ in 1:nthreads()]
    @threads for i in 1:L
        X = rand!(distx, tmpX[threadid()])
        Y = rand!(disty, tmpY[threadid()])
        pval[i] = pvalue(TestFunc(X, Y))
    end
    pval
    ecdf(pval)
end

distname(dist) = replace(string(dist), r"\{[^\}]*\}"=>"")

function plot_ecdf(ecdf_pval, TestFunc, distx, disty, m, n, a; kwargs...)
    plot(p -> ecdf_pval(p), 0, 0.1; label="ecdf of P-values")
    plot!([0, 0.1], [0, 0.1]; label="", ls=:dot, c=:black)
    plot!(legend=:topleft)
    plot!(xtick=0:0.01:0.1, ytick=0:0.01:1)
    plot!(xguide="nominal significance level α", 
        yguide="probability of P-value < α")
    s = (a < 0 ? "" : "+") * string(round(a; digits=3))
    title!("$TestFunc\n\
        X: $(distname(distx)), m=$m\n\
        Y: $(distname(disty))$s, n=$n")
    plot!(size=(400, 450))
    plot!(; kwargs...)
end

function plot_pvals(
        TestFunc1 = MannWhitneyUTest,
        TestFunc2 = UnequalVarianceTTest;
        distx = Normal(0, 1), disty = Normal(0, 4), m = 100, n = 50,
        L = 10^6, a = nothing, kwargs...)
    
    if isnothing(a)
        @show a = tieshift(distx, disty)
        @show prob_x_le_y(distx, disty + a)
    else
        @show a
        @show median(distx) - median(disty)
    end
    ecdf_pval1 = @time sim(TestFunc1;
        distx = distx,
        disty = disty + a,
        m, n, L, kwargs...)
    ymax1 = ecdf_pval1(0.1)
    @show Δμ = mean(distx) - mean(disty)
    @show mean(distx), mean(disty + Δμ)
    ecdf_pval2 = @time sim(TestFunc2;
        distx = distx,
        disty = disty + Δμ,
        m, n, L, kwargs...)
    ymax2 = ecdf_pval2(0.1)
    ymax = max(ymax1, ymax2)
    P1 = plot_ecdf(ecdf_pval1, TestFunc1, distx, disty, m, n, a;
        ylim=(-0.002, 1.02*ymax), kwargs...)
    P2 = plot_ecdf(ecdf_pval2, TestFunc2, distx, disty, m, n, Δμ;
        ylim=(-0.002, 1.02*ymax), kwargs...)
    plot(P1, P2; size=(800, 450), topmargin=3.5Plots.mm)
end

# %% [markdown]
# ## 正規分布

# %%
plot_pvals(distx = Normal(0, 1), disty = Normal(0, 1), m = 25, n = 25)

# %%
plot_pvals(distx = Normal(0, 1), disty = Normal(0, 1), m = 50, n = 25)

# %%
plot_pvals(distx = Normal(0, 1), disty = Normal(0, 2), m = 25, n = 25)

# %%
plot_pvals(distx = Normal(0, 1), disty = Normal(0, 2), m = 50, n = 25)

# %%
plot_pvals(distx = Normal(0, 1), disty = Normal(0, 1), m = 10, n = 10)

# %%
plot_pvals(distx = Normal(0, 1), disty = Normal(0, 1), m = 20, n = 10)

# %%
plot_pvals(distx = Normal(0, 1), disty = Normal(0, 2), m = 10, n = 10)

# %%
plot_pvals(distx = Normal(0, 1), disty = Normal(0, 2), m = 20, n = 10)

# %% [markdown]
# ## ガンマ分布

# %%
plot_pvals(distx = Gamma(2, 3), disty = Gamma(2, 3), m = 25, n = 25)

# %%
plot_pvals(distx = Gamma(2, 3), disty = Gamma(2, 3), m = 50, n = 25)

# %%
plot_pvals(distx = Gamma(6, 1), disty = Gamma(2, 3), m = 25, n = 25)

# %%
plot_pvals(distx = Gamma(6, 1), disty = Gamma(2, 3), m = 50, n = 25)

# %%
plot_pvals(distx = Gamma(2, 3), disty = Gamma(2, 3), m = 10, n = 10)

# %%
plot_pvals(distx = Gamma(2, 3), disty = Gamma(2, 3), m = 20, n = 10)

# %%
plot_pvals(distx = Gamma(6, 1), disty = Gamma(2, 3), m = 10, n = 10)

# %%
plot_pvals(distx = Gamma(6, 1), disty = Gamma(2, 3), m = 20, n = 10)

# %% [markdown]
# ## 雑多

# %%
plot_pvals(distx = Laplace(0, 1), disty = Laplace(0, 2), m = 10, n = 10)

# %%
plot_pvals(distx = Laplace(0, 1), disty = Laplace(0, 2), m = 25, n = 25)

# %%
plot_pvals(distx = Laplace(0, 1), disty = Laplace(0, 2), m = 50, n = 50)

# %%
plot_pvals(distx = Laplace(0, 1), disty = Laplace(0, 2), m = 20, n = 10)

# %%
plot_pvals(distx = Laplace(0, 1), disty = Laplace(0, 2), m = 50, n = 25)

# %%
plot_pvals(distx = Uniform(-1, 1), disty = Uniform(-2, 2), m = 10, n = 10)

# %%
plot_pvals(distx = Uniform(-1, 1), disty = Uniform(-2, 2), m = 20, n = 10)

# %%
plot_pvals(distx = Uniform(-1, 1), disty = Uniform(-2, 2), m = 50, n = 25)

# %%
@doc LogNormal

# %%
plot(LogNormal(0, 1), 0, 10)
plot!(LogNormal(0, 1.25), 0, 10)

# %%
std(LogNormal(0, 1)), std(LogNormal(0, 1.25))

# %%
plot_pvals(distx = LogNormal(0, 1), disty = LogNormal(0, 1), m = 500, n = 500, L=10^5)

# %%
plot_pvals(distx = LogNormal(0, 1), disty = LogNormal(0, 1.25), m = 500, n = 500, L=10^5)

# %%
plot_pvals(distx = LogNormal(0, 1), disty = LogNormal(0, 1.25), m = 1000, n = 500, L=10^5)

# %%
plot(TDist(2), -10, 10)
plot!(2TDist(2), -10, 10)

# %%
var(TDist(2)), var(TDist(3)), var(TDist(4)), var(TDist(5))

# %%
plot_pvals(distx = TDist(2), disty = TDist(2), m = 100, n = 100)

# %%
plot_pvals(distx = TDist(2), disty = TDist(2), m = 200, n = 100)

# %%
plot_pvals(distx = TDist(2), disty = 2TDist(2), m = 100, n = 100)

# %%
plot_pvals(distx = TDist(2), disty = 2TDist(2), m = 200, n = 100)

# %% [markdown]
# ## median matching vs. tie shifting

# %%
distx, disty = Uniform(-1, 1), Exponential()
m, n, = 100, 100

@show a = tieshift(distx, disty)
ecdf_pval1 = @time sim(; distx = distx, disty = disty + a, m, n)
P1 = plot_ecdf(ecdf_pval1, MannWhitneyUTest, distx, disty, m, n, a)

@show a = median(distx) - median(disty)
ecdf_pval2 = @time sim(; distx = distx, disty = disty + a, m, n)
P2 = plot_ecdf(ecdf_pval2, MannWhitneyUTest, distx, disty, m, n, a)

plot(P1, P2; size=(800, 450), topmargin=4Plots.mm)

# %%
distx, disty = Uniform(-1, 1), Exponential()
@show distx, std(distx)
@show disty, std(disty)

a = @show tieshift(distx, disty)
P1 = plot(distx, -2, 6; label="distx")
plot!(disty + a, -2, 6; label="disty + ($(round(a; digits=3)))", ls=:dash)
title!("case of tie shifting")

a = @show median(distx) - median(disty)
P2 = plot(distx, -2, 6; label="distx")
plot!(disty + a, -2, 6; label="disty + ($(round(a; digits=3)))", ls=:dash)
vline!([median(distx)]; label="median(distx)", ls=:dot, lw=1.5)
title!("case of matching medians")

plot(P1, P2; size=(800, 250))

# %%
distx, disty = Uniform(-1, 1), Exponential(4)
m, n, = 100, 100

@show a = tieshift(distx, disty)
ecdf_pval1 = @time sim(; distx = distx, disty = disty + a, m, n)
P1 = plot_ecdf(ecdf_pval1, MannWhitneyUTest, distx, disty, m, n, a)

@show a = median(distx) - median(disty)
ecdf_pval2 = @time sim(; distx = distx, disty = disty + a, m, n)
P2 = plot_ecdf(ecdf_pval2, MannWhitneyUTest, distx, disty, m, n, a)

plot(P1, P2; size=(800, 450), topmargin=4Plots.mm)

# %%
distx, disty = Uniform(-1, 1), Exponential(4)
@show distx, std(distx)
@show disty, std(disty)

a = @show tieshift(distx, disty)
P1 = plot(distx, -4, 10; label="distx")
plot!(disty + a, -4, 10; label="disty + ($(round(a; digits=3)))", ls=:dash)
title!("case of tie shifting")

a = @show median(distx) - median(disty)
P2 = plot(distx, -4, 10; label="distx")
plot!(disty + a, -4, 10; label="disty + ($(round(a; digits=3)))", ls=:dash)
vline!([median(distx)]; label="median(distx)", ls=:dot, lw=1.5)
title!("case of matching medians")

plot(P1, P2; size=(800, 250))

# %%
distx, disty = Uniform(-1, 1), Exponential(0.5773502691896257)
m, n, = 100, 100

@show a = tieshift(distx, disty)
ecdf_pval1 = @time sim(; distx = distx, disty = disty + a, m, n)
P1 = plot_ecdf(ecdf_pval1, MannWhitneyUTest, distx, disty, m, n, a)

@show a = median(distx) - median(disty)
ecdf_pval2 = @time sim(; distx = distx, disty = disty + a, m, n)
P2 = plot_ecdf(ecdf_pval2, MannWhitneyUTest, distx, disty, m, n, a)

plot(P1, P2; size=(800, 450), topmargin=4Plots.mm)

# %%
distx, disty = Uniform(-1, 1), Exponential(0.5773502691896257)
@show distx, std(distx)
@show disty, std(disty)

a = @show tieshift(distx, disty)
P1 = plot(distx, -2, 4; label="distx")
plot!(disty + a, -2, 4; label="disty + ($(round(a; digits=3)))", ls=:dash)
title!("case of tie shifting")

a = @show median(distx) - median(disty)
P2 = plot(distx, -2, 4; label="distx")
plot!(disty + a, -2, 4; label="disty + ($(round(a; digits=3)))", ls=:dash)
vline!([median(distx)]; label="median(distx)", ls=:dot, lw=1.5)
title!("case of matching medians")

plot(P1, P2; size=(800, 250))

# %%
