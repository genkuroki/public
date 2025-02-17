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
#     display_name: Julia 1.9.3
#     language: julia
#     name: julia-1.9
# ---

# %%
using Distributions
using HypothesisTests
using Random
using StatsPlots
default(fmt=:png, tickfontsize=6, titlefontsize=12)

ECDF(A, x) = count(≤(x), A)/length(A)

function degree_of_freedom_Welch(m, sx², n, sy²)
    (sx²/m + sy²/n)^2 / ((sx²/m)^2/(m-1) + (sy²/n)^2/(n-1))
end

function plot_t_tests(;
        distx = Normal(),
        disty = Normal(),
        m = 10, 
        n = 10,
        s = 1,
        L = 10^6,
        ytick = 0:0.01:1,
        ylim = :auto,
        legend = :bottom,
        bin = :auto
    )
    println("(distribution of X) = distx = ", distx) 
    println("(distribution of Y) = disty = ", disty)
    println("(size of sample from X) = m = ", m)
    println("(size of sample from Y) = n = ", n)
    if s != 1 println("multiple of disty = ", s) end
    println()
    
    @show mean(distx)
    if s == 1 @show(mean(disty)) else @show(s*mean(disty)) end
    @show std(distx)
    if s == 1 @show(std(disty)) else @show(abs(s)*std(disty)) end
    println()
    
    T_Student = Vector{Float64}(undef, L)
    T_Welch = Vector{Float64}(undef, L)
    
    df_Student = m + n - 2
    DF_Welch = Vector{Float64}(undef, L)
    
    Xtmp = [Vector{eltype(distx)}(undef, m) for _ in 1:Threads.nthreads()]
    Ytmp = [Vector{eltype(disty)}(undef, n) for _ in 1:Threads.nthreads()]
    Threads.@threads for i in 1:L
        tid = Threads.threadid()
        X = rand!(distx, Xtmp[tid])
        Y = rand!(disty, Ytmp[tid])
    
        X̄ = mean(X)
        Ȳ = s*mean(Y)
        SX = std(X)
        SY = std(Y)
        S = √(((m-1)*SX^2 + (n-1)*SY^2)/(m+n-2))

        T_Student[i] = (X̄ - Ȳ)/(S*√(1/m + 1/n))
        T_Welch[i] = (X̄ - Ȳ)/√(SX^2/m + SY^2/n)
        
        DF_Welch[i] = degree_of_freedom_Welch(m, SX^2, n, SY^2)
    end
    @show T_Student ≈ T_Welch
    println()
    
    @show df_Student
    @show mean(DF_Welch) std(DF_Welch)
    println()

    pval_Student = @. 2ccdf(TDist(df_Student), abs(T_Student))
    pval_Welch = @. 2ccdf(TDist(DF_Welch), abs(T_Welch))

    @show ECDF(pval_Student, 0.05) ECDF(pval_Welch, 0.05)
    println()
    @show ECDF(pval_Student, 0.01) ECDF(pval_Welch, 0.01)
    println()
    
    P = plot(α -> ECDF(pval_Student, α), 0, 0.1; label="Student t-test")
    plot!(α -> ECDF(pval_Welch, α), 0, 0.1; label="Welch t-test", ls=:dash)
    plot!(identity, 0, 0.1; label="", ls=:dot, c=:black)
    plot!(; xtick=0:0.01:1, ytick, xrotation=90)
    plot!(; xguide="nominal significance level α",
        yguide="probability of p-value ≤ α")
    plot!(; legend=:outertop, legendfontsize=10)

    Q = stephist(pval_Student; norm=true, bin, label="Student t-test")
    stephist!(pval_Welch; norm=true, bin, ls=:dash, label="Welch t-test")
    plot!(; xguide="P-value", yguide="probability density")
    plot!(; xtick=0:0.05:1, ytick=0:0.1:5, xrotation=90)
    plot!(; legend=:outertop, legendfontsize=10, ylim)
    
    plot(P, Q; size=(1000, 460), layout=@layout [a{0.4w} b])
    plot!(; bottommargin=4Plots.mm, leftmargin=4Plots.mm)
end

function plot_df_and_tstatratio(;
        m = 9,
        n = 11,
        tick = Any[0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100]
    )
    log10tick = (tick, string.(tick))
    as = 10.0 .^ range(log10(minimum(tick)), log10(maximum(tick)), 1000)

    P = plot(as, a -> degree_of_freedom_Welch(m, 1, n, a^2); label="Welch t-test")
    hline!([(m-1)+(n-1)]; label="Student t-test", ls=:dash)
    hline!([min(m-1, n-1)]; label="min(m-1, n-1)", ls=:dashdot)
    hline!([max(m-1, n-1)]; label="max(m-1, n-1)", ls=:dashdotdot)
    vline!([√((n*(n-1))/(m*(m-1)))]; label="√((n*(n-1))/(m*(m-1)))", ls=:dot)
    vline!([1]; label="", ls=:dot, c=:gray)
    plot!(xscale=:log10, xtick=log10tick)
    plot!(xguide="std(sample of Y) / std(sample of X)")
    title!("degree of freedom for m = $m, n = $n")

    Q = plot(as, a -> √(1/m + a^2/n) / √((((m-1)+(n-1)*a^2)/(m-1+n-1))*(1/m+1/n));
        label="t_Student / t_Welch")
    hline!([√(m/n), √(n/m)]; label="√(m/n), √(n/m)", ls=:dash)
    vline!([1]; label="", ls=:dot, c=:gray)
    plot!(xscale=:log10, xtick=log10tick)
    plot!(xguide="std(sample of Y) / std(sample of X)")
    plot!(ylim=(0.9, 1.1) .* minmax(√(m/n), √(n/m)))
    title!("t_Student / t_Welch for m = $m,  n = $n")

    plot(P, Q; size=(1000, 400))
    plot!(bottommargin=4Plots.mm)
end

# %%
plot_t_tests(; distx = Normal(), disty = Normal(), m = 10, n = 10)

# %%
plot_t_tests(; distx = Normal(), disty = Normal(), m = 9, n = 11)

# %%
plot_t_tests(; distx = Normal(), disty = Normal(), m = 8, n = 12)

# %%
plot_t_tests(; distx = Normal(1.3, 1), disty = Normal(0, 1), m = 10, n = 10,
    ytick=0:0.1:1, ylim=(-0.06, 3))

# %%
plot_t_tests(; distx = Normal(0, 2), disty = Normal(0, 1), m = 10, n = 10)

# %%
plot_df_and_tstatratio(; m = 10, n = 10)

# %%
plot_t_tests(; distx = Normal(0, 2), disty = Normal(0, 1), m = 9, n = 11)

# %%
plot_df_and_tstatratio(; m = 9, n = 11)

# %%
plot_t_tests(; distx = Normal(0, 2), disty = Normal(0, 1), m = 8, n = 12)

# %%
plot_df_and_tstatratio(; m = 8, n = 12)

# %%
plot_t_tests(; distx = Normal(0, 2), disty = Normal(0, 1), m = 80, n = 120)

# %%
plot_t_tests(; distx = Normal(0, 2), disty = Normal(0, 1), m = 800, n = 1200)

# %%
plot_t_tests(; distx = Normal(0, 3), disty = Normal(0, 1), m = 9, n = 11)

# %%
plot_t_tests(; distx = Normal(0, 3), disty = Normal(0, 1), m = 90, n = 110)

# %%
plot_t_tests(; distx = Normal(0, 3), disty = Normal(0, 1), m = 11, n = 9)

# %%
plot_t_tests(; distx = Normal(2.94, 3), disty = Normal(0, 1), m = 11, n = 9,
    ytick=0:0.1:1, ylim=(-0.06, 3))

# %%
plot(Normal(2.83, 3), -6, 12; label="Normal(0, 3)")
plot!(Normal(0, 1), -6, 12; label="Normal(0, 1)", ls=:dash)

# %%
plot_t_tests(; distx = Normal(0, 3), disty = Normal(0, 1), m = 12, n = 8)

# %%
plot_t_tests(; distx = Normal(2.83, 3), disty = Normal(0, 1), m = 12, n = 8,
    ytick=0:0.1:1, ylim=(-0.06, 3))

# %%
plot_t_tests(; distx = Normal(2.83, 3), disty = Normal(0, 1), m = 15, n = 10,
    ytick=0:0.1:1, ylim=(-0.06, 3))

# %%
distx = Gamma(8, 4)
disty = Gamma(32, 1)

plot(distx; label="distx")
plot!(disty; label="disty", ls=:dash)
plot!(legendfontsize=12, size=(400, 250)) |> display

plot_t_tests(; distx, disty, m = 20, n = 12)

# %%
distx = Gamma(8, 4) + 8.6
disty = Gamma(32, 1)

plot(distx; label="distx")
plot!(disty; label="disty", ls=:dash)
plot!(legendfontsize=12, size=(400, 250)) |> display

plot_t_tests(; distx, disty, m = 20, n = 12, ytick=0:0.1:1, ylim=(-0.06, 3))

# %%
plot_t_tests(; distx = Normal(0, 1.1), disty = Normal(0, 1), m = 10, n = 50)

# %%
plot_df_and_tstatratio(; m = 10, n = 50)

# %%
plot_t_tests(; distx = Normal(0, 1.1), disty = Normal(0, 1), m = 20, n = 100)

# %%
plot_df_and_tstatratio(; m = 20, n = 100,
    tick = Any[0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500])

# %%
plot_t_tests(; distx = Normal(0, 1.2), disty = Normal(0, 1), m = 10, n = 30)

# %%
plot_df_and_tstatratio(; m = 10, n = 30,
    tick = Any[0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500])

# %%
plot_t_tests(; distx = Normal(0, 1.2), disty = Normal(0, 1), m = 20, n = 60)

# %%
plot_df_and_tstatratio(; m = 20, n = 60,
    tick = Any[0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500])

# %%
plot_t_tests(; distx = Normal(0, 1.2), disty = Normal(0, 1), m = 20, n = 60)

# %%
plot_t_tests(; distx = Gamma(1, 4), disty = Gamma(4, 1), m = 100, n = 100)

# %%
plot_t_tests(; distx = Gamma(1, 4), disty = Gamma(4, 1), m = 90, n = 110)

# %%
plot_t_tests(; distx = Gamma(1, 4), disty = Gamma(4, 1), m = 80, n = 120)

# %%
distx = LogNormal(0.1, 0.7)
disty = LogNormal(0.2, 0.7)
plot_t_tests(; distx, disty, m = 10^3, n = 10^3, ytick=0:0.1:1, ylim=(-0.06, 3))

# %%
plot_df_and_tstatratio(; m=10, n=10)

# %%
plot_df_and_tstatratio(; m=10, n=11)

# %%
plot_df_and_tstatratio(; m=10, n=12)

# %%
plot_df_and_tstatratio(; m=10, n=14)

# %%
plot_df_and_tstatratio(; m=10, n=20)

# %%
plot_df_and_tstatratio(; m=100, n=100)

# %%
plot_df_and_tstatratio(; m=100, n=110)

# %%
plot_df_and_tstatratio(; m=100, n=120)

# %%
plot_df_and_tstatratio(; m=100, n=140)

# %%
plot_df_and_tstatratio(; m=100, n=200)

# %%
z = 1.96
plot(n -> 2ccdf(TDist(n-1), z), 5, 100; label="2ccdf(TDist(n-1), $z)")
hline!([2ccdf(Normal(), z)]; label="2ccdf(Normal(), $z)", ls=:dash)
plot!(xguide="n")
plot!(xtick=0:5:100, ytick=0:0.01:1)

# %%
z = 2.576
plot(n -> 2ccdf(TDist(n-1), z), 5, 100; label="2ccdf(TDist(n-1), $z)")
hline!([2ccdf(Normal(), z)]; label="2ccdf(Normal(), $z)", ls=:dash)
plot!(xguide="n")
plot!(xtick=0:5:100, ytick=0:0.01:1)

# %% tags=[]
@show dist = Normal()
@show c = cquantile(dist, 0.05/2)
for p in 0.07:-0.01:0.03
    @eval @show cquantile(dist, $p/2) / c
end
plot(x -> 2ccdf(dist, x), c/1.1, 1.1c; label="P-value")
plot!(xtick=0:0.1:4, ytick=0:0.01:1)
vline!([cquantile(dist, p/2) for p in 0.04:0.01:0.06]; label="", ls=:dot)
hline!(0.04:0.01:0.06; label="", ls=:dot)
plot!(xguide="test statistic", yguide="P-value")
title!("$dist")

# %% tags=[]
@show dist = TDist(10)
@show c = cquantile(dist, 0.05/2)
for p in 0.07:-0.01:0.03
    @eval @show cquantile(dist, $p/2) / c
end
plot(x -> 2ccdf(dist, x), c/1.1, 1.1c; label="P-value")
plot!(xtick=0:0.1:4, ytick=0:0.01:1)
vline!([cquantile(dist, p/2) for p in 0.04:0.01:0.06]; label="", ls=:dot)
hline!(0.04:0.01:0.06; label="", ls=:dot)
plot!(xguide="test statistic", yguide="P-value")
title!("$dist")

# %%
@show dist = TDist(4)
@show c = cquantile(dist, 0.05/2)
for p in 0.07:-0.01:0.03
    @eval @show cquantile(dist, $p/2) / c
end
plot(x -> 2ccdf(dist, x), c/1.1, 1.1c; label="P-value")
plot!(xtick=0:0.1:4, ytick=0:0.01:1)
vline!([cquantile(dist, p/2) for p in 0.04:0.01:0.06]; label="", ls=:dot)
hline!(0.04:0.01:0.06; label="", ls=:dot)
plot!(xguide="test statistic", yguide="P-value")
title!("$dist")

# %%
plot_df_and_tstatratio(; m=9, n=10)

# %%
plot_t_tests(; distx=Normal(0, 3), disty=Normal(0, 1), m = 9, n = 10)

# %%
plot_t_tests(; distx=Normal(0, 3), disty=Normal(0, 1), m = 27, n = 30)

# %%
plot_t_tests(; distx=Normal(0, 3), disty=Normal(0, 1), m = 90, n = 100)

# %%
plot_t_tests(; distx=Normal(0, 3), disty=Normal(0, 1), m = 10, n = 9)

# %%
plot_t_tests(; distx=Normal(0, 3), disty=Normal(0, 1), m = 30, n = 27)

# %%
plot_t_tests(; distx=Normal(0, 3), disty=Normal(0, 1), m = 90, n = 100)

# %%
plot_df_and_tstatratio(; m=8, n=10)

# %%
plot_t_tests(; distx = Normal(0, 3), disty = Normal(0, 1), m = 8, n = 10)

# %%
plot_t_tests(; distx = Normal(0, 3), disty = Normal(0, 1), m = 24, n = 30)

# %%
plot_t_tests(; distx = Normal(0, 3), disty = Normal(0, 1), m = 80, n = 100)

# %% [markdown]
# https://tanigaku.jp/wp/?p=2556

# %%
X = [153, 153, 152, 156, 158, 160, 162, 151, 151, 150]
Y = [153, 146, 138, 152, 140, 128, 116, 146, 156, 142]

@show X
@show Y
@show m = length(X)
@show n = length(Y)
@show X̄ = mean(X)
@show Ȳ = mean(Y)
@show SX = std(X)
@show SY = std(Y)
@show SY^2/SX^2
@show S = √(((m-1)*SX^2 + (n-1)*SY^2)/(m+n-2))
@show SE_Student = S * √(1/m + 1/n)
@show SE_Welch = √(SX^2/m + SY^2/n)
@show T_Student = (X̄ - Ȳ)/SE_Student
@show T_Welch = (X̄ - Ȳ)/SE_Welch
@show df_Student = m + n - 2
@show DF_Welch = degree_of_freedom_Welch(m, SX^2, n, SY^2)
@show pvalue_Student = 2ccdf(TDist(df_Student), abs(T_Student))
@show pvalue_Welch = 2ccdf(TDist(DF_Welch), abs(T_Welch))
;

# %% tags=[]
@show SY^2/SX^2
@show 2ccdf(FDist(9, 9), SY^2/SX^2)
@show VarianceFTest(Y, X)

SY2overSX2 = [var(randn(10))/var(randn(10)) for _ in 1:10^6]
stephist(SY2overSX2; norm=true, label="var(size-10 normal sample) / var(size-10 normal sample)")
plot!(x -> pdf(FDist(9, 9), x), ls=:dash, label="FDist(9, 9)")
plot!(xlim=(-0.2, 9), legend=:outertop)
vline!([SY^2/SX^2]; ls=:dash, label="SY²/SX²", c=:red)

# %%
plot_t_tests(; distx=Normal(X̄, SX), disty=Normal(X̄, SY), m = 10, n = 10)

# %%
