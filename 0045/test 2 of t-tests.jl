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
    #@show T_Student ≈ T_Welch
    #println()
    
    @show df_Student
    @show mean(DF_Welch) std(DF_Welch)
    println()

    pval_Student = @. 2ccdf(TDist(df_Student), abs(T_Student))
    pval_Welch = @. 2ccdf(TDist(DF_Welch), abs(T_Welch))

    @show ECDF(pval_Welch, 0.05)
    #@show ECDF(pval_Student, 0.05)
    println()
    #@show ECDF(pval_Welch, 0.01)
    #@show ECDF(pval_Student, 0.01)
    #println()
    
    P = plot()
    plot!(α -> ECDF(pval_Welch, α), 0, 0.1; label="Welch t-test")
    #plot!(α -> ECDF(pval_Student, α), 0, 0.1; label="Student t-test", ls=:dot)
    plot!(identity, 0, 0.1; label="", ls=:dot, c=:black)
    plot!(; xtick=0:0.01:1, ytick, xrotation=90)
    plot!(; xguide="nominal significance level α",
        yguide="probability of p-value ≤ α")
    plot!(; legend=:outertop, legendfontsize=10)

    Q = plot()
    stephist!(pval_Welch; norm=true, bin, label="Welch t-test")
    #stephist!(pval_Student; norm=true, bin, label="Student t-test", ls=:dot)
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
ff = [0, 0, 100, 0, 100, 0, 0]
ff = [50, 0, 0, 0, 150, 0, 0]
ff = [5, 5, 30, 105, 55, 0, 0]
ff = [5, 10, 40, 70, 75, 0, 0]
ff = [10, 15, 30, 45, 90, 0, 0]
#ff = reverse(ff)
pp = ff / sum(ff)
dist = Categorical(pp)
@show mean(dist) std(dist)
bar(dist; alpha=0.3, label="", size=(300, 200))

# %%
fx = [5, 10, 40, 70, 75, 0, 0]
px = fx / sum(fx)
distx = Categorical(px)

fy = [10, 15, 30, 45, 90, 0, 0]
py = fy / sum(fy)
disty = Categorical(py)

P = bar(distx; alpha=0.3, label="distx", ylim=(-0.05, 0.5))
Q = bar(disty; alpha=0.3, label="disty", ylim=(-0.05, 0.5))
plot(P, Q; size=(600, 200)) |> display

plot_t_tests(; distx, disty, m=50, n=35, ytick=0:0.01:1, ylim=(-0.05, 2), bin=50)

# %%
fx = [5, 10, 40, 70, 75, 0, 0]
px = fx / sum(fx)
distx = Categorical(px)

fy = [10, 15, 30, 45, 90, 0, 0]
py = fy / sum(fy)
disty = Categorical(py)

P = bar(distx; alpha=0.3, label="distx", ylim=(-0.05, 0.5))
Q = bar(disty; alpha=0.3, label="disty", ylim=(-0.05, 0.5))
plot(P, Q; size=(600, 200)) |> display

plot_t_tests(; distx, disty, m=20, n=14, ytick=0:0.01:1, ylim=(-0.05, 2), bin=50)

# %% tags=[]
fx = [5, 10, 40, 70, 75, 0, 0]
px = fx / sum(fx)
distx = Categorical(px)

fy = reverse([10, 15, 30, 45, 90, 0, 0])
py = fy / sum(fy)
disty = Categorical(py)

P = bar(distx; alpha=0.3, label="distx", ylim=(-0.05, 0.5))
Q = bar(disty; alpha=0.3, label="disty", ylim=(-0.05, 0.5))
plot(P, Q; size=(600, 200)) |> display

plot_t_tests(; distx, disty, m=20, n=14, ytick=0:0.01:1, ylim=(-0.05, 2), bin=0:0.01:1.01)

# %% tags=[]
fx = [5, 10, 40, 70, 75, 0, 0]
px = fx / sum(fx)
distx = Categorical(px)

fy = reverse([10, 15, 30, 45, 90, 0, 0])
py = fy / sum(fy)
disty = Categorical(py)

P = bar(distx; alpha=0.3, label="distx", ylim=(-0.05, 0.5))
Q = bar(disty; alpha=0.3, label="disty", ylim=(-0.05, 0.5))
plot(P, Q; size=(600, 200)) |> display

plot_t_tests(; distx, disty, m=50, n=35, ytick=0:0.01:1, ylim=(-0.05, 2), bin=0:0.01:1.01)

# %% tags=[]
fx = [5, 10, 40, 70, 75, 0, 0]
px = fx / sum(fx)
distx = Categorical(px)

fy = reverse([10, 15, 30, 45, 90, 0, 0])
py = fy / sum(fy)
disty = Categorical(py)

P = bar(distx; alpha=0.3, label="distx", ylim=(-0.05, 0.5))
Q = bar(disty; alpha=0.3, label="disty", ylim=(-0.05, 0.5))
plot(P, Q; size=(600, 200)) |> display

plot_t_tests(; distx, disty, m=100, n=70, ytick=0:0.01:1, ylim=(-0.05, 2), bin=0:0.01:1.01)

# %%
distx = Gamma(2, 1/2)
disty = Exponential()

P = plot(distx; label="distx", xlim=(-0.2, 6), ylim=(-0.05, 1.1))
Q = plot(disty; label="disty", xlim=(-0.2, 6), ylim=(-0.05, 1.1))
plot(P, Q; size=(600, 200)) |> display

plot_t_tests(; distx, disty, m=50, n=35, ytick=0:0.01:1, ylim=(-0.05, 2), bin=0:0.01:1.01)

# %%
distx = Gamma(2, 1/2) + 0.585
disty = Exponential()

P = plot(distx, -0.2, 6; label="distx", xlim=(-0.2, 6), ylim=(-0.05, 1.1))
Q = plot(disty, -0.2, 6; label="disty", xlim=(-0.2, 6), ylim=(-0.05, 1.1))
plot(P, Q; size=(600, 200)) |> display

plot_t_tests(; distx, disty, m=50, n=35, ytick=0:0.1:1, ylim=(-0.05, 2), bin=0:0.01:1.01)

# %%
distx = Gamma(2, 1/2)
distx = distx - 1
disty = Exponential()
disty = disty - 1

P = plot(distx, -5, 5; label="distx", xlim=(-5, 5), ylim=(-0.05, 1.1))
Q = plot(x -> pdf(disty, -x), -5, 5; label="disty", xlim=(-5, 5), ylim=(-0.05, 1.1))
plot(P, Q; size=(600, 200)) |> display

plot_t_tests(; distx, disty, m=50, n=35, ytick=0:0.01:1, ylim=(-0.05, 2), bin=0:0.01:1.01, s=-1)

# %%
distx = Gamma(2, 1/2)
distx = distx - 1
disty = Exponential()
disty = disty - 1

P = plot(distx, -5, 5; label="distx", xlim=(-5, 5), ylim=(-0.05, 1.1))
Q = plot(x -> pdf(disty, -x), -5, 5; label="disty", xlim=(-5, 5), ylim=(-0.05, 1.1))
plot(P, Q; size=(600, 200)) |> display

plot_t_tests(; distx, disty, m=100, n=70, ytick=0:0.01:1, ylim=(-0.05, 2), bin=0:0.01:1.01, s=-1)

# %%
distx = Gamma(2, 1/2)
distx = distx - 1
disty = Exponential()
disty = disty - 1

P = plot(distx, -5, 5; label="distx", xlim=(-5, 5), ylim=(-0.05, 1.1))
Q = plot(x -> pdf(disty, -x), -5, 5; label="disty", xlim=(-5, 5), ylim=(-0.05, 1.1))
plot(P, Q; size=(600, 200)) |> display

plot_t_tests(; distx, disty, m=200, n=140, ytick=0:0.01:1, ylim=(-0.05, 2), bin=0:0.01:1.01, s=-1)

# %%
distx = Gamma(2, 1/2)
distx = distx - 1 + 0.5218
disty = Exponential()
disty = disty - 1

P = plot(distx, -5, 5; label="distx", xlim=(-5, 5), ylim=(-0.05, 1.1))
Q = plot(x -> pdf(disty, -x), -5, 5; label="disty", xlim=(-5, 5), ylim=(-0.05, 1.1))
plot(P, Q; size=(600, 200)) |> display

plot_t_tests(; distx, disty, m=50, n=35, ytick=0:0.1:1, ylim=(-0.05, 2), bin=0:0.01:1.01, s=-1)

# %%
