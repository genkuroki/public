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
#     display_name: Julia 1.9.4
#     language: julia
#     name: julia-1.9
# ---

# %%
using Printf
using Random
Random.seed!(4649373)

using Distributions
using StatsPlots
default(fmt=:png, titlefontsize=10, size=(500, 350))

ECDF(A, x) = count(≤(x), A) / length(A)

# %%
distx = Normal(0, 1)
m = 200
disty = Normal(0, 4)
n = 50
L = 4000
Nshuffles = 4000

@show se = √(var(distx)/m + var(disty)/n)
@show se_shuffle = √(var(disty)/m + var(distx)/n)

nth = Threads.nthreads()
XYtmp = [zeros(m+n) for _ in 1:nth]

pval = zeros(L)
pval2 = zeros(L)
pval3 = zeros(L)
diff = zeros(L)
@time Threads.@threads for i in eachindex(pval)
    tid = Threads.threadid()
    XY = XYtmp[tid]
    X = rand!(distx, @view XY[1:m])
    Y = rand!(disty, @view XY[m+1:m+n])
    diff[i] = mean(X) - mean(Y)
    diff_shuffle = zeros(Nshuffles)
    for j in 1:Nshuffles
        shuffle!(XY)
        diff_shuffle[j] = @views mean(XY[1:m]) - mean(XY[m+1:m+n])
    end
    P = 2ECDF(diff_shuffle, diff[i])
    pval[i] = min(1, P, 2-P)
    P2 = 2cdf(Normal(0, se_shuffle), diff[i])
    pval2[i] = min(1, P2, 2-P2)
    P3 = 2cdf(Normal(0, se), diff[i])
    pval3[i] = min(1, P3, 2-P3)
end

@show mean(diff) std(diff)

@show ECDF.((pval, pval2, pval3), 0.10)
@show ECDF.((pval, pval2, pval3), 0.05)
plot()
stephist!(pval; norm=true, label="", bin=0:0.025:1.025)
stephist!(pval2; norm=true, label="", bin=0:0.025:1.025, ls=:dash)
stephist!(pval3; norm=true, label="", bin=0:0.025:1.025, ls=:dashdot)
plot!(xtick=0:0.1:1)

# %%
distx, m = Normal(0, 1), 200
disty, n = Normal(0, 4), 50
X = rand(distx, m)
Y = rand(disty, n)
@show diff = mean(X) - mean(Y)
XY = [X; Y]
Nshuffles = 10^4
diff_shuffle = zeros(Nshuffles)
for j in 1:Nshuffles
    shuffle!(XY)
    diff_shuffle[j] = @views mean(XY[1:m]) - mean(XY[m+1:m+n])
end
@show se2_true = var(distx)/m + var(disty)/n
@show se2_shuffle = var(distx)/n + var(disty)/m
@show mean(diff_shuffle), var(diff_shuffle)
stephist(diff_shuffle; norm=true, label="")
plot!(Normal(0, √se2_shuffle); label="", ls=:dash)
plot!(Normal(0, √se2_true); label="", ls=:dot)
vline!([diff]; label="", ls=:dot)

# %%
function welch_t_test(X, Y; μ = 0.0)
    m, X̄, SX2 = length(X), mean(X), var(X)
    n, Ȳ, SY2 = length(Y), mean(Y), var(Y)
    sehat2 = SX2/m + SY2/n
    tvalue = (X̄ - Ȳ - μ) / √sehat2
    df = sehat2^2 / ((SX2/m)^2/(m-1) + (SY2/n)^2/(n-1))
    pvalue = 2ccdf(TDist(df), abs(tvalue))
    (; pvalue, tvalue, sehat2, df)
end

pvalue_welch_t_test(X, Y; μ = 0.0) =welch_t_test(X, Y; μ).pvalue

using HypothesisTests
X = randn(100)
Y = randn(100)
pvalue_welch_t_test(X, Y), pvalue(UnequalVarianceTTest(X, Y))

# %%
function student_t_test(X, Y; μ = 0.0)
    m, X̄, SX2 = length(X), mean(X), var(X)
    n, Ȳ, SY2 = length(Y), mean(Y), var(Y)
    S2 = ((m-1)*SX2 + (n-1)*SY2) / (m+n-2)
    sehat2 = S2 * (1/m + 1/n)
    tvalue = (X̄ - Ȳ - μ) / √sehat2
    df = m + n - 2
    pvalue = 2ccdf(TDist(df), abs(tvalue))
    (; pvalue, tvalue, sehat2, df)
end

pvalue_student_t_test(X, Y; μ = 0.0) = student_t_test(X, Y; μ).pvalue

using HypothesisTests
X = randn(100)
Y = randn(100)
pvalue_student_t_test(X, Y), pvalue(EqualVarianceTTest(X, Y))

# %%
"""
gammadist(σ, β)

returns the gamma distribution with standard deviation `σ` and skewness `β`.
"""
gammadist(σ, β) = β == 0 ? Normal(0, σ) : Gamma(4/β^2, β*σ/2)

gam = gammadist.(1:5, 2:2:10)
[std.(gam), skewness.(gam), kurtosis.(gam), shape.(gam)]

# %%
"""
inversegammadist(σ, β)

returns the inverse gamma distribution with standard deviation `σ` and skewness `β`.
"""
function inversegammadist(σ, β)
    β == 0 && return Normal(0, σ)
    α = 3 + 8/β^2 + 4/β * √(1 + 4/β^2)
    θ = σ * (α - 1) * √(α - 2)
    InverseGamma(α, θ)
end

igam = inversegammadist.(1:5, 2:2:10)
[std.(igam), skewness.(igam), kurtosis.(igam), shape.(igam)]

# %%
function plot_sim_diff_shuffle(;
        distx = Normal(0, 1), m = 200,
        disty = Normal(0, 4), n = 50,
        difffunc1 = (X, Y) -> mean(X) - mean(Y),
        name1 = "perm(Δmean)",
        difffunc2 = (X, Y) -> (mean(X) - mean(Y)) / √(var(X)/length(X) + var(Y)/length(Y)),
        name2 = "perm(T-stat)",
        L = 5000,
        Nshuffles = 5000,
        α = 0.05,
        xtick = 0:0.01:1,
        ytick = 0:0.01:1,
        kwargs...,
    )

    @show distx disty m n
#     @show se_true = √(var(distx)/m + var(disty)/n)
#     @show se_shuffle = √(var(distx)/n + var(disty)/m)
    println()

    nth = Threads.nthreads()
    XYtmp = [zeros(m+n) for _ in 1:nth]
    diffshuffle1tmp = [zeros(Nshuffles) for _ in 1:nth]
    diffshuffle2tmp = [zeros(Nshuffles) for _ in 1:nth]

    diff1 = zeros(L)
    diff2 = zeros(L)
    pval1 = zeros(L)
    pval2 = zeros(L)
    pval_welch = zeros(L)
    pval_student = zeros(L)
    Threads.@threads for i in 1:L
        tid = Threads.threadid()
        XY = XYtmp[tid]
        diffshuffle1 = diffshuffle1tmp[tid]
        diffshuffle2 = diffshuffle2tmp[tid]
        X = rand!(distx, @view XY[1:m])
        Y = rand!(disty, @view XY[m+1:m+n])
        diff1[i] = difffunc1(X, Y)
        diff2[i] = difffunc2(X, Y)
        pval_welch[i] = pvalue_welch_t_test(X, Y)
        pval_student[i] = pvalue_student_t_test(X, Y)
        for j in 1:Nshuffles
            shuffle!(XY)
            @views X, Y = XY[1:m], XY[m+1:m+n]
            diffshuffle1[j] = difffunc1(X, Y)
            diffshuffle2[j] = difffunc2(X, Y)
        end
        P1 = 2ECDF(diffshuffle1, diff1[i])
        pval1[i] = min(1, P1, 2-P1)
        P2 = 2ECDF(diffshuffle2, diff2[i])
        pval2[i] = min(1, P2, 2-P2)
    end

#     @show mean(diff1) std(diff1)
#     @show mean(diff2) std(diff2)
#     println()
    
    er1 = ECDF(pval1, α)
    er2 = ECDF(pval2, α)
    er_student = ECDF(pval_student, α)
    er_welch = ECDF(pval_welch, α)
    @printf "Probabilities of P-value ≤ %3.1f%%\n" 100α
    @printf "  %-15s %4.1f%%\n" name1*":" 100er1
    @printf "  %-15s %4.1f%%\n" "Student:" 100er_student
    @printf "  %-15s %4.1f%%\n" name2*":" 100er2
    @printf "  %-15s %4.1f%%\n" "Welch:" 100er_welch
    println()
    
    αs = range(0, 0.1, 1001)
    Q = plot()
    plot!(αs, α -> ECDF(pval1, α), label=name1, bin=0:0.025:1.025)
    plot!(αs, α -> ECDF(pval_student, α), label="Student", bin=0:0.025:1.025, ls=:dash)
    plot!(αs, α -> ECDF(pval2, α), label=name2, bin=0:0.025:1.025, ls=:dashdot)
    plot!(αs, α -> ECDF(pval_welch, α), label="Welch", bin=0:0.025:1.025, ls=:dashdotdot)
    plot!(αs, identity; label="", ls=:dot, c=:black, alpha=0.7)
    plot!(; xtick, ytick)
    plot!(; xguide="α", yguide="probablity of P-value ≤ α")
    plot!(; size=(400, 400))
    plot!(; kwargs...)
    
#    plot(P, Q; size=(400, 400))
end

# %% [markdown]
# # 頑健性/脆弱性

# %%
distx, disty = Normal(0, 2), Normal(0, 1)
m, n = 50, 50
L = Nshuffles = 5000
@time plot_sim_diff_shuffle(; distx, m, disty, n, L, Nshuffles)

# %%
distx, disty = Normal(0, 2), Normal(0, 1)
m, n = 25, 100
L = Nshuffles = 5000
@time plot_sim_diff_shuffle(; distx, m, disty, n, L, Nshuffles)

# %%
distx, disty = Normal(0, 2), Normal(0, 1)
m, n = 100, 25
L = Nshuffles = 5000
@time plot_sim_diff_shuffle(; distx, m, disty, n, L, Nshuffles)

# %%
distx, disty = Normal(0, 2), Normal(0, 1)
m, n = 70, 100
L = Nshuffles = 5000
@time plot_sim_diff_shuffle(; distx, m, disty, n, L, Nshuffles)

# %%
distx, disty = Normal(0, 1.5), Normal(0, 1)
m, n = 50, 100
L = Nshuffles = 5000
@time plot_sim_diff_shuffle(; distx, m, disty, n, L, Nshuffles)

# %%
distx, disty = inversegammadist(1, 3), inversegammadist(1, 3)
m, n = 50, 50
L = Nshuffles = 5000
@time plot_sim_diff_shuffle(; distx, m, disty, n, L, Nshuffles)

# %%
distx, disty = inversegammadist(1, 3), inversegammadist(1, 3)
m, n = 25, 100
L = Nshuffles = 5000
@time plot_sim_diff_shuffle(; distx, m, disty, n, L, Nshuffles)

# %%
distx, disty = inversegammadist(1, 3), inversegammadist(1, 3)
m, n = 100, 400
L = Nshuffles = 5000
@time plot_sim_diff_shuffle(; distx, m, disty, n, L, Nshuffles)

# %%
distx, disty = inversegammadist(2, 3), inversegammadist(1, 3)
distx = distx + mean(disty) - mean(distx)
m, n = 50, 50
L = Nshuffles = 5000
@time plot_sim_diff_shuffle(; distx, m, disty, n, L, Nshuffles)

# %%
distx, disty = inversegammadist(2, 3), inversegammadist(1, 3)
distx = distx + mean(disty) - mean(distx)
m, n = 100, 100
L = Nshuffles = 5000
@time plot_sim_diff_shuffle(; distx, m, disty, n, L, Nshuffles)

# %%
distx, disty = inversegammadist(2, 3), inversegammadist(1, 3)
distx = distx + mean(disty) - mean(distx)
m, n = 25, 100
L = Nshuffles = 5000
@time plot_sim_diff_shuffle(; distx, m, disty, n, L, Nshuffles)

# %%
distx, disty = inversegammadist(2, 3), inversegammadist(1, 3)
distx = distx + mean(disty) - mean(distx)
m, n = 100, 400
L = Nshuffles = 5000
@time plot_sim_diff_shuffle(; distx, m, disty, n, L, Nshuffles)

# %%
distx, disty = inversegammadist(2, 3), inversegammadist(1, 3)
distx = distx + mean(disty) - mean(distx)
m, n = 100, 25
L = Nshuffles = 5000
@time plot_sim_diff_shuffle(; distx, m, disty, n, L, Nshuffles)

# %% [markdown]
# # 検出力

# %%
distx, disty = Normal(1.5, 1.8), Normal(0, 1)
plot(distx; label="distx")
plot!(disty; label="disty", ls=:dash)
plot!(size=(300, 180))

# %%
Random.seed!(4649373)
distx, disty = Normal(1.5, 1.8), Normal(0, 1)
m, n = 20, 10
L = Nshuffles = 5000
ytick = 0:0.1:1
legend = :bottomright
@time plot_sim_diff_shuffle(; distx, m, disty, n, L, Nshuffles, ytick, legend)

# %%
distx, disty = Normal(0, 1.8), Normal(0, 1)
m, n = 20, 10
L = Nshuffles = 5000
ytick = 0:0.01:1
legend = :bottomright
@time plot_sim_diff_shuffle(; distx, m, disty, n, L, Nshuffles, ytick, legend)

# %%
distx, disty = gammadist(2, 1.12), gammadist(1, 1)
@show mean(distx), mean(disty)
@show std(distx), std(disty)
plot(distx; label="distx")
plot!(disty; label="disty", ls=:dash)
plot!(size=(300, 180), xlim=(0, 12))

# %%
Random.seed!(4649373)
distx, disty = gammadist(2, 1.12), gammadist(1, 1)
m, n = 20, 10
L = Nshuffles = 5000
ytick = 0:0.1:1
legend = :bottomright
@time plot_sim_diff_shuffle(; distx, m, disty, n, L, Nshuffles, ytick, legend)

# %%
Random.seed!(4649373)
distx, disty = gammadist(2, 1.12), gammadist(1, 1)
distx = distx + mean(disty) - mean(distx)
m, n = 10, 20
L = Nshuffles = 5000
ytick = 0:0.01:1
legend = :bottomright
@time plot_sim_diff_shuffle(; distx, m, disty, n, L, Nshuffles, ytick, legend)

# %% [markdown]
# # ベンチマーク

# %%
distx, disty = Normal(0, 2), Normal(0, 1)
m, n = 25, 100
L = Nshuffles = 500
@time plot_sim_diff_shuffle(; distx, m, disty, n, L, Nshuffles)

# %%
distx, disty = Normal(0, 2), Normal(0, 1)
m, n = 25, 100
L = Nshuffles = 1000
@time plot_sim_diff_shuffle(; distx, m, disty, n, L, Nshuffles)

# %%
distx, disty = Normal(0, 2), Normal(0, 1)
m, n = 25, 100
L = Nshuffles = 2000
@time plot_sim_diff_shuffle(; distx, m, disty, n, L, Nshuffles)

# %%
distx, disty = Normal(0, 2), Normal(0, 1)
m, n = 25, 100
L = Nshuffles = 3000
@time plot_sim_diff_shuffle(; distx, m, disty, n, L, Nshuffles)

# %%
distx, disty = Normal(0, 2), Normal(0, 1)
m, n = 25, 100
L = Nshuffles = 4000
@time plot_sim_diff_shuffle(; distx, m, disty, n, L, Nshuffles)

# %%
