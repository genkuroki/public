# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:hydrogen
#     text_representation:
#       extension: .jl
#       format_name: hydrogen
#       format_version: '1.3'
#       jupytext_version: 1.19.6
#   kernelspec:
#     display_name: Julia
#     language: julia
#     name: julia
# ---

# %% [markdown]
# # 2標本検定達のαエラーについて (正規母集団の場合)
#
# Studentのt検定, Welchのt検定, Mann-WhitneyのU検定, Brunner-Munzel検定の正規母集団に関するαエラーを比較する.
#
# * 黒木 玄
# * 2026-11-02
# * [Repository](https://github.com/genkuroki/public/blob/main/0057/alpha-errors%20of%202-sample%20tests.ipynb)
# * [Google Colab で実行](https://colab.research.google.com/github/genkuroki/public/blob/main/0057/alpha-errors%20of%202-sample%20tests.ipynb)

# %%
# Google Colabと自分のパソコンの両方で使えるようにするための工夫

haskey(ENV, "COLAB_GPU") && (ENV["JULIA_PKG_PRECOMPILE_AUTO"] = "0")
using Pkg

"""すでにPkg.add済みのパッケージのリスト"""
_packages_added = [sort!(readdir(Sys.STDLIB));
    sort!([info.name for (uuid, info) in Pkg.dependencies() if info.is_direct_dep])]

"""_packages_added内にないパッケージをPkg.addする"""
add_pkg_if_not_added_yet(pkg) = if isnothing(Base.find_package(pkg))
    println(stderr, "# $(pkg).jl is not added yet, so let's add it.")
    Pkg.add(pkg)
end

"""expr::Exprからusing内の`.`を含まないモジュール名を抽出"""
function find_using_pkgs(expr::Expr)
    pkgs = String[]
    function traverse(expr::Expr)
        if expr.head == :using
            for arg in expr.args
                if arg.head == :. && length(arg.args) == 1
                    push!(pkgs, string(arg.args[1]))
                elseif arg.head == :(:) && length(arg.args[1].args) == 1
                    push!(pkgs, string(arg.args[1].args[1]))
                end
            end
        else
            for arg in expr.args arg isa Expr && traverse(arg) end
        end
    end
    traverse(expr)
    pkgs
end

"""必要そうなPkg.addを追加するマクロ"""
macro autoadd(expr)
    pkgs = find_using_pkgs(expr)
    :(add_pkg_if_not_added_yet.($(pkgs)); $expr)
end

# %%
mynthreads() = Threads.nthreads() + Threads.nthreads(:interactive)
using Random # for rand!

@autoadd begin
using Distributions
# using Roots
using Plots
end
default(fmt=:png, titlefontsize=12)

ecdf_(A, x) = count(≤(x), A) / length(A)

# %%
function tvalue_welch(m, x̄, sx², n, ȳ, sy²; Δμ=0)
    (x̄ - ȳ - Δμ) / √(sx²/m + sy²/n)
end

function tvalue_welch(x, y; Δμ=0)
    m, x̄, sx² = length(x), mean(x), var(x)
    n, ȳ, sy² = length(y), mean(y), var(y)
    tvalue_welch(m, x̄, sx², n, ȳ, sy²; Δμ)
end

function degree_of_freedom_welch(m, sx², n, sy²)
    (sx²/m + sy²/n)^2 / ((sx²/m)^2/(m-1) + (sy²/n)^2/(n-1))
end

function degree_of_freedom_welch(x, y)
    m, sx² = length(x), var(x)
    n, sy² = length(y), var(y)
    degree_of_freedom_welch(m, sx², n, sy²)
end

function pvalue_welch(m, x̄, sx², n, ȳ, sy²; Δμ=0)
    t = tvalue_welch(m, x̄, sx², n, ȳ, sy²; Δμ)
    ν = degree_of_freedom_welch(m, sx², n, sy²)
    2ccdf(TDist(ν), abs(t))
end

function pvalue_welch(x, y; Δμ=0)
    m, x̄, sx² = length(x), mean(x), var(x)
    n, ȳ, sy² = length(y), mean(y), var(y)
    pvalue_welch(m, x̄, sx², n, ȳ, sy²; Δμ)
end

function confint_welch(m, x̄, sx², n, ȳ, sy²; α=0.05)
    ν = degree_of_freedom_welch(m, sx², n, sy²)
    c = quantile(TDist(ν), 1-α/2)
    SEhat = √(sx²/m + sy²/n)
    [x̄-ȳ-c*SEhat, x̄-ȳ+c*SEhat]
end

function confint_welch(x, y; α=0.05)
    m, x̄, sx² = length(x), mean(x), var(x)
    n, ȳ, sy² = length(y), mean(y), var(y)
    confint_welch(m, x̄, sx², n, ȳ, sy²; α)
end

# %%
s²_student(m, sx², n, sy²) = ((m-1)*sx² + (n-1)*sy²)/(m+n-2)

function tvalue_student(m, x̄, sx², n, ȳ, sy²; Δμ=0)
    s² = s²_student(m, sx², n, sy²)
    (x̄ - ȳ - Δμ) / √(s²*(1/m + 1/n))
end

function tvalue_student(x, y; Δμ=0)
    m, x̄, sx² = length(x), mean(x), var(x)
    n, ȳ, sy² = length(y), mean(y), var(y)
    tvalue_student(m, x̄, sx², n, ȳ, sy²; Δμ)
end

function pvalue_student(m, x̄, sx², n, ȳ, sy²; Δμ=0)
    t = tvalue_student(m, x̄, sx², n, ȳ, sy²; Δμ)
    2ccdf(TDist(m+n-2), abs(t))
end

function pvalue_student(x, y; Δμ=0)
    m, x̄, sx² = length(x), mean(x), var(x)
    n, ȳ, sy² = length(y), mean(y), var(y)
    pvalue_student(m, x̄, sx², n, ȳ, sy²; Δμ)
end

function confint_student(m, x̄, sx², n, ȳ, sy²; α=0.05)
    c = quantile(TDist(m+n-2), 1-α/2)
    s² = s²_student(m, sx², n, sy²)
    SEhat = √(s²*(1/m + 1/n))
    [x̄-ȳ-c*SEhat, x̄-ȳ+c*SEhat]
end

function confint_student(x, y; α=0.05)
    m, x̄, sx² = length(x), mean(x), var(x)
    n, ȳ, sy² = length(y), mean(y), var(y)
    confint_student(m, x̄, sx², n, ȳ, sy²; α)
end

# %%
hodges_lehmann(X, Y) = median(x - y for x in X, y in Y)

winrate(X, Y) = mean((x < y) + (x == y)/2 for x in X, y in Y)

function brunner_munzel_test(X, Y; p=1/2, α=0.05)
    phat = winrate(X, Y)
    m, n = length(X), length(Y)
    sx2 = 1/(m-1) * sum(x -> (mean(y -> (y < x) + (y == x)/2, Y) - (1 - phat))^2, X)
    sy2 = 1/(n-1) * sum(y -> (mean(x -> (x < y) + (x == y)/2, X) - phat)^2, Y)
    sehat = √(sx2/m + sy2/n)
    tvalue = (phat - p)/sehat
    df = (sx2/m + sy2/n)^2 / ((sx2/m)^2/(m-1) + (sy2/n)^2/(n-1))
    pvalue = sehat > 0 ? 2ccdf(TDist(df), abs(tvalue)) : phat ≈ p ? 1.0 : 0.0
    c = sehat > 0 ? cquantile(TDist(df), α/2) : 0.0
    confint_p = (phat - c*sehat, phat + c*sehat)
    (; p, phat, sehat, tvalue, df, pvalue, α, confint_p)
end

pvalue_brunner_munzel(X, Y; p=1/2) = brunner_munzel_test(X, Y; p).pvalue

function _pvalue_brunner_munzel(p, phat, sehat, df)
    tvalue = (phat - p)/sehat
    sehat > 0 ? 2ccdf(TDist(df), abs(tvalue)) : phat ≈ p ? 1.0 : 0.0
end

# function confint_bm_p_roots(X, Y; α=0.05)
#     f(p) = brunner_munzel_test(X, Y; p).pvalue - α
#     find_zeros(f, -1, 2)
# end

function aminamax(X, Y)
    xmin, xmax = extrema(X)
    ymin, ymax = extrema(Y)
    width = max(xmax, ymax) - min(xmin, ymin)
    xmin - ymax - max(0.1, 0.05width), xmax - ymin + max(0.1, 0.05width)
end

# function tieshift(X, Y; p=1/2)
#     f(a) = winrate(X, Y .+ a) - p
#     amin, amax = aminamax(X, Y)
#     find_zero(f, (amin, amax))
# end

# function confint_bm_tieshift(X, Y; α=0.05)
#     f(a) = brunner_munzel_test(X, Y .+ a).pvalue - α
#     amin, amax = aminamax(X, Y)
#     find_zeros(f, amin, amax)
# end

# %%
function mann_whitney_u_test(X, Y; correct=true)
    m, n = length(X), length(Y)
    phat = mean((x < y) + (x == y)/2 for x in X, y in Y)
    sehat = √((m+n+1)/(12m*n))
    zvalue = (phat - 1/2)/sehat
    correction = correct/(2m*n*sehat)
    pvalue = 2ccdf(Normal(), max(0, abs(zvalue) - correction))
    (; phat, sehat, zvalue, pvalue)
end

pvalue_mann_whitney_u_test(X, Y; correct=true) = mann_whitney_u_test(X, Y; correct).pvalue

# using HypothesisTests
# X = randn(100)
# Y = randn(100)
# @show pvalue_mann_whitney_u_test(X, Y) pvalue(ApproximateMannWhitneyUTest(X, Y))
# X = randn(10)
# Y = randn(10)
# @show pvalue_mann_whitney_u_test(X, Y) pvalue(ApproximateMannWhitneyUTest(X, Y));

# %%
function sim_pvalues(;
        distx = Normal(0, 1),
        disty = Normal(0, 1),
        m = 50,
        n = 50,
        L = 10^5,
        np = true
    )
    pval_s = zeros(L)
    pval_w = zeros(L)
    if np
        pval_mw = zeros(L)
        pval_bm = zeros(L)
    else
        pval_mw = zeros(0)
        pval_bm = zeros(0)
    end
    Xtmp = [zeros(m) for _ in 1:mynthreads()]
    Ytmp = [zeros(n) for _ in 1:mynthreads()]
    Threads.@threads :static for i in 1:L
        tid = Threads.threadid()
        X = rand!(distx, Xtmp[tid])
        Y = rand!(disty, Ytmp[tid])
        pval_s[i] = pvalue_student(X, Y)
        pval_w[i] = pvalue_welch(X, Y)
        if np
            pval_mw[i] = pvalue_mann_whitney_u_test(X, Y)
            pval_bm[i] = pvalue_brunner_munzel(X, Y)
        end
    end
    (; pval_s, pval_w, pval_mw, pval_bm)
end

function powers(;
        distx = Normal(0, 1),
        disty = Normal(0, 1),
        m = 50,
        n = 50,
        L = 10^5,
        α = 0.05,
        np = false,
    )
    (; pval_s, pval_w, pval_mw, pval_bm) = sim_pvalues(; distx, disty, m, n, L, np)
    power_s = ecdf_(pval_s, α)
    power_w = ecdf_(pval_w, α)
    if np
        power_mw = ecdf_(pval_mw, α)
        power_bm = ecdf_(pval_bm, α)
    else
        power_mw = NaN
        power_bm = NaN
    end
    (; power_s, power_w, power_mw, power_bm)
end

function plot_powers(;
        Δμ = 0.0,
        σratios = [0.33, 0.5, 0.7, 1.0, 1.4, 2.0, 3.0],
        xticks = (σratios, string.(σratios)),
        xscale = :log10,
        m = 50,
        n = 50,
        L = 10^5,
        α = 0.05,
        np = false,
        μ₁ = 0.0,
        σ₁ = 1.0,
        kwargs...
    )
    powers_s = zeros(length(σratios))
    powers_w = similar(powers_s)
    if np
        powers_mw = similar(powers_s)
        powers_bm = similar(powers_s)
    else
        powers_mw = zeros(0)
        powers_bm = zeros(0)
    end
    for (i, σratio) in enumerate(σratios)
        distx = Normal(0, 1)
        disty = Normal(Δμ, σratio)
        (; power_s, power_w, power_mw, power_bm) = powers(; distx, disty, m, n, L, α, np)
        powers_s[i] = power_s
        powers_w[i] = power_w
        if np
            powers_mw[i] = power_mw
            powers_bm[i] = power_bm
        end
    end
    plot(σratios, powers_s; marker=:o, msc=:auto, label="Student", c=1)
    plot!(σratios, powers_w; marker=:star, msc=:auto, label="Welch", ls=:dash, c=2)
    np && plot!(σratios, powers_mw; marker=:utriangle, msc=:auto, label="MW", ls=:dashdot, c=3)
    np && plot!(σratios, powers_bm; marker=:diamond, msc=:auto, label="BM", ls=:dashdotdot, c=4)
    plot!(; ylim=(-0.03, 1.03))
    plot!(; xticks, yticks=0:0.05:1, xscale)
    yguide = Δμ == 0 ? "α-error rate" : "power"
    plot!(; xguide="s₂/s₁", yguide)
    Δμstr = Δμ == 0 ? "0" : Δμ != 1 ? "$(Δμ)s₁" : "s₁"
    title!("Δμ=$Δμstr,  n₁=$m,  n₂=$n   (sim_iters=$L)")
    plot!(; size=(500, 300))
    plot!(; kwargs...)
end

function plot_α_error_rates(;
        distx = Normal(0, 1),
        disty = Normal(0, 1),
        m = 50,
        n = 50,
        L = 10^6,
        np = true,
        kwargs...
    )
    @show distx disty
    (; pval_s, pval_w, pval_mw, pval_bm) = sim_pvalues(; distx, disty, m, n, L, np)
    plot(identity, 0, 0.1; label="", lw=0.5, c=:black, ls=:dot)
    plot!(α -> ecdf_(pval_s, α), 0, 0.1; label="Student", c=1)
    plot!(α -> ecdf_(pval_w, α), 0, 0.1; label="Welch", c=2, ls=:dash)
    np && plot!(α -> ecdf_(pval_mw, α), 0, 0.1; label="MW", c=3, ls=:dashdot)
    np && plot!(α -> ecdf_(pval_bm, α), 0, 0.1; label="BM", c=4, ls=:dashdotdot)
    plot!(xtick=0:0.01:1, ytick=0:0.01:1)
    plot!(xguide="α", yguide="probability of P-value ≤ α")
    title!("n₁=$m,  n₂=$n   (sim_iters=$L)")
    plot!(; size=(400, 400), titlefontsize=10)
    plot!(; kwargs...)
end

# %%
@time plot_powers(; Δμ=0.4)

# %%
@time plot_powers(; Δμ=0.5)

# %%
@time plot_powers(; Δμ=0.5, np=true)

# %%
@time plot_powers(; Δμ=0.5, np=true)

# %%
@time plot_powers(; Δμ=0, ylim=(0.045, 0.055), L=10^7, ytick=0:0.001:1)

# %%
@time plot_powers(; Δμ=0, ylim=(0.045, 0.08), L=10^6, ytick=0:0.005:1, np=true)

# %%
PP = []
for n in (50, 60, 70, 80, 100, 150)
    P = plot_powers(; Δμ=0, m=50, n, ylim=(0.00, 0.21), L=10^6, ytick=0:0.01:1)
    push!(PP, P)
end
plot(PP...; size=(1000, 1000), layout=(3, 2))

# %%
PP = []
for n in (50, 60, 70, 80, 100, 150)
    P = plot_powers(; Δμ=0, m=50, n, ylim=(0.00, 0.21), L=10^5, ytick=0:0.01:1, np=true)
    push!(PP, P)
end
plot(PP...; size=(1000, 1000), layout=(3, 2))

# %%
PP = []
for n in (50, 60, 70, 80, 100, 150)
    P = plot_powers(; Δμ=0.5, m=50, n)
    push!(PP, P)
end
plot(PP...; size=(1000, 1000), layout=(3, 2))

# %%
@time plot_powers(; Δμ=0, m=10, n=10, ylim=(0.045, 0.06), L=10^7, ytick=0:0.005:1)

# %%
@time plot_powers(; Δμ=0, m=10, n=10, ylim=(0.04, 0.062), L=10^6, ytick=0:0.005:1, np=true)

# %%
PP = []
for n in (10, 12, 14, 16, 20, 30)
    P = plot_powers(; Δμ=0, m=10, n, ylim=(0.00, 0.22), L=10^6, ytick=0:0.01:1)
    push!(PP, P)
end
plot(PP...; size=(1000, 1000), layout=(3, 2))

# %%
PP = []
for n in (10, 12, 14, 16, 20, 30)
    P = plot_powers(; Δμ=0, m=10, n, ylim=(0.00, 0.22), L=10^5, ytick=0:0.01:1, np=true)
    push!(PP, P)
end
plot(PP...; size=(1000, 1000), layout=(3, 2))

# %%
PP = []
for n in (10, 12, 14, 16, 20, 30)
    P = plot_powers(; Δμ=1, m=10, n, ylim=(0.00, 1.00), L=10^5, ytick=0:0.05:1, np=true)
    push!(PP, P)
end
plot(PP...; size=(1000, 1000), layout=(3, 2))

# %%
plot_α_error_rates(; distx=Normal(0, 1), disty=Normal(0, 1), m=50, n=50)

# %%
plot_α_error_rates(; distx=Normal(0, 1), disty=Normal(0, 4), m=50, n=50)

# %%
plot_α_error_rates(; distx=Normal(0, 1), disty=Normal(0, 4), m=70, n=50)

# %%
plot_α_error_rates(; distx=Normal(0, 1), disty=Normal(0, 4), m=100, n=50)

# %%
distx = Exponential()
disty = Exponential()
@show mean(disty) - mean(distx)
@show std(disty) / std(disty)
plot_α_error_rates(; distx, disty, m=50, n=50, L=10^6, np=false)

# %%
distx = Exponential()
disty = Exponential()
@show mean(disty) - mean(distx)
@show std(disty) / std(disty)
plot_α_error_rates(; distx, disty, m=50, n=100, L=10^6, np=false)

# %%
distx = Exponential()
disty = Exponential()
@show mean(disty) - mean(distx)
@show std(disty) / std(disty)
plot_α_error_rates(; distx, disty, m=50, n=250, L=10^6, np=false)

# %%
distx = Normal(1)
disty = Exponential()
@show mean(disty) - mean(distx)
@show std(disty) / std(disty)
plot_α_error_rates(; distx, disty, m=50, n=50, L=10^6, np=false)

# %%
distx = Normal(1)
disty = Exponential()
@show mean(disty) - mean(distx)
@show std(disty) / std(disty)
plot_α_error_rates(; distx, disty, m=50, n=100, L=10^6, np=false)

# %%
distx = Normal(1)
disty = Exponential()
@show mean(disty) - mean(distx)
@show std(disty) / std(disty)
plot_α_error_rates(; distx, disty, m=100, n=50, L=10^6, np=false)

# %%
distx = Normal(1)
disty = Exponential()
@show mean(disty) - mean(distx)
@show std(disty) / std(disty)
plot_α_error_rates(; distx, disty, m=50, n=250, L=10^6, np=false)

# %%
distx = Normal(1)
disty = Exponential()
@show mean(disty) - mean(distx)
@show std(disty) / std(disty)
plot_α_error_rates(; distx, disty, m=10, n=50, L=10^6, np=false)

# %%
distx = Normal(1)
disty = Exponential()
@show mean(disty) - mean(distx)
@show std(disty) / std(disty)
plot_α_error_rates(; distx, disty, m=250, n=50, L=10^6, np=false)

# %%
