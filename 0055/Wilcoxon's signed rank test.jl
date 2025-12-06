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

# %% [markdown]
# # Wilcoxonの符号順位検定の頑健性チェック
#
# * 黒木玄
# * 作成: 2025-12-06
# * 文脈: https://x.com/genkuroki/status/1996556375800025396
# * Colab: https://colab.research.google.com/github/genkuroki/public/blob/main/0055/Wilcoxon's%20signed%20rank%20test.ipynb

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
9
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
using Random

_nthreads() = Threads.nthreads(:interactive) + Threads.nthreads(:default)

using Random
using Statistics

@autoadd begin
using Distributions
using HypothesisTests
using Plots
end

default(fmt=:png, legend=false, size=(400, 400), titlefontsize=10, tickfontsize=7)

Base.show(io::IO, d::InverseGamma) =
    print(io, "InverseGamma(α=", d.invd.α, ", θ=", d.θ, ")")
Base.show(io::IO, ::MIME"text/plain", d::InverseGamma) = Base.show(io, d)
@eval Distributions begin
function logpdf(d::InverseGamma, x::Real)
    x ≤ 0 && return -Inf
    (α, θ) = params(d)
    α * log(θ) - loggamma(α) - (α + 1) * log(x) - θ / x
end
end

distname(dist) = replace(string(dist), r"{[^\}]*}"=>"")

_ecdf(A, x) = count(≤(x), A) / length(A)

function median_x1_plus_x2_monte_carlo(dist::ContinuousUnivariateDistribution; niters=10^8, naves=1)
    s = 0.0
    X1_plus_X2 = zeros(niters)
    for j in 1:naves
        Threads.@threads for i in 1:niters
            X1_plus_X2[i] = rand(dist) + rand(dist)
        end
        s += median!(X1_plus_X2)
    end
    s / naves
end

function prob_x1_plus_x2_gt_zero_monte_carlo(dist; niters=10^8, naves=1)
    s = 0.0
    X1_plus_X2 = zeros(niters)
    for j in 1:naves
        Threads.@threads for i in 1:niters
            X1_plus_X2[i] = rand(dist) + rand(dist)
        end
        s += count(>(0), X1_plus_X2) / length(X1_plus_X2)
    end
    s / naves
end

@show dist = Gamma(0.03, 10000)
@time m1 = median_x1_plus_x2_monte_carlo(dist)
@time m2 = median_x1_plus_x2_monte_carlo(dist)
@time m3 = median_x1_plus_x2_monte_carlo(dist)
@show m1 m2 m3
m = m1
@time p1 = prob_x1_plus_x2_gt_zero_monte_carlo(dist - m/2)
@time p2 = prob_x1_plus_x2_gt_zero_monte_carlo(dist - m/2)
@time p3 = prob_x1_plus_x2_gt_zero_monte_carlo(dist - m/2)
@show p1 p2 p3

# %%
function sim_signed_rank_test(dist::ContinuousUnivariateDistribution, n;
        niters = 10^5,
        distshift = median_x1_plus_x2_monte_carlo(dist) / 2
    )
    dist_null = dist - distshift
    nth = _nthreads()
    Xtmp = [zeros(n) for _ in 1:nth]
    pval = zeros(niters)
    Threads.@threads :static for i in 1:niters
        tid = Threads.threadid()
        X = rand!(dist_null, Xtmp[tid])
        pval[i] = pvalue(SignedRankTest(X))
    end
    pval
end

function plot_sim_signed_rank_test(dist::ContinuousUnivariateDistribution, n;
        niters = 10^5,
        str_dist = distname(dist), 
        distshift = median_x1_plus_x2_monte_carlo(dist) / 2,
        dist_null = dist - distshift,
        μ = mean(dist_null),
        σ = std(dist_null),
        xmin = μ - 5σ,
        xmax = μ + 5σ,
        xs = range(xmin, xmax, 2000),
        ymax = min(10, maximum(pdf(dist_null, x) for x in xs)),
        amax = 0.1,
        r = x -> round(x; digits=4),
        rsd = x -> round(x; sigdigits=4),
    )
    str_dist_null = str_dist * (distshift + 1 ≈ 1 ? "" : (distshift > 0 ? " - " : " + ") * "$(rsd(abs(distshift)))")
    
    println("Xᵢ ~ ", str_dist_null) 
    println("P(X₁ + X₂ > 0) = ", r(prob_x1_plus_x2_gt_zero_monte_carlo(dist_null)))
    println("P(X₁ > 0) = ", r(ccdf(dist_null, 0)))
    println("median(X₁)) = ", r(median(dist_null)))
    println("mean(X₁)) = ", r(mean(dist_null)))
    
    pval = sim_signed_rank_test(dist, n; niters, distshift)
    P = plot(xs, x -> pdf(dist_null, x))
    plot!(ylim=(-0.05ymax, 1.03ymax))
    title!("$str_dist_null,  n=$n")
    Q = plot(α -> _ecdf(pval, α), 0.0, amax)
    plot!(identity; ls=:dot, c=:black, alpha=0.5, lw=0.8)
    xtick = ytick = amax > 0.14 ? (0:0.1:1) : (0:0.01:1)
    plot!(; xtick, ytick)
    plot!(xguide="α", yguide="probability of P-value ≤ α")
    plot!(P, Q; size=(400, 560), layout=@layout[a{0.3h}; b])
end

# %%
plot_sim_signed_rank_test(Normal(), 10)

# %%
plot_sim_signed_rank_test(InverseGamma(2.01), 10; xmin=-2, xmax=5)

# %%
plot_sim_signed_rank_test(LogNormal(), 10; xmin=-2, xmax=6)

# %%
plot_sim_signed_rank_test(Exponential(), 10; xmin=-1, xmax=3.5, amax=1)

# %%
plot_sim_signed_rank_test(Exponential(), 10; xmin=-1, xmax=3.5, amax=1)

# %%
plot_sim_signed_rank_test(Exponential(), 10; xmin=-1, xmax=3.5)

# %%
plot_sim_signed_rank_test(Exponential(), 10; xmin=-1, xmax=3.5)

# %%
plot_sim_signed_rank_test(Gamma(0.1, 10), 40; xmin=-0.12, xmax=0.02, ymax=20, amax=1)

# %%
plot_sim_signed_rank_test(Gamma(0.1, 10), 40; xmin=-0.12, xmax=0.02, ymax=20, amax=1)

# %%
plot_sim_signed_rank_test(Gamma(0.1, 10), 40; xmin=-0.12, xmax=0.02, ymax=20)

# %%
plot_sim_signed_rank_test(Gamma(0.1, 10), 40; xmin=-0.12, xmax=0.02, ymax=20)

# %%
plot_sim_signed_rank_test(Gamma(0.05, 20), 40; xmin=-0.03, xmax=0.11, amax=1)

# %%
plot_sim_signed_rank_test(Gamma(0.05, 20), 40; xmin=-0.03, xmax=0.11)

# %%
plot_sim_signed_rank_test(Gamma(0.05, 20), 40; xmin=-0.03, xmax=0.11)

# %%
