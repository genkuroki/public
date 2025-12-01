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
# https://okumuralab.org/~okumura/python/t_u_test.html はヒストグラムを使っている点が非常によろしくない。
#
# 関連ノートブック: https://github.com/genkuroki/public/blob/main/0051/Student%20vs.%20Welch.ipynb

# %%
using Random
using Distributions
using Plots
using HypothesisTests
using QuadGK
using Roots

distname(dist) = replace(string(dist), r"{[^}]*}"=>"")
distname(dist::InverseGamma) = "InverseGamma(α=$(shape(dist)), θ=$(scale(dist)))"
safe_pdf(dist, x) = minimum(dist) < x < maximum(dist) ? pdf(dist, x) : zero(eltype(dist))

default(fmt=:png, legend=false, titlefontsize=10, plot_titlefontsize=10, size=(320, 200))

# approximate Mann-Whitney U-testではデフォルトでは連続性補正がかかっている。
# https://github.com/JuliaStats/HypothesisTests.jl/blob/master/src/mann_whitney.jl#L236-L249
my_pvalue(x::ApproximateMannWhitneyUTest; continuity_correction=true) =
    continuity_correction ? pvalue(x) : 2ccdf(Normal(), abs(x.mu)/x.sigma)

function win_rate_of_y_against_x(distx::UnivariateDistribution, disty::ContinuousUnivariateDistribution)
    f(y) = cdf(distx, y) * pdf(disty, y)
    quadgk(f, extrema(disty)...)[1]
end

function equalizing_handicap(distx::UnivariateDistribution, disty::ContinuousUnivariateDistribution; p=1/2)
    f(s) = win_rate_of_y_against_x(distx + s, disty) - p
    find_zero(f, 0.0)
end

function p_values(distx, m, disty, n; niters=10^6, continuity_correction=true,
        handicap_mean = mean(disty) - mean(distx),
        handicap_win_rate = equalizing_handicap(distx, disty),
    )
    exact_mwu = m + n ≤ 66 ? true : false
    name_mwu = exact_mwu ? "exact Mann-Whitney U-test" : "approximate Mann-Whitney U-test"
    pvalue_student_t = zeros(niters)
    pvalue_welch_t = zeros(niters)
    pvalue_mann_whitney_u = zeros(niters)
    pvalue_mann_whitney_u_approximate = zeros(niters)
    nth = Threads.nthreads(:interactive) + Threads.nthreads(:default)
    tmpx = [zeros(m) for _ in 1:nth]
    tmpy = [zeros(n) for _ in 1:nth]
    Threads.@threads :static for i in 1:niters
        tid = Threads.threadid()
        datax = rand!(distx, tmpx[tid])
        datay = rand!(disty, tmpy[tid])
        datax .+= handicap_mean
        pvalue_student_t[i] = pvalue(EqualVarianceTTest(datax, datay))
        pvalue_welch_t[i] = pvalue(UnequalVarianceTTest(datax, datay))
        datax .-= handicap_mean
        datax .+= handicap_win_rate
        if exact_mwu
            pvalue_mann_whitney_u[i] = pvalue(ExactMannWhitneyUTest(datax, datay))
        else
            pvalue_mann_whitney_u[i] = pvalue(ApproximateMannWhitneyUTest(datax, datay))
        end
        pvalue_mann_whitney_u_approximate[i] = my_pvalue(ApproximateMannWhitneyUTest(datax, datay); continuity_correction)
    end
    (; name_mwu, pvalue_student_t, pvalue_welch_t, pvalue_mann_whitney_u, pvalue_mann_whitney_u_approximate)
end
    
function hist_p_values(distx, m, disty, n; niters=10^6, bin=0:0.05:1.05, continuity_correction=true,
        handicap_mean = mean(disty) - mean(distx),
        handicap_win_rate = equalizing_handicap(distx, disty),
    )
    @time (;name_mwu,  pvalue_student_t, pvalue_welch_t, pvalue_mann_whitney_u, pvalue_mann_whitney_u_approximate) =
        p_values(distx, m, disty, n; niters, continuity_correction, handicap_mean, handicap_win_rate)
    P1 = histogram(pvalue_student_t; norm=true, c=1, alpha=0.5, title="Student t-test", bin)
    P2 = histogram(pvalue_welch_t; norm=true, c=2, alpha=0.5, title="Welch t-test", bin)
    P3 = histogram(pvalue_mann_whitney_u; norm=true, c=3, alpha=0.5, title=name_mwu, bin)
    P4 = histogram(pvalue_mann_whitney_u_approximate; norm=true, c=4, alpha=0.5, title="approximate Mann-Whitney U-test", bin)
    plot(P1, P2, P3, P4; size=(640, 400), layout=(2, 2))
    plot!(plot_title="$(distname(distx)), m=$m vs. $(distname(disty)), n=$n")
end

_ecdf(A, x) = count(≤(x), A) / length(A)

function ecdf_p_values(distx, m, disty, n; niters=10^6, maxalpha=1, continuity_correction=true,
        handicap_mean = mean(disty) - mean(distx),
        handicap_win_rate = equalizing_handicap(distx, disty),
    )
    @time (;name_mwu,  pvalue_student_t, pvalue_welch_t, pvalue_mann_whitney_u, pvalue_mann_whitney_u_approximate) =
        p_values(distx, m, disty, n; niters, continuity_correction, handicap_mean, handicap_win_rate)
    P1 = plot(x -> _ecdf(pvalue_student_t, x), 0, maxalpha; c=1, title="Student t-test")
    plot!(identity; ls=:dot, c=:gray, alpha=0.8)
    P2 = plot(x -> _ecdf(pvalue_welch_t, x), 0, maxalpha; c=2, title="Welch t-test")
    plot!(identity; ls=:dot, c=:gray, alpha=0.8)
    P3 = plot(x -> _ecdf(pvalue_mann_whitney_u, x), 0, maxalpha; c=3, title=name_mwu)
    plot!(identity; ls=:dot, c=:gray, alpha=0.8)
    P4 = plot(x -> _ecdf(pvalue_mann_whitney_u_approximate, x), 0, maxalpha; c=4, title="approximate Mann-Whitney U-test")
    plot!(identity; ls=:dot, c=:gray, alpha=0.8)
    plot(P1, P2, P3, P4; size=(640, 640), layout=(2, 2))
    plot!(plot_title="$(distname(distx)), m=$m vs. $(distname(disty)), n=$n")
end

# %% [markdown]
# ## P値の経験分布のヒストグラム

# %%
hist_p_values(Normal(), 10, Normal(1), 10)

# %%
hist_p_values(Normal(), 10, Normal(), 10; continuity_correction=false)

# %%
hist_p_values(Normal(), 10, Normal(), 10; bin=0:0.02:1.02)

# %%
hist_p_values(Normal(), 10, Normal(), 10; bin=0:0.02:1.02, continuity_correction=false)

# %%
hist_p_values(Uniform(), 10, Uniform(), 10)

# %%
hist_p_values(Uniform(), 10, Uniform(), 10; continuity_correction=false)

# %%
hist_p_values(Uniform(), 10, Uniform(), 10; bin=0:0.02:1.02)

# %%
hist_p_values(Uniform(), 10, Uniform(), 10; bin=0:0.02:1.02, continuity_correction=false)

# %%
dist = InverseGamma(3.1, 1)
@show skewness(dist)
plot(x -> pdf(dist, x), 0, 5std(dist); title=distname(dist))

# %%
hist_p_values(InverseGamma(3.1, 1), 10, InverseGamma(3.1, 1), 10)

# %%
hist_p_values(InverseGamma(3.1, 1), 10, InverseGamma(3.1, 1), 10; continuity_correction=false)

# %%
hist_p_values(InverseGamma(3.1, 1), 10, InverseGamma(3.1, 1), 10; bin=0:0.02:1)

# %%
hist_p_values(InverseGamma(3.1, 1), 10, InverseGamma(3.1, 1), 10; bin=0:0.02:1, continuity_correction=false)

# %% [markdown]
# ## P値の経験累積分布関数のグラフ

# %% [markdown]
# ### 等分布の場合

# %%
ecdf_p_values(Normal(), 10, Normal(), 10)

# %%
ecdf_p_values(Normal(), 10, Normal(), 10; continuity_correction=false)

# %%
ecdf_p_values(Uniform(), 10, Uniform(), 10)

# %%
ecdf_p_values(Uniform(), 10, Uniform(), 10; continuity_correction=false)

# %%
dist = InverseGamma(3.1, 1)
@show skewness(dist)
plot(x -> pdf(dist, x), 0, 5std(dist); title=distname(dist))

# %%
ecdf_p_values(InverseGamma(3.1, 1), 10, InverseGamma(3.1, 1), 10)

# %%
ecdf_p_values(InverseGamma(3.1, 1), 10, InverseGamma(3.1, 1), 10; continuity_correction=false)

# %%
ecdf_p_values(InverseGamma(3.1, 1), 10, InverseGamma(3.1, 1), 10; maxalpha=0.1)

# %%
ecdf_p_values(InverseGamma(3.1, 1), 10, InverseGamma(3.1, 1), 10; maxalpha=0.1, continuity_correction=false)

# %%
ecdf_p_values(InverseGamma(3.1, 1), 5, InverseGamma(3.1, 1), 5; maxalpha=0.1)

# %%
ecdf_p_values(InverseGamma(3.1, 1), 5, InverseGamma(3.1, 1), 5; maxalpha=0.1, continuity_correction=false)

# %%
ecdf_p_values(InverseGamma(3.1, 1), 5, InverseGamma(3.1, 1), 10; maxalpha=0.1)

# %%
ecdf_p_values(InverseGamma(3.1, 1), 5, InverseGamma(3.1, 1), 10; maxalpha=0.1, continuity_correction=false)

# %%
ecdf_p_values(InverseGamma(3.1, 1), 5, InverseGamma(3.1, 1), 20; maxalpha=0.1)

# %%
ecdf_p_values(InverseGamma(3.1, 1), 5, InverseGamma(3.1, 1), 20; maxalpha=0.1, continuity_correction=false)

# %% [markdown]
# ### 非等分布の場合

# %% [markdown]
# #### 分散が異なる2つの正規分布の場合

# %%
ecdf_p_values(Normal(0, 1), 10, Normal(0, 3), 10; maxalpha=0.1, continuity_correction=false)

# %%
ecdf_p_values(Normal(0, 1), 30, Normal(0, 3), 10; maxalpha=0.1, continuity_correction=false)

# %%
ecdf_p_values(Normal(0, 1), 10, Normal(0, 3), 30; maxalpha=0.1, continuity_correction=false)

# %%
ecdf_p_values(Normal(0, 1), 100, Normal(0, 3), 100; maxalpha=0.1, continuity_correction=false)

# %%
ecdf_p_values(Normal(0, 1), 300, Normal(0, 3), 100; maxalpha=0.1, continuity_correction=false)

# %%
ecdf_p_values(Normal(0, 1), 100, Normal(0, 3), 300; maxalpha=0.1, continuity_correction=false)

# %% [markdown]
# #### 分散が異なる2つの逆ガンマ分布の場合

# %%
distx = InverseGamma(3.1, 1)
disty = InverseGamma(3.1, 3)
igx, igy = distx, disty

@show distname(distx) distname(disty)
@show std(disty)/std(distx)
@show skewness(distx) skewness(disty)
@show handicap_mean = mean(disty) - mean(distx)
@show handicap_win_rate = equalizing_handicap(distx, disty)
xmin = min(0, handicap_mean, handicap_win_rate)
xmax = max(5std(distx) + handicap_win_rate, 5std(disty))
xs = range(xmin, xmax, 1000)
P1 = plot(xs, x -> safe_pdf(distx+handicap_mean, x); label="distx + handicap")
plot!(xs, x -> pdf(disty, x); label="disty")
title!("fair for t-tests")
P2 = plot(xs, x -> safe_pdf(distx+handicap_win_rate, x); label="distx + handicap")
plot!(xs, x -> pdf(disty, x); label="disty")
title!("fair for Mann-Whitney U-tests")
plot(P1, P2; size=(640, 200), legend=true)

# %%
ecdf_p_values(igx, 10, igx, 10; maxalpha=0.1, continuity_correction=false)

# %%
ecdf_p_values(igx, 30, igx, 30; maxalpha=0.1, continuity_correction=false)

# %%
ecdf_p_values(igx, 100, igx, 100; maxalpha=0.1, continuity_correction=false)

# %%
ecdf_p_values(igx, 10, igy, 10; maxalpha=0.1, continuity_correction=false)

# %%
ecdf_p_values(igx, 30, igy, 10; maxalpha=0.1, continuity_correction=false)

# %%
ecdf_p_values(igx, 10, igy, 30; maxalpha=0.1, continuity_correction=false)

# %%
ecdf_p_values(igx, 100, igy, 100; maxalpha=0.1, continuity_correction=false)

# %%
ecdf_p_values(igx, 300, igy, 100; maxalpha=0.1, continuity_correction=false)

# %%
ecdf_p_values(igx, 100, igy, 300; maxalpha=0.1, continuity_correction=false)

# %%
ecdf_p_values(igx, 300, igy, 300; maxalpha=0.1, continuity_correction=false)

# %%
ecdf_p_values(igx, 900, igy, 300; maxalpha=0.1, continuity_correction=false, niters=10^5)

# %%
ecdf_p_values(igx, 300, igy, 900; maxalpha=0.1, continuity_correction=false, niters=10^5)

# %%
ecdf_p_values(igx, 900, igy, 900; maxalpha=0.1, continuity_correction=false, niters=10^5)

# %%
ecdf_p_values(igx, 2700, igy, 900; maxalpha=0.1, continuity_correction=false, niters=10^5)

# %%
ecdf_p_values(igx, 2700, igy, 2700; maxalpha=0.1, continuity_correction=false, niters=10^5)

# %%
