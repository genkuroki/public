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
#     display_name: Julia 1.11.2
#     language: julia
#     name: julia-1.11
# ---

# %%
using Distributions
using Random
using StatsPlots
default(fmt=:png, titlefontsize=10, size=(500, 300))

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
using Distributions
using Roots

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

function confint_bm_p_roots(X, Y; α=0.05)
    f(p) = brunner_munzel_test(X, Y; p).pvalue - α
    find_zeros(f, -1, 2)
end

function aminamax(X, Y)
    xmin, xmax = extrema(X)
    ymin, ymax = extrema(Y)
    width = max(xmax, ymax) - min(xmin, ymin)
    xmin - ymax - max(0.1, 0.05width), xmax - ymin + max(0.1, 0.05width)
end

function tieshift(X, Y; p=1/2)
    f(a) = winrate(X, Y .+ a) - p
    amin, amax = aminamax(X, Y)
    find_zero(f, (amin, amax))
end

function confint_bm_tieshift(X, Y; p=1/2, α=0.05)
    f(a) = brunner_munzel_test(X, Y .+ a; p).pvalue - α
    amin, amax = aminamax(X, Y)
    find_zeros(f, amin, amax)
end

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

function confint_mw_tieshift(X, Y; α=0.05, correct=true)
    f(a) = pvalue_mann_whitney_u_test(X, Y .+ a; correct) - α
    amin, amax = aminamax(X, Y)
    find_zeros(f, amin, amax)
end

using HypothesisTests
X = randn(100)
Y = randn(100)
@show pvalue_mann_whitney_u_test(X, Y) pvalue(ApproximateMannWhitneyUTest(X, Y))
X = randn(10)
Y = randn(10)
@show pvalue_mann_whitney_u_test(X, Y) pvalue(ApproximateMannWhitneyUTest(X, Y));

# %%
function plot_brunner_munzel_pvalue_function(X, Y; α=0.05,
        as=range(-60, 60, 500), ps=range(0, 1, 400), size=(500, 400), kwargs...)
    f(a, p) = (y = pvalue_brunner_munzel(X, Y .+ a; p); y < α ? NaN : y)
    heatmap(as, ps, f; c=:turbo, cbar_title="P-value", clim=(0, 1))
    hline!([0.5]; label="", c=(α == 0 ? :white : :black), lw=0.5)
    vline!([0.0]; label="", c=(α == 0 ? :white : :black), lw=0.5)
    plot!(xtick=-100:10:100, ytick=0:0.05:1)
    plot!(xguide="a", yguide="p")
    title!("Brunner-Munzel " * (α == 0 ? "P-value function" : "$(100(1-α))% confidence region"))
    plot!(; size, kwargs...)
end

function plot_brunner_munzel_svalue_function(X, Y; α=0.05,
        as=range(-60, 60, 500), ps=range(0, 1, 400), size=(500, 400), kwargs...)
    f(a, p) = (y = pvalue_brunner_munzel(X, Y .+ a; p); y < α ? NaN : y)
    heatmap(as, ps, (a, p) -> -log2(f(a, p)); c=reverse(cgrad(:turbo)), 
        cbar_title="S-value = −log₂(P-value) [bit]", clim=(0, min(10, -log2(α))))
    hline!([0.5]; label="", c=(α == 0 ? :white : :black), lw=0.5)
    vline!([0.0]; label="", c=(α == 0 ? :white : :black), lw=0.5)
    plot!(xtick=-100:10:100, ytick=0:0.05:1)
    plot!(xguide="a", yguide="p")
    title!("Brunner-Munzel " * (α == 0 ? "S-value function" : "$(100(1-α))% confidence region"))
    plot!(; size, kwargs...)
end

# %%
X = [24, 24, 36, 36, 37, 38, 40, 44, 45, 48, 50, 53, 55, 57, 58, 67, 72, 75, 105, 127]
Y = [7, 7, 9, 10, 13, 13, 14, 17, 23, 24, 25, 27, 32, 35, 38, 40, 40, 41, 41, 42, 44, 45, 50, 53, 63, 74, 76, 84, 156, 727]

P = dotplot([X, Y]; label="", xtick=(1:2, ["X", "Y"]), msc=:auto, ms=3, ma=0.7, title="all data")
Q = dotplot([X, Y]; label="", xtick=(1:2, ["X", "Y"]), msc=:auto, ms=3, ma=0.7, title="data < 200", ylim=(-5, 180))
R = dotplot(vec(X .- Y'); xtick=(1:1, ["X − Y"]), label="", c=3, msc=:auto, ms=1.5, ma=0.7)
title!("differences of pairs")
plot(P, Q, R; layout=(1, 3), size=(600, 300))

# %%
r(x) = round.(x; sigdigits=3)

P_BM = pvalue_brunner_munzel(X, Y)
HL = hodges_lehmann(X, Y)
CI_BM = confint_bm_tieshift(X, Y)

P_W = pvalue_welch(X, Y)
DM = mean(X) - mean(Y)
CI_W = confint_welch(X, Y)

P_MW = pvalue_mann_whitney_u_test(X, Y)
CI_MW = confint_mw_tieshift(X, Y)

P_S = pvalue_student(X, Y)
CI_S = confint_student(X, Y)

println("Brunner-Munzel: ", "null P-value = ", r(P_BM), ",  point estimate = ", r(HL), ",  95% confidence interval = ", r.(CI_BM))
println("Welch t-test: ", "null P-value = ", r(P_W), ",  point estimate = ", r(DM), ",  95% confidence interval = ", r.(CI_W))
println("Mann-Whitney: ", "null P-value = ", r(P_MW), ",  point estimate = ", r(HL), ",  95% confidence interval = ", r.(CI_MW))
println("Student: ", "null P-value = ", r(P_S), ",  point estimate = ", r(DM), ",  95% confidence interval = ", r.(CI_S))

plot(a -> pvalue_brunner_munzel(X, Y .+ a), -100, 100; label="Brunner-Munzel", c=1)
vline!([HL]; label="Hodges-Lehmann", ls=:dot, c=1)
plot!(a -> pvalue_welch(X, Y .+ a); label="Welch", c=2, ls=:dash)
vline!([DM]; label="difference of means", ls=:dot, c=2)
plot!(a -> pvalue_mann_whitney_u_test(X, Y .+ a); label="Mann-Whitney", c=3, ls=:dashdot)
plot!(a -> pvalue_student(X, Y .+ a); label="Student", c=4, ls=:dashdotdot)
plot!(xtick=-100:20:100, ytick=0:0.05:1)
plot!(xguide="a", yguide="P-value")
title!("comparison between X and Y+a with p=1/2")
plot!(size=(600, 400))

# %%
for α in (0.05, 0.01, 0.0)
    plot_brunner_munzel_pvalue_function(X, Y; α) |> display
end

# %%
for α in (0.05, 0.01, 0.0)
    plot_brunner_munzel_svalue_function(X, Y; α) |> display
end

# %%
@show confint_bm_p_roots(X, Y .+ 30) .|> r;

# %%
@show confint_bm_tieshift(X, Y; p=0.6) .|> r;

# %%
X = [46, 51, 54, 54, 55, 61, 62, 63, 65, 65, 68, 70, 71, 71, 73, 74, 74, 74, 76, 76, 83, 84, 86, 86, 86, 87, 88, 88, 90, 91]
Y = [25, 26, 28, 28, 28, 35, 42, 47, 52, 60, 70, 91, 122, 145, 446]

P = dotplot([X, Y]; label="", xtick=(1:2, ["X", "Y"]), msc=:auto, ms=3, ma=0.7, title="all data")
Q = dotplot([X, Y]; label="", xtick=(1:2, ["X", "Y"]), msc=:auto, ms=3, ma=0.7, title="data < 200", ylim=(-5, 180))
R = dotplot(vec(X .- Y'); xtick=(1:1, ["X − Y"]), label="", c=3, msc=:auto, ms=1.5, ma=0.7)
title!("differences of pairs")
plot(P, Q, R; layout=(1, 3), size=(600, 300))

# %%
r(x) = round.(x; sigdigits=3)

P_BM = pvalue_brunner_munzel(X, Y)
HL = hodges_lehmann(X, Y)
CI_BM = confint_bm_tieshift(X, Y)

P_W = pvalue_welch(X, Y)
DM = mean(X) - mean(Y)
CI_W = confint_welch(X, Y)

P_MW = pvalue_mann_whitney_u_test(X, Y)
CI_MW = confint_mw_tieshift(X, Y)

P_S = pvalue_student(X, Y)
CI_S = confint_student(X, Y)

println("Brunner-Munzel: ", "null P-value = ", r(P_BM), ",  point estimate = ", r(HL), ",  95% confidence interval = ", r.(CI_BM))
println("Welch t-test: ", "null P-value = ", r(P_W), ",  point estimate = ", r(DM), ",  95% confidence interval = ", r.(CI_W))
println("Mann-Whitney: ", "null P-value = ", r(P_MW), ",  point estimate = ", r(HL), ",  95% confidence interval = ", r.(CI_MW))
println("Student: ", "null P-value = ", r(P_S), ",  point estimate = ", r(DM), ",  95% confidence interval = ", r.(CI_S))

plot(a -> pvalue_brunner_munzel(X, Y .+ a), -100, 100; label="Brunner-Munzel", c=1)
vline!([HL]; label="Hodges-Lehmann", ls=:dot, c=1)
plot!(a -> pvalue_welch(X, Y .+ a); label="Welch", c=2, ls=:dash)
vline!([DM]; label="difference of means", ls=:dot, c=2)
plot!(a -> pvalue_mann_whitney_u_test(X, Y .+ a); label="Mann-Whitney", c=3, ls=:dashdot)
plot!(a -> pvalue_student(X, Y .+ a); label="Student", c=4, ls=:dashdotdot)
plot!(xtick=-100:20:100, ytick=0:0.05:1)
plot!(xguide="a", yguide="P-value")
title!("comparison between X and Y+a with p=1/2")
plot!(size=(600, 400))

# %%
@show pvalue_welch(X, Y .+ 30) .|> r
@show pvalue_student(X, Y .+ 30) .|> r;

# %%
for α in (0.05, 0.01, 0.0)
    plot_brunner_munzel_pvalue_function(X, Y; α) |> display
end

# %%
for α in (0.05, 0.01, 0.0)
    plot_brunner_munzel_svalue_function(X, Y; α) |> display
end

# %% [markdown]
# https://x.com/genkuroki/status/1868839785118548371

# %% [markdown]
# <img src="IMG_7529.jpeg">

# %%
Y = [fill(0, 94); fill(1, 4); fill(2, 4); 3]
X = [fill(0, 65); fill(1, 14); fill(2, 10); fill(3, 6); fill(4, 4); 5; 5; 6; 6; 7; 7; 8; 9; 10; 11; 11; 19]
@show X Y
@show length(Y) length(X)

P = dotplot([X, Y]; label="", xtick=(1:2, ["X", "Y"]), msc=:auto, ms=1.5, ma=0.7, title="all data")
#Q = dotplot([X, Y]; label="", xtick=(1:2, ["X", "Y"]), msc=:auto, ms=1.5, ma=0.7, title="data < 200", ylim=(-5, 180))
R = dotplot(vec(X .- Y'); xtick=(1:1, ["X − Y"]), label="", c=3, msc=:auto, ms=0.75, ma=0.1)
title!("differences of pairs")
#plot(P, Q, R; layout=(1, 3), size=(600, 300))
plot(P, R; layout=@layout[a{0.33w} b], size=(1000, 300))

# %%
@show confint_bm_p_roots(X, Y) .|> r;

# %%
@show confint_bm_tieshift(X, Y) .|> r;

# %%
r(x) = round.(x; sigdigits=3)

P_BM = pvalue_brunner_munzel(X, Y)
HL = hodges_lehmann(X, Y)
CI_BM = confint_bm_tieshift(X, Y)

P_W = pvalue_welch(X, Y)
DM = mean(X) - mean(Y)
CI_W = confint_welch(X, Y)

P_MW = pvalue_mann_whitney_u_test(X, Y)
CI_MW = confint_mw_tieshift(X, Y)

P_S = pvalue_student(X, Y)
CI_S = confint_student(X, Y)

println("Brunner-Munzel: ", "null P-value = ", r(P_BM), ",  point estimate = ", r(HL), ",  95% confidence interval = ", r.(CI_BM))
println("Welch t-test: ", "null P-value = ", r(P_W), ",  point estimate = ", r(DM), ",  95% confidence interval = ", r.(CI_W))
println("Mann-Whitney: ", "null P-value = ", r(P_MW), ",  point estimate = ", r(HL), ",  95% confidence interval = ", r.(CI_MW))
println("Student: ", "null P-value = ", r(P_S), ",  point estimate = ", r(DM), ",  95% confidence interval = ", r.(CI_S))

as = range(-1, 3, 401)
plot(as, a -> pvalue_brunner_munzel(X, Y .+ a); label="Brunner-Munzel", c=1)
vline!([HL]; label="Hodges-Lehmann", ls=:dot, c=1)
#plot!(a -> pvalue_welch(X, Y .+ a); label="Welch", c=2, ls=:dash)
vline!([DM]; label="difference of means", ls=:dot, c=2)
plot!(a -> pvalue_mann_whitney_u_test(X, Y .+ a); label="Mann-Whitney", c=3, ls=:dashdot)
#plot!(a -> pvalue_student(X, Y .+ a); label="Student", c=4, ls=:dashdotdot)
#plot!(xtick=-100:20:100, ytick=0:0.05:1)
plot!(xguide="a", yguide="P-value")
title!("comparison between X and Y+a with p=1/2")
plot!(size=(600, 400))

# %%
for α in (0.05, 0.01, 0.0)
    plot_brunner_munzel_svalue_function(X, Y; α, as=range(-1, 3, 401), xtick=:auto) |> display
end

# %%
log(100)

# %%
