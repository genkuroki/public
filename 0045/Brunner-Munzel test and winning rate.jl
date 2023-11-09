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
using Random
Random.seed!(4649373)

# %%
using Distributions
using QuadGK
using Roots

function winningrate(
        distx::ContinuousUnivariateDistribution,
        disty::ContinuousUnivariateDistribution
    )
    quadgk(y -> cdf(distx, y)*pdf(disty, y), extrema(disty)...)[1]
end

function tieshift(
        distx::ContinuousUnivariateDistribution,
        disty::ContinuousUnivariateDistribution;
        p = 1/2
    )
    find_zero(0.0) do a
        winningrate(distx + a, disty) - p
    end
end

@show distx = Exponential(1)
@show disty = Exponential(1.1)
@show distx = distx - mean(distx)
@show disty = disty - mean(disty)
@show winningrate(distx, disty)
@show a = tieshift(distx, disty)
@show winningrate(distx + a, disty);9

# %%
@show distx = Exponential(1)
@show disty = Exponential(1.1)
@show distx = distx - median(distx)
@show disty = disty - median(disty)
@show winningrate(distx, disty)
@show a = tieshift(distx, disty)
@show winningrate(distx + a, disty);

# %%
@show distx = LogNormal()
@show disty = 1.1LogNormal()
@show distx = distx - mean(distx)
@show disty = disty - mean(disty)
@show winningrate(distx, disty)
@show a = tieshift(distx, disty)
@show winningrate(distx + a, disty);

# %%
@show distx = Exponential(1)
@show disty = Exponential(1.1)
@show distx = distx - median(distx)
@show disty = disty - median(disty)
@show winningrate(distx, disty)
@show a = tieshift(distx, disty)
@show winningrate(distx + a, disty);

# %%
c(x...) = [x...]

@show sigma1 = c(1, 2, 3, 4, 4, 4, 4)
@show sigma2 = c(4, 4, 4, 4, 3, 2, 1)
@show sigmas1 = @. sqrt(log((1 + sqrt((4 * sigma1^2)/25 + 1))/2))
@show sigmas2 = @. sqrt(log((1 + sqrt((4 * sigma2^2)/25 + 1))/2))
@show mu1 = log(5)
@show mu2 = log(5)
@show lognormals1 = LogNormal.(mu1, sigmas1)
@show lognormals2 = LogNormal.(mu2, sigmas2)
@show median.(lognormals1)
@show median.(lognormals2)
@show winningrate.(lognormals1, lognormals2)
@show tieshift.(lognormals1, lognormals2);

# %%
using Distributions

safediv(x, y) = x == 0 ? x : isinf(y) ? zero(y) : x/y

"""
    h_brunner_munzel(x, y)

この函数は, x < y のとき 1.0 を, x = y のとき 0.5 を, それら以外のとき 0.0 返す.
"""
h_brunner_munzel(x, y) = (x < y) + (x == y)/2

@doc raw"""
    brunner_munzel_test(X, Y; p = 1/2)

この函数は数値のベクトルのデータ `X`, `Y` について, Brunner-Munzel検定関係の統計量達を計算する. 詳細は以下の通り.

函数 $H(x, y)$ と $\hat{p}$, $H^x_i$, $H^y_j$, $\bar{H}^x$, $\bar{H}^y$ を次のように定める:

```math
\begin{aligned}
&
m = \mathrm{length}(X), \quad
n = \mathrm{length}(Y), \quad
x_i = X[i], \quad
y_j = Y[j],
\\ &
\hat{p} = \frac{1}{mn}\sum_{i=1}^m \sum_{j=1}^n H(x_i, y_j),
\\ &
H(x, y) = \begin{cases} 1 & (x < y) \\ 1/2 & (x = y) \\ 0 & (x > y), \end{cases}
\\ &
H^x_i = \sum_{j=1}^n H(y_j, x_i), \quad
H^y_j = \sum_{i=1}^m H(x_i, y_j),
\\ &
\bar{H}^x = \frac{1}{m} \sum_{i=1}^m H^x_i = n - n\hat{p},
\\ &
\bar{H}^y = \frac{1}{n} \sum_{j=1}^n H^y_j = m\hat{p}.
\end{aligned}
```

この函数は以下達の named tuple で返す:

```math
\begin{aligned}
&
\mathrm{phat} = 
\hat{p} = \frac{\bar{H}^y - \bar{H}^x + n}{m + n},
\\ &
\mathrm{sx2} =
\hat{\sigma}_x^2 = \frac{1}{n^2}\frac{1}{m-1}\sum_{i=1}^m (H^x_i - \bar{H}^x)^2,
\\ &
\mathrm{sy2} =
\hat{\sigma}_y^2 = \frac{1}{m^2}\frac{1}{n-1}\sum_{j=1}^n (H^y_j - \bar{H}^y)^2,
\\ &
\mathrm{sehat} = 
\widehat{\mathrm{se}} = \sqrt{\frac{\hat{\sigma}_x^2}{m} + \frac{\hat{\sigma}_y^2}{n}}, 
\\ &
\mathrm{tvalue} = t = \frac{\hat{p} - p}{\widehat{\mathrm{se}}},
\\ &
\mathrm{df} =
\nu = 
\frac
{\left(\hat{\sigma}_x^2/m + \hat{\sigma}_y^2/n\right)^2}
{
\dfrac{\left(\hat{\sigma}_x^2/m\right)^2}{m-1} +
\dfrac{\left(\hat{\sigma}_y^2/n\right)^2}{n-1}
},
\\ &
\mathrm{pvalue} =
2\mathrm{ccdf}(\mathrm{TDist}(\nu), |t|),
\\ &
\mathrm{p} = p.
\end{aligned}
```
"""
function brunner_munzel_test(X, Y; p=1/2)
    m, n = length(X), length(Y)
    phat = mean(h_brunner_munzel(x, y) for x in X, y in Y)
    Hbarx = n*(1 - phat)
    Hbary = m*phat
    sx2 = 1/n^2 * 1/(m-1) * sum(x -> (sum(h_brunner_munzel(y, x) for y in Y) - Hbarx)^2, X)
    sy2 = 1/m^2 * 1/(n-1) * sum(y -> (sum(h_brunner_munzel(x, y) for x in X) - Hbary)^2, Y)
    sehat = √(sx2/m + sy2/n)
    tvalue = (phat - p)/sehat
    df = (sx2/m + sy2/n)^2 / ((sx2/m)^2/(m-1) + (sy2/n)^2/(n-1))
    pvalue = 2ccdf(TDist(df), abs(tvalue))
    (; phat, sehat, tvalue, df, pvalue, p)
end

@doc raw"""
    pvalue_brunner_munzel_test(X, Y; p = 1/2)

この函数はBrunner-Munzel検定のP値 `pvalue` を返す.
"""
pvalue_brunner_munzel_test(X, Y; p=1/2) = brunner_munzel_test(X, Y; p).pvalue

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

# %%
using HypothesisTests

X = rand(30, 20)
Y = rand(30, 20)
pval1 = [HypothesisTests.pvalue(MannWhitneyUTest(x, y)) for (x, y) in zip(eachcol(X), eachcol(Y))]
pval2 = [pvalue_mann_whitney_u_test(x, y) for (x, y) in zip(eachcol(X), eachcol(Y))]
pval3 = [pvalue_mann_whitney_u_test(x, y; correct=false) for (x, y) in zip(eachcol(X), eachcol(Y))]
[pval1 pval2 pval3 pval1-pval2 pval1-pval3]

# %%
using HypothesisTests

X = rand(1000, 20)
Y = rand(1000, 20)
pval1 = [HypothesisTests.pvalue(MannWhitneyUTest(x, y)) for (x, y) in zip(eachcol(X), eachcol(Y))]
pval2 = [pvalue_mann_whitney_u_test(x, y) for (x, y) in zip(eachcol(X), eachcol(Y))]
pval3 = [pvalue_mann_whitney_u_test(x, y; correct=false) for (x, y) in zip(eachcol(X), eachcol(Y))]
[pval1 pval2 pval3 pval1-pval2 pval1-pval3]

# %%
X = rand(100)
Y = rand(100)
@time @show brunner_munzel_test(X, Y)
@time @show mann_whitney_u_test(X, Y);

# %%
distx = Exponential(1)
disty = Exponential(3)
distx = distx - mean(distx)
disty = disty - mean(disty)
a = minimum(quantile.((distx, disty), 0.01))
b = maximum(quantile.((distx, disty), 0.99))
a, b

# %%
using HypothesisTests
ECDF(A, x) = count(≤(x), A)/length(A)
using StatsPlots
default(fmt=:png)
using Random

function plot_sim_bm_mw(;
        distx = Exponential(1),
        disty = Exponential(3),
        xlim = :auto,
        ylim = :auto,
        m = 100,
        n = 100,
        L = 10^6,
        shifttype = :median,
        α = 0.05,
        legend = :bottomright,
        correct = true,
    )
    @show distx
    @show disty
    @show tiesh = tieshift(distx, disty)
    @show mediansh = median(disty) - median(distx)
    @show meansh = mean(disty) - mean(distx)
    println()
    @show shifttype
    if shifttype == :tie
        @show distx = distx + tiesh
    elseif shifttype == :median
        @show distx = distx + mediansh
    elseif shifttype == :mean
        @show distx = distx + meansh
    end
    @show winningrate(distx, disty)
    @show median(distx) median(disty)
    @show mean(distx) mean(disty)

    if xlim == :auto
        a = minimum(quantile.((distx, disty), 0.01))
        b = maximum(quantile.((distx, disty), 0.99))
        l = b - a
        a = a - l/10
        b = b + l/10
    else
        a, b = xlim
    end
    P1 = plot(distx, a, b; label="distx")
    plot!(disty, a, b; label="disty", ls=:dash)
    plot!(; ylim)
    title!("m = $m,  n = $n"; titlefontsize=10)

    @show m
    @show n
    println()

    pval_bm = zeros(L)
    pval_mw = zeros(L)

    nth = Threads.nthreads()
    Xtmp = [zeros(m) for _ in 1:nth]
    Ytmp = [zeros(n) for _ in 1:nth]
    
    @time Threads.@threads for i in 1:L
        tid = Threads.threadid()
        X = rand!(distx, Xtmp[tid])
        Y = rand!(disty, Ytmp[tid])
        pval_bm[i] = pvalue_brunner_munzel_test(X, Y)
        pval_mw[i] = pvalue_mann_whitney_u_test(X, Y; correct)
    end

    @show α
    @show ECDF(pval_bm, α)
    @show ECDF(pval_mw, α)

    _tick = [0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1]
    xtick = ytick = (_tick, string.(_tick))

    αs = exp.(range(log(0.002), log(1), 1000))
    P2 = plot(αs, α -> ECDF(pval_bm, α); label="Brunner-Munzel test")
    plot!(αs, α -> ECDF(pval_mw, α); label="Mann-Whitney U-test", ls=:dash)
    plot!(αs, identity; label="", ls=:dot, c=:gray)
    plot!(αs, x->0.8x; label="", ls=:dot, c=:gray)
    plot!(αs, x->1.2x; label="", ls=:dot, c=:gray)
    plot!(; xscale=:log10, yscale=:log10, xtick, ytick)
    plot!(; legend)
    
    plot(P1, P2; size=(400, 570), layout=@layout [a{0.3h}; b])
    plot!(leftmargin=4Plots.mm)
end

# %%
plot_sim_bm_mw(; m=100, n=100, shifttype=:tie)

# %%
plot_sim_bm_mw(; m=100, n=100, shifttype=:median)

# %%
plot_sim_bm_mw(; m=100, n=100, shifttype=:mean)

# %%
plot_sim_bm_mw(; m=200, n=200, shifttype=:median, L=10^5)

# %%
plot_sim_bm_mw(; m=500, n=500, shifttype=:median, L=5*10^4)

# %%
plot_sim_bm_mw(; m=1000, n=1000, shifttype=:median, L=10^4)

# %%
"""
gammadist(σ, β)

returns the gamma distribution with standard deviation `σ` and skewness `β`.
"""
gammadist(σ, β) = Gamma(4/β^2, β*σ/2)

# %%
gam = gammadist.(1:5, 2:2:10)
std.(gam), skewness.(gam)

# %%
[(β, shape(gammadist(1, β)), gammadist(1, β)) for β in 0.2:0.2:3]

# %%
P = plot()
for β in [0.4; 1:3]
    dist = gammadist(1, β)
    @eval @show (mean(gammadist(1, $β)), median(gammadist(1, $β)), shape(gammadist(1, $β)))
    dist = dist - median(dist)
    plot!(dist, -3, 3; label="$β", ls=:auto)
    plot!(ylim=(-0.05, 2.05))
end
plot!()

# %%
distx = gammadist(1, 2.5)
disty = gammadist(1, 3)
xlim = (-1.0, 1.5)
ylim = (-0.2, 5.2)
@show std(distx) std(disty)
@show skewness(distx) skewness(disty)
plot_sim_bm_mw(; distx, disty, xlim, ylim, m=100, n=100, shifttype=:median, L=10^5)

# %%
distx = gammadist(1, 2.5)
disty = gammadist(1, 3)
xlim = (-1.0, 1.5)
ylim = (-0.2, 5.2)
@show std(distx) std(disty)
@show skewness(distx) skewness(disty)
plot_sim_bm_mw(; distx, disty, xlim, ylim, m=100, n=100, shifttype=:tie, L=10^5)

# %%
distx = gammadist(1, 2.5)
disty = gammadist(1.5, 3)
xlim = (-1.0, 1.5)
ylim = (-0.2, 5.2)
@show std(distx) std(disty)
@show skewness(distx) skewness(disty)
plot_sim_bm_mw(; distx, disty, xlim, ylim, m=100, n=100, shifttype=:median, L=10^5)

# %%
distx = gammadist(1, 2.5)
disty = gammadist(1.5, 3)
xlim = (-1.0, 1.5)
ylim = (-0.2, 5.2)
@show std(distx) std(disty)
@show skewness(distx) skewness(disty)
plot_sim_bm_mw(; distx, disty, xlim, ylim, m=100, n=100, shifttype=:tie, L=10^5)

# %%
distx = gammadist(1, 2.5)
disty = gammadist(2, 3)
xlim = (-0.5, 1.5)
ylim = (-0.2, 5.2)
@show std(distx) std(disty)
@show skewness(distx) skewness(disty)
plot_sim_bm_mw(; distx, disty, xlim, ylim, m=100, n=100, shifttype=:median, L=10^5)

# %%
distx = gammadist(1, 2.5)
disty = gammadist(2, 3)
xlim = (-0.5, 1.5)
ylim = (-0.2, 5.2)
@show std(distx) std(disty)
@show skewness(distx) skewness(disty)
plot_sim_bm_mw(; distx, disty, xlim, ylim, m=100, n=100, shifttype=:tie, L=10^5)

# %%
distx = gammadist(1.5, 2.5)
disty = gammadist(1, 3)
xlim = (-1, 1.5)
ylim = (-0.2, 5.2)
@show std(distx) std(disty)
@show skewness(distx) skewness(disty)
plot_sim_bm_mw(; distx, disty, xlim, ylim, m=100, n=100, shifttype=:median, L=10^5)

# %%
distx = gammadist(1.5, 2.5)
disty = gammadist(1, 3)
xlim = (-1, 1.5)
ylim = (-0.2, 5.2)
@show std(distx) std(disty)
@show skewness(distx) skewness(disty)
plot_sim_bm_mw(; distx, disty, xlim, ylim, m=100, n=100, shifttype=:tie, L=10^5)

# %%
distx = gammadist(2, 2.5)
disty = gammadist(1, 3)
xlim = (-1, 1.5)
ylim = (-0.2, 5.2)
@show std(distx) std(disty)
@show skewness(distx) skewness(disty)
plot_sim_bm_mw(; distx, disty, xlim, ylim, m=100, n=100, shifttype=:median, L=10^5)

# %%
distx = gammadist(2, 2.5)
disty = gammadist(1, 3)
xlim = (-1.0, 1.5)
ylim = (-0.2, 5.2)
@show std(distx) std(disty)
@show skewness(distx) skewness(disty)
plot_sim_bm_mw(; distx, disty, xlim, ylim, m=100, n=100, shifttype=:tie, L=10^5)

# %%
distx = gammadist(3, 2.5)
disty = gammadist(1, 3)
xlim = (-1.5, 1.5)
ylim = (-0.2, 5.2)
@show std(distx) std(disty)
@show skewness(distx) skewness(disty)
plot_sim_bm_mw(; distx, disty, xlim, ylim, m=100, n=100, shifttype=:median, L=10^5)

# %%
distx = gammadist(3, 2.5)
disty = gammadist(1, 3)
xlim = (-1.5, 1.5)
ylim = (-0.2, 5.2)
@show std(distx) std(disty)
@show skewness(distx) skewness(disty)
plot_sim_bm_mw(; distx, disty, xlim, ylim, m=100, n=100, shifttype=:tie, L=10^5)

# %%
distx = gammadist(1, 3)
disty = gammadist(2, 3)
xlim = (-1.0, 2.0)
ylim = (-0.2, 5.2)
@show std(distx) std(disty)
@show skewness(distx) skewness(disty)
plot_sim_bm_mw(; distx, disty, xlim, ylim, m=100, n=100, shifttype=:median, L=10^6, correct=false)

# %%
distx = gammadist(1, 3)
disty = gammadist(2, 3)
xlim = (-1.0, 2.0)
ylim = (-0.2, 5.2)
@show std(distx) std(disty)
@show skewness(distx) skewness(disty)
plot_sim_bm_mw(; distx, disty, xlim, ylim, m=100, n=100, shifttype=:median, L=10^6, correct=true)

# %%
function inversegammadist(σ, β)
    α = 3 + 8/β^2 + 4/β * √(1 + 4/β^2)
    θ = σ * (α - 1) * √(α - 2)
    InverseGamma(α, θ)
end

# %%
igam = inversegammadist.(1:5, 2:2:10)
std.(igam), skewness.(igam)

# %%
@show igam = inversegammadist.(1, Inf)
std.(igam), skewness.(igam)

# %%
P = plot()
for β in [0.4; 1:3]
    dist = inversegammadist(1, β)
    @show (β, mean(dist), median(dist), shape(dist))
    dist = dist - median(dist)
    plot!(dist, -3, 6; label="$β", ls=:auto)
end
plot!()

# %%
distx = inversegammadist(1, 1)
disty = inversegammadist(1, 10)
@show std(distx) std(disty)
@show skewness(distx) skewness(disty)
plot_sim_bm_mw(; distx, disty, m=100, n=100, shifttype=:median, L=10^5)

# %%
distx = inversegammadist(1, 1)
disty = inversegammadist(1, 10)
@show std(distx) std(disty)
@show skewness(distx) skewness(disty)
plot_sim_bm_mw(; distx, disty, m=100, n=100, shifttype=:mean, L=10^5)

# %%
distx = inversegammadist(1, 1)
disty = inversegammadist(1, 10)
@show std(distx) std(disty)
@show skewness(distx) skewness(disty)
plot_sim_bm_mw(; distx, disty, m=100, n=100, shifttype=:tie, L=10^5)

# %%
distx = inversegammadist(1, 1)
disty = inversegammadist(1, 10)
@show std(distx) std(disty)
@show skewness(distx) skewness(disty)
plot_sim_bm_mw(; distx, disty, m=50, n=100, shifttype=:median, L=10^5)

# %%
distx = inversegammadist(1, 1)
disty = inversegammadist(1, 10)
@show std(distx) std(disty)
@show skewness(distx) skewness(disty)
plot_sim_bm_mw(; distx, disty, m=50, n=100, shifttype=:mean, L=10^5)

# %%
distx = inversegammadist(1, 1)
disty = inversegammadist(1, 10)
@show std(distx) std(disty)
@show skewness(distx) skewness(disty)
plot_sim_bm_mw(; distx, disty, m=50, n=100, shifttype=:tie, L=10^5)

# %%
distx = inversegammadist(2, 1)
disty = inversegammadist(1, 10)
@show std(distx) std(disty)
@show skewness(distx) skewness(disty)
plot_sim_bm_mw(; distx, disty, m=50, n=100, shifttype=:median, L=10^5)

# %%
distx = inversegammadist(2, 1)
disty = inversegammadist(1, 10)
@show std(distx) std(disty)
@show skewness(distx) skewness(disty)
plot_sim_bm_mw(; distx, disty, m=50, n=100, shifttype=:mean, L=10^5)

# %%
distx = inversegammadist(2, 1)
disty = inversegammadist(1, 10)
@show std(distx) std(disty)
@show skewness(distx) skewness(disty)
plot_sim_bm_mw(; distx, disty, m=50, n=100, shifttype=:tie, L=10^5)

# %%
