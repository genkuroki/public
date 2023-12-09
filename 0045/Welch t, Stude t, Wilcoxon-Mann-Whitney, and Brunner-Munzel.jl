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

# %% [markdown]
# See https://github.com/genkuroki/public/blob/main/0045/Welch%20t%2C%20Student%2C%20Wilcoxon-Mann-Whitney%2C%20and%20Brunner-Munzel.ipynb

# %%
using Random
Random.seed!(4649373)

using Distributions
using QuadGK
using Roots
using StatsPlots
default(fmt=:png)

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
@show disty = Exponential(1.5)
@show winningrate(distx, disty)
println()
@show tiesh = tieshift(distx, disty)
@show mediansh = median(disty) - median(distx)
@show meansh = mean(disty) - mean(distx)
println()
@show winningrate(distx + tiesh, disty)
@show winningrate(distx + mediansh, disty)
@show winningrate(distx + meansh, disty)
println()
@show median(distx + tiesh)
@show median(distx + mediansh)
@show median(distx + meansh)
@show median(disty)
println()
@show mean(distx + tiesh)
@show mean(distx + mediansh)
@show mean(distx + meansh)
@show mean(disty)
;

# %%
using Distributions

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
\bar{H}^x = \frac{1}{m} \sum_{i=1}^m H^x_i = n - n\hat{p}, \quad
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
    phat = mean((x < y) + (x == y)/2 for x in X, y in Y)
    Hbarx = n*(1 - phat)
    Hbary = m*phat
    sx2 = 1/n^2 * 1/(m-1) * sum(x -> (sum((y < x) + (y == x)/2 for y in Y) - Hbarx)^2, X)
    sy2 = 1/m^2 * 1/(n-1) * sum(y -> (sum((x < y) + (x == y)/2 for x in X) - Hbary)^2, Y)
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

using RCall
@rimport lawstat
X = rand(100)
Y = rand(100)
pvalue_brunner_munzel_test(X, Y), rcopy(lawstat.brunner_munzel_test(X, Y))[:p_value]

# %%
@doc brunner_munzel_test

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

using HypothesisTests
X = randn(100)
Y = randn(100)
pvalue_mann_whitney_u_test(X, Y), pvalue(ApproximateMannWhitneyUTest(X, Y))

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
"""
gammadist(σ, β)

returns the gamma distribution with standard deviation `σ` and skewness `β`.
"""
gammadist(σ, β) = Gamma(4/β^2, β*σ/2)

gam = gammadist.(1:5, 2:2:10)
std.(gam), skewness.(gam), kurtosis.(gam)

# %%
@doc gammadist

# %%
"""
inversegammadist(σ, β)

returns the inverse gamma distribution with standard deviation `σ` and skewness `β`.
"""
function inversegammadist(σ, β)
    γ = 1/β
    α = 3 + 8γ^2 + 4γ * √(1 + 4γ^2)
    θ = σ * (α - 1) * √(α - 2)
    InverseGamma(α, θ)
end

igam = inversegammadist.(1:5, 2:2:10)
std.(igam), skewness.(igam), kurtosis.(igam)

# %%
@doc inversegammadist

# %%
function prc(s, x, α)
    if x < α/2
        printstyled(s, x; color=:blue, bold=true)
    elseif x < α/1.3
        printstyled(s, x; color=:blue)
    elseif x < 1.3α
        print(s, x)
    elseif x < 2α
        printstyled(s, x; color=:red)
    else
        printstyled(s, x; color=:red, bold=true)
    end
    println()
end

for k in 0.01:0.01:0.11
    prc("hoge: ", k, 0.05)
end

# %%
using Printf
undefvector(T::Type, m) = Vector{T}(undef, m)
undefvector(m) = undefvector(Float64, m)
ECDF(A, x) = count(≤(x), A)/length(A)
rd(x) = round(x; digits=1)

function plot_sim(;
        distx = inversegammadist(2, 5),
        disty = Normal(0, 1),
        xlim = :auto,
        ylim = :auto,
        m = 25,
        n = 100,
        L = 10^6,
        shifttype = :auto,
        α = 0.05,
        legend = :bottomright,
        correct = true,
        verbose = false,
    )
    
    tiesh = tieshift(distx, disty)
    mediansh = median(disty) - median(distx)
    meansh = mean(disty) - mean(distx)
    if verbose
        @show tiesh mediansh meansh
    end
    println()
    
    distx_tie = distx + tiesh
    distx_median = distx + mediansh
    distx_mean = distx + meansh
    if shifttype == :tie
        distx = distx_tie
    elseif shifttype == :median
        distx = distx_median
    elseif shifttype == :mean || shifttype == :auto
        distx = distx_mean
    end
    if verbose
        @show shifttype
        @show distx
        @show winningrate(distx, disty)
        @show median(distx) median(disty)
        @show mean(distx) mean(disty)
        println()
    end

    if xlim == :auto
        a = minimum(quantile.((distx_mean, distx_median, distx_tie, disty), 0.01))
        b = maximum(quantile.((distx_mean, distx_median, distx_tie, disty), 0.99))
        l = b - a
        a = a - l/10
        b = b + l/10
    else
        a, b = xlim
    end
    P1 = plot()
    if shifttype == :auto
        plot!(distx_mean, a, b; label="distx_mean", c=:blue)
        plot!(distx_tie, a, b; label="distx_tie", c=:darkcyan, ls=:dash)
    else
        plot!(distx, a, b; label="distx", c=:blue)
    end
    plot!(disty, a, b; label="disty", ls=:dashdot, c=:red)
    plot!(; ylim)
    title!("m = $m,  n = $n,  Niters = $L"; titlefontsize=10)

    pval_we = undefvector(L)
    pval_st = undefvector(L)
    tval_we = undefvector(L)
    tval_st = undefvector(L)
    sehat2_we = undefvector(L)
    sehat2_st = undefvector(L)
    
    pval_bm = undefvector(L)
    pval_mw = undefvector(L)
    tval_bm = undefvector(L)
    tval_mw = undefvector(L)
    sehat2_bm = undefvector(L)
    sehat2_mw = undefvector(L)
    phat = undefvector(L)

    nth = Threads.nthreads()
    Xtmp = [undefvector(m) for _ in 1:nth]
    Ytmp = [undefvector(n) for _ in 1:nth]
    
    @time Threads.@threads for i in 1:L
        tid = Threads.threadid()
        X = rand!(distx, Xtmp[tid])
        Y = rand!(disty, Ytmp[tid])
        we = welch_t_test(X, Y)
        st = student_t_test(X, Y)
        pval_we[i], tval_we[i], sehat2_we[i] = we.pvalue, we.tvalue, we.sehat2
        pval_st[i], tval_st[i], sehat2_st[i] = st.pvalue, st.tvalue, st.sehat2
        if shifttype == :auto
            @. X = X - meansh + tiesh
        end
        bm = brunner_munzel_test(X, Y)
        mw = mann_whitney_u_test(X, Y; correct)
        pval_bm[i], tval_bm[i], sehat2_bm[i] = bm.pvalue, bm.tvalue, bm.sehat^2
        pval_mw[i], tval_mw[i], sehat2_mw[i] = mw.pvalue, mw.zvalue, mw.sehat^2
        phat[i] = bm.phat
    end
    println()
    @printf "%6s  %6s  %10s  %10s\n" "" "std" "skewness" "kurtosis"
    println("-"^40)
    @printf "%6s  %6.1f  %10.1f  %10.1f\n" "distx" std(distx) skewness(distx) kurtosis(distx)
    @printf "%6s  %6.1f  %10.1f  %10.1f\n" "disty" std(disty) skewness(disty) kurtosis(disty)
    println()
    
    e_we = ECDF(pval_we, α)
    e_st = ECDF(pval_st, α)
    e_bm = ECDF(pval_bm, α)
    e_mw = ECDF(pval_mw, α)
    println("true α-error rate for α = $α")
    prc("  Welch t-test:        ", e_we, α)
    prc("  Student t-test:      ", e_st, α)
    prc("  Brunner-Munzel test: ", e_bm, α)
    prc("  WMW test:            ", e_mw, α)
    println()

    _tick = [0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1]
    xtick = ytick = (_tick, string.(_tick))
    αs = exp.(range(log(_tick[begin]), log(_tick[end]), 1000))
    P2a = plot()
    plot!(αs, α -> ECDF(pval_we, α); label="Welch ($(rd(100e_we))%)", c=1)
    plot!(αs, α -> ECDF(pval_st, α); label="Student ($(rd(100e_st))%)", ls=:dash, c=2)
    plot!(αs, α -> ECDF(pval_bm, α); label="BM ($(rd(100e_bm))%)", ls=:dashdot, c=3)
    plot!(αs, α -> ECDF(pval_mw, α); label="WMW ($(rd(100e_mw))%)", ls=:dashdotdot, c=4)
    scatter!([α], [e_we]; ms=3, msc=:auto, label="", c=1)
    scatter!([α], [e_st]; ms=3, msc=:auto, label="", c=2)
    scatter!([α], [e_bm]; ms=3, msc=:auto, label="", c=3)
    scatter!([α], [e_mw]; ms=3, msc=:auto, label="", c=4)
    plot!(αs, identity; label="", ls=:dot, c=:gray)
    plot!(αs, x->0.8x; label="", ls=:dot, c=:gray)
    plot!(αs, x->1.2x; label="", ls=:dot, c=:gray)
    plot!(; xscale=:log10, yscale=:log10, xtick, ytick)
    plot!(; xguide="α", yguide="probability of P-value ≤ α")
    plot!(; legend)
    
    _tick = round.([0.2, 0.4, 1, 2, 4] * α; sigdigits=2)
    xtick = ytick = (_tick, string.(_tick))
    αs = exp.(range(log(_tick[begin]), log(_tick[end]), 1000))
    P2b = plot()
    plot!(αs, α -> ECDF(pval_we, α); label="Welch ($(rd(100e_we))%)", c=1)
    plot!(αs, α -> ECDF(pval_st, α); label="Student ($(rd(100e_st))%)", ls=:dash, c=2)
    plot!(αs, α -> ECDF(pval_bm, α); label="BM ($(rd(100e_bm))%)", ls=:dashdot, c=3)
    plot!(αs, α -> ECDF(pval_mw, α); label="WMW ($(rd(100e_mw))%)", ls=:dashdotdot, c=4)
    scatter!([α], [e_we]; ms=3, msc=:auto, label="", c=1)
    scatter!([α], [e_st]; ms=3, msc=:auto, label="", c=2)
    scatter!([α], [e_bm]; ms=3, msc=:auto, label="", c=3)
    scatter!([α], [e_mw]; ms=3, msc=:auto, label="", c=4)
    plot!(αs, identity; label="", ls=:dot, c=:gray)
    plot!(αs, x->0.8x; label="", ls=:dot, c=:gray)
    plot!(αs, x->1.2x; label="", ls=:dot, c=:gray)
    plot!(; xscale=:log10, yscale=:log10, xtick, ytick)
    plot!(; xguide="α", yguide="probability of P-value ≤ α")
    plot!(; legend)

    P3 = plot()
    stephist!(tval_we; norm=true, label="We")
    stephist!(tval_st; norm=true, label="St", ls=:dash)
    plot!(Normal(), -5, 5; label="", ls=:dot, c=:gray)
    vline!([0.0]; label="", ls=:dot, c=:gray)
    plot!(; yscale=:log10, ylim=(1e-4, 1))
    title!("t-values"; titlefontsize=10)

    se2 = var(distx)/m + var(disty)/n
    if verbose
        @show (mean(sehat2_we)/se2, std(sehat2_we)/se2)
        @show (mean(sehat2_st)/se2, std(sehat2_st)/se2)
    end
    P4 = plot()
    stephist!(sehat2_we/se2; norm=true, label="We", c=1)
    vline!([mean(sehat2_we)/se2]; label="", ls=:dot, c=1)
    stephist!(sehat2_st/se2; norm=true, label="St", ls=:dash, c=2)
    vline!([mean(sehat2_st)/se2]; label="", ls=:dot, c=2)
    plot!(; xlim = quantile.(([sehat2_we/se2; sehat2_st/se2],), (0.0, 0.99)))
    vline!([1.0]; label="", ls=:dot, c=:gray)
    title!("SEhat²/SE²"; titlefontsize=10)
    
    P5 = plot()
    stephist!(tval_bm; norm=true, label="BM", ls=:dashdot, c=3)
    stephist!(tval_mw; norm=true, label="WMW", ls=:dashdotdot, c=4)
    plot!(Normal(), -5, 5; label="", ls=:dot, c=:gray)
    vline!([0.0]; label="", ls=:dot, c=:gray)
    plot!(; yscale=:log10, ylim=(1e-4, 1))
    title!("t-values"; titlefontsize=10)

    se2 = var(phat)
    if verbose
        @show (mean(sehat2_bm)/se2, std(sehat2_bm)/se2)
        @show mean(sehat2_mw)/se2
    end
    P6 = plot()
    stephist!(sehat2_bm/se2; norm=true, label="BM", ls=:dashdot, c=3)
    vline!([mean(sehat2_bm)/se2]; label="", ls=:dot, c=3)
    vline!([mean(sehat2_mw)/se2]; label="WMW", ls=:dashdotdot, c=4)
    vline!([1.0]; label="", ls=:dot, c=:gray)
    title!("SEhat²/SE²"; titlefontsize=10)
    
    layout = @layout [[a{0.2h} ; b; c] [d; e; f; g]]
    plot(P1, P2a, P2b, P3, P4, P5, P6; size=(800, 1000), layout)
    plot!(leftmargin=4Plots.mm, bottommargin=4Plots.mm)
end

plot_sim(; L=10^5)

# %%
@printf "%f %f" π exp(1)

# %%
plot_sim(; shifttype=:mean)

# %%
plot_sim(; shifttype=:tie)

# %%
plot_sim(; shifttype=:median)

# %%
distx = inversegammadist(1, 5)
disty = inversegammadist(1, 5)
plot_sim(; distx, disty, m=25, n=100)

# %%
distx = inversegammadist(1, 5)
disty = inversegammadist(1, 5)
plot_sim(; distx, disty, m=100, n=400, L=10^5)

# %%
distx = inversegammadist(1, 5)
disty = inversegammadist(1, 5)
plot_sim(; distx, disty, m=400, n=1600, L=10^5)

# %%
distx = inversegammadist(1.2, 5)
disty = inversegammadist(1, 5)
plot_sim(; distx, disty, m=25, n=100)

# %%
distx = inversegammadist(1.2, 5)
disty = inversegammadist(1, 5)
plot_sim(; distx, disty, m=100, n=400, L=10^5)

# %%
distx = inversegammadist(1.2, 5)
disty = inversegammadist(1, 5)
plot_sim(; distx, disty, m=400, n=1600, L=10^5)

# %%
distx = inversegammadist(2, 5)
disty = inversegammadist(1, 5)
plot_sim(; distx, disty, m=25, n=100)

# %%
distx = inversegammadist(1, 5)
disty = inversegammadist(1.2, 5)
plot_sim(; distx, disty, m=25, n=100)

# %%
distx = inversegammadist(1, 5)
disty = inversegammadist(1.5, 5)
plot_sim(; distx, disty, m=25, n=100)

# %%
distx = inversegammadist(1, 5)
disty = inversegammadist(1.5, 5)
plot_sim(; distx, disty, m=25, n=100, shifttype=:none)

# %%
distx = gammadist(1, 1.5)
disty = gammadist(1, 1.5)
plot_sim(; distx, disty, m=25, n=100)

# %%
distx = gammadist(1.2, 1.5)
disty = gammadist(1, 1.5)
plot_sim(; distx, disty, m=25, n=100)

# %%
distx = gammadist(1.5, 1.5)
disty = gammadist(1, 1.5)
plot_sim(; distx, disty, m=25, n=100)

# %%
distx = gammadist(1, 1.5)
disty = gammadist(1.2, 1.5)
plot_sim(; distx, disty, m=25, n=100)

# %%
distx = gammadist(1, 1.5)
disty = gammadist(1.5, 1.5)
plot_sim(; distx, disty, m=25, n=100)

# %%
distx = gammadist(1, 1.5)
disty = gammadist(1.5, 1.5)
plot_sim(; distx, disty, m=25, n=100, shifttype=:none)

# %%
distx = inversegammadist(2, 5)
disty = inversegammadist(1, 5)
plot_sim(; distx, disty, m=100, n=100)

# %%
distx = inversegammadist(2, 5)
disty = inversegammadist(1, 5)
plot_sim(; distx, disty, m=90, n=100)

# %%
distx = inversegammadist(3, 5)
disty = inversegammadist(1, 5)
plot_sim(; distx, disty, m=90, n=100)

# %%
distx = inversegammadist(3, 5)
disty = inversegammadist(1, 5)
plot_sim(; distx, disty, m=360, n=400, L=10^5)

# %%
distx = inversegammadist(2, 5)
disty = inversegammadist(1, 5)
plot_sim(; distx, disty, m=80, n=100)

# %%
distx = inversegammadist(2, 5)
disty = inversegammadist(1, 5)
plot_sim(; distx, disty, m=320, n=400, L=10^5)

# %%
distx = inversegammadist(3, 5)
disty = inversegammadist(1, 5)
plot_sim(; distx, disty, m=80, n=100)

# %% tags=[]
distx = inversegammadist(3, 5)
disty = inversegammadist(1, 5)
plot_sim(; distx, disty, m=320, n=400, L=10^5)

# %%
distx = inversegammadist(1, 5)
disty = inversegammadist(3, 5)
plot_sim(; distx, disty, m=80, n=100)

# %% tags=[]
distx = inversegammadist(1, 5)
disty = inversegammadist(3, 5)
plot_sim(; distx, disty, m=320, n=400, L=10^5)

# %%
distx = gammadist(3, 1.5)
disty = gammadist(1, 1.5)
plot_sim(; distx, disty, m=80, n=100)

# %%
distx = gammadist(3, 1.5)
disty = gammadist(1, 1.5)
plot_sim(; distx, disty, m=320, n=400, L=10^5)

# %%
distx = gammadist(1, 1.5)
disty = gammadist(3, 1.5)
plot_sim(; distx, disty, m=80, n=100)

# %%
distx = gammadist(1, 1.5)
disty = gammadist(3, 1.5)
plot_sim(; distx, disty, m=320, n=400, L=10^5)

# %% [markdown]
# * https://onlinelibrary.wiley.com/doi/10.1002/sim.3561
#   * https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fsim.3561&file=sim_3561_sm_SupplMat.pdf

# %%
distx = gammadist(2, 3)
disty = gammadist(1, 3)
plot_sim(; distx, disty, m=25, n=25, ylim=(-0.1, 3.1), shifttype=:mean)

# %%
distx = gammadist(2, 3)
disty = gammadist(1, 3)
plot_sim(; distx, disty, m=25, n=25, ylim=(-0.1, 3.1), shifttype=:median)

# %%
distx = gammadist(2, 3)
disty = gammadist(1, 3)
plot_sim(; distx, disty, m=25, n=25, ylim=(-0.1, 3.1), shifttype=:tie)

# %%
distx = gammadist(2, 3)
disty = gammadist(1, 3)
plot_sim(; distx, disty, m=25, n=100, ylim=(-0.1, 3.1), shifttype=:mean)

# %%
distx = gammadist(2, 3)
disty = gammadist(1, 3)
plot_sim(; distx, disty, m=25, n=100, ylim=(-0.1, 3.1), shifttype=:median, verbose=true)

# %%
distx = gammadist(2, 3) - 0.2711292144816104
disty = gammadist(1, 3)
plot_sim(; distx, disty, m=25, n=100, ylim=(-0.1, 3.1), shifttype=:none, verbose=true)

# %%
distx = gammadist(2, 3) - 0.2711292144816104 - 0.042
disty = gammadist(1, 3)
plot_sim(; distx, disty, m=25, n=100, ylim=(-0.1, 3.1), shifttype=:none, verbose=true)

# %%
distx = gammadist(2, 3) - 0.2711292144816104 + 0.01
disty = gammadist(1, 3)
plot_sim(; distx, disty, m=100, n=25, ylim=(-0.1, 3.1), shifttype=:none, verbose=true)

# %%
distx = gammadist(2, 3)
disty = gammadist(1, 3)
plot_sim(; distx, disty, m=100, n=25, ylim=(-0.1, 3.1), shifttype=:mean, verbose=true)

# %%
