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
#     display_name: Julia 1.9.0-beta3
#     language: julia
#     name: julia-1.9
# ---

# %%
using Distributions
using HypothesisTests
using Random
using StatsBase
using StatsPlots
default(fmt=:png, titlefontsize=10, tickfontsize=6, guidefontsize=8)

rd(x; d=3) = round(x; digits=d)

function simMWU(distX, distY, m, n; L = 10^6)
    tmpX = [Vector{Int}(undef, m) for _ in 1:Threads.nthreads()]
    tmpY = [Vector{Int}(undef, n) for _ in 1:Threads.nthreads()]
    pval = Vector{Float64}(undef, L)
    Threads.@threads for i in 1:L
        X = rand!(distX, tmpX[Threads.threadid()])
        Y = rand!(distY, tmpY[Threads.threadid()])
        pval[i] = pvalue(MannWhitneyUTest(X, Y))
    end
    pval
    _ecdf_pval = ecdf(pval)
    ecdf_pval(x) = _ecdf_pval(x)
    ecdf_pval
end

# %%
distX = DiscreteUniform(1, 5)
distY = Categorical(0.05, 0.2, 0.5, 0.2, 0.05)
m, n = 100, 200
@show m, n
@show rd.((mean(distX), mean(distY)))
@show rd.((var(distX), var(distY)))
@show rd.((skewness(distX), skewness(distY)))
@show rd.((kurtosis(distX), kurtosis(distY)))
println()

ymax = maximum([probs(distX); probs(distY)]) + 0.02
P1 = bar(1:5, x -> pdf(distX, round(Int, x)); label="", title="distX, m=$m")
plot!(ylim=(-0.02, ymax))
P2 = bar(1:5, x -> pdf(distY, round(Int, x)); label="", title="distY, n=$n", c=2)
plot!(ylim=(-0.02, ymax))

ecdf_pval = simMWU(distX, distY, m, n)
@show ecdf_pval(0.05)
@show ecdf_pval(0.01)
P3 = plot(ecdf_pval, 0, 0.1; label="")
plot!([0, 0.1], [0, 0.1]; label="", c=:black, ls=:dot, lw=0.5)
plot!(xguide="nominal significance level α", yguide="probability of pvalue ≤ α")
plot!(xtick=0:0.01:1, ytick=0:0.01:1, xrotation=30)
title!("MannWhitneyUTest")
plot!(size=(450, 450))

plot(P1, P2, P3; size=(600, 400), layout=@layout[[a; b] c{0.7w}])

# %%
safediv(x, y) = x == 0 ? x : isinf(y) ? zero(y) : x/y

"""
    h_brunner_munzel(x, y)

この函数は, x < y のとき 1.0 を, x = y のとき 0.5 を返す.
"""
h_brunner_munzel(x, y) = (x < y) + (x == y)/2

@doc raw"""
    phat_brunner_munzel(X, Y)

まず以下のようにおく:

```math
\begin{aligned}
&
H(x, y) = \begin{cases} 1 & (x < y) \\ 1/2 & (x = y), \end{cases}
\\ &
m = \mathrm{length}(X), \quad
n = \mathrm{length}(Y), \quad
x_i = X[i], \quad
y_j = Y[j]
\end{aligned}
```

この函数は次の $\hat{p}$ を返す:

```math
\hat{p} = \frac{1}{mn}\sum_{i=1}^m \sum_{j=1}^n H(x_i, y_j).
```
"""
phat_brunner_munzel(X, Y) = mean(h_brunner_munzel(x, y) for x in X, y in Y)

@doc raw"""
    statistics_brunner_munzel(X, Y,
        Hx = similar(X, Float64),
        Hy = similar(Y, Float64);
        p = 1/2
    )

この函数はデータ `X`, `Y` について, Brunner-Munzel検定関係の統計量達を計算する. 詳細は以下の通り.

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
H(x, y) = \begin{cases} 1 & (x < y) \\ 1/2 & (x = y), \end{cases}
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
\hat{p} = \frac{\bar{H}^x - \bar{H}^y + n}{m + n},
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
2\mathrm{ccdf}(\mathrm{TDist}(\nu), |t|).
\end{aligned}
```
"""
function statistics_brunner_munzel(X, Y,
        Hx = similar(X, Float64),
        Hy = similar(Y, Float64);
        p = 1/2
    )
    m, n = length(X), length(Y)
    for (i, x) in pairs(X)
        Hx[i] = sum(h_brunner_munzel(y, x) for y in Y)
    end
    for (j, y) in pairs(Y)
        Hy[j] = sum(h_brunner_munzel(x, y) for x in X)
    end
    phat = (mean(Hy) - mean(Hx) + n)/(m + n)
    sx2, sy2 = var(Hx)/n^2, var(Hy)/m^2
    sehat = √(sx2/m + sy2/n)
    tvalue = (phat - p)/sehat
    df = safediv((sx2/m + sy2/n)^2, (sx2/m)^2/(m-1) + (sy2/n)^2/(n-1))
    pvalue = (df != 0 && isfinite(df)) ? 2ccdf(TDist(df), abs(tvalue)) : zero(df)
    (; phat, sx2, sy2, sehat, tvalue, df, pvalue)
end

@doc raw"""
    pvalue_brunner_munzel(X, Y,
        Hx = similar(X, Float64),
        Hy = similar(Y, Float64);
        p = 1/2
    )

この函数はBrunner-Munzel検定のP値 `pvalue` を返す.
"""
function pvalue_brunner_munzel(X, Y,
        Hx = similar(X, Float64),
        Hy = similar(Y, Float64);
        p = 1/2
    )
    statistics_brunner_munzel(X, Y, Hx, Hy; p).pvalue
end

"""
    tieshift(X::AbstractVector, Y::AbstractVector; p = 1/2)

この函数は `phat_brunner_munzel(X, Y .+ a)` の値が `p` に等しくなる `a` を返す.
"""
function tieshift(X::AbstractVector, Y::AbstractVector; p = 1/2)
    shiftmin = minimum(X) - maximum(Y) - 0.1
    shiftmax = maximum(X) - minimum(Y) + 0.1
    find_zero(a -> phat_brunner_munzel(X, Y .+ a) - p, (shiftmin, shiftmax))
end

@doc raw"""
    brunner_munzel(X, Y,
        Hx = similar(X, Float64),
        Hy = similar(Y, Float64),
        Ytmp = similar(Y, Float64);
        p = 1/2,
        α = 0.05,
        maxsplit = 30
    )

この函数はBrunner-Munzel検定を実行する. 詳細は以下の通り.

この函数は `phat`, `sehat`, `tvalue`, `df`, `p`, `pvalue`, `α` および\
以下達の named tuple を返す.

```math
\begin{aligned}
&
\mathrm{confint\_p} = (\text{$p$ の信頼度 $1-\alpha$ の信頼区間}),
\\ &
\mathrm{confint\_shift} = (\text{2つの集団が互角になるようなシフトの信頼度 $1-\alpha$ の信頼区間}),
\\ &
\mathrm{pvalue\_shift} = ($\mathrm{confint\_shift}$ の計算で使われたP値函数),
\\ &
\mathrm{shifthat} = (\text{2つの集団が互角になるようなシフトの点推定値}).
\end{aligned}
```

さらに, $\mathrm{shiftmin}$, $\mathrm{shiftmax}$ はデータから推定されるシフトの下限と上限.

"""
function brunner_munzel(X, Y,
        Hx = similar(X, Float64),
        Hy = similar(Y, Float64),
        Ytmp = similar(Y, Float64);
        p = 1/2,
        α = 0.05,
        maxsplit = 30
    )
    (; phat, sehat, tvalue, df, pvalue) = statistics_brunner_munzel(X, Y, Hx, Hy; p)
    
    c = df == 0 ? Inf : quantile(TDist(df), 1 - α/2)
    confint_p = [max(0, phat - c*sehat), min(1, phat + c*sehat)]
    
    function pvalue_shift(a)
        @. Ytmp = Y + a
        pvalue_brunner_munzel(X, Ytmp, Hx, Hy; p)
    end
    shiftmin = minimum(X) - maximum(Y) - 0.1
    shiftmax = maximum(X) - minimum(Y) + 0.1
    shifthat = tieshift(X, Y; p)
    confint_shift = [
        find_zero(a -> pvalue_shift(a) - α, (shiftmin, shifthat))
        find_zero(a -> pvalue_shift(a) - α, (shifthat, shiftmax))
    ]
    
    (; phat, sehat, tvalue, df, p, pvalue, α, confint_p,
        confint_shift, pvalue_shift, shifthat, shiftmin, shiftmax)
end

function show_plot_brunner_munzel(X, Y,
        Hx = similar(X, Float64),
        Hy = similar(Y, Float64),
        Ytmp = similar(Y, Float64);
        p = 1/2,
        α = 0.05,
        showXY = false,
        kwargs...
    )
    showXY && (@show X Y)
    (; phat, sehat, tvalue, df, p, pvalue, α, confint_p, 
        confint_shift, pvalue_shift, shifthat, shiftmin, shiftmax) =
        brunner_munzel(X, Y, Hx, Hy, Ytmp; p, α)
    pprint((; phat, sehat, tvalue, df, p, pvalue, α, confint_p,
            confint_shift, shifthat))
    println()
    @show median(X) median(Y)
    plot(pvalue_shift, shiftmin, shiftmax; label="")
    vline!([tieshift(X, Y)]; label="", ls=:dash)    
    title!("P-value function of shift")
    plot!(ytick=0:0.05:1)
    plot!(; kwargs...)
end

# %%
function simBM(distX, distY, m, n; L = 10^6)
    tmpX = [Vector{Int}(undef, m) for _ in 1:Threads.nthreads()]
    tmpY = [Vector{Int}(undef, n) for _ in 1:Threads.nthreads()]
    pval = Vector{Float64}(undef, L)
    Threads.@threads for i in 1:L
        X = rand!(distX, tmpX[Threads.threadid()])
        Y = rand!(distY, tmpY[Threads.threadid()])
        pval[i] = pvalue_brunner_munzel(X, Y)
    end
    pval
    _ecdf_pval = ecdf(pval)
    ecdf_pval(x) = _ecdf_pval(x)
    ecdf_pval
end

# %% tags=[]
m, n = 100, 200
ecdf_pval = simBM(distX, distY, m, n)
@show ecdf_pval(0.05)
@show ecdf_pval(0.01)
P4 = plot(ecdf_pval, 0, 0.1; label="")
plot!([0, 0.1], [0, 0.1]; label="", c=:black, ls=:dot, lw=0.5)
plot!(xguide="nominal significance level α", yguide="probability of pvalue ≤ α")
plot!(xtick=0:0.01:1, ytick=0:0.01:1, xrotation=30)
title!("Brunner-Munzel test")# for distX m=$m vs. distY n=$n")
plot!(size=(450, 450))

plot(P1, P2, P3, P4; size=(640, 512), layout=@layout[[a b]; [c{0.65h} d{0.65h}]])

# %%
function simT(distX, distY, m, n; L = 10^6)
    tmpX = [Vector{Int}(undef, m) for _ in 1:Threads.nthreads()]
    tmpY = [Vector{Int}(undef, n) for _ in 1:Threads.nthreads()]
    pval = Vector{Float64}(undef, L)
    Threads.@threads for i in 1:L
        X = rand!(distX, tmpX[Threads.threadid()])
        Y = rand!(distY, tmpY[Threads.threadid()])
        pval[i] = pvalue(UnequalVarianceTTest(X, Y))
    end
    pval
    _ecdf_pval = ecdf(pval)
    ecdf_pval(x) = _ecdf_pval(x)
    ecdf_pval
end

# %%
m, n = 100, 200
ecdf_pval = simT(distX, distY, m, n)
@show ecdf_pval(0.05)
@show ecdf_pval(0.01)

P5 = plot(ecdf_pval, 0, 0.1; label="")
plot!([0, 0.1], [0, 0.1]; label="", c=:black, ls=:dot, lw=0.5)
plot!(xguide="nominal significance level α", yguide="probability of pvalue ≤ α")
plot!(xtick=0:0.01:1, ytick=0:0.01:1, xrotation=30)
title!("Welch t-test")# for distX m=$m vs. distY n=$n")
plot!(size=(450, 450))

plot(P1, P2, P3, P5; size=(640, 512), layout=@layout[[a b]; [c{0.65h} d{0.65h}]])

# %%
p = 1/12
distX = Categorical(p, 2p, 1-6p, 2p, p)
q = 3p/2
distY = Categorical(q, 0, 1-2q, 0, q)
m, n = 50, 200
@show m, n
@show rd.((mean(distX), mean(distY)))
@show rd.((var(distX), var(distY)))
@show rd.((skewness(distX), skewness(distY)))
@show rd.((kurtosis(distX), kurtosis(distY)))
println()

ymax = maximum([probs(distX); probs(distY)]) + 0.02
P1 = bar(1:5, x -> pdf(distX, round(Int, x)); label="", title="distX, m=$m")
plot!(ylim=(-0.02, ymax))
P2 = bar(1:5, x -> pdf(distY, round(Int, x)); label="", title="distY, n=$n", c=2)
plot!(ylim=(-0.02, ymax))

ecdf_pval = simMWU(distX, distY, m, n)
@show ecdf_pval(0.05)
@show ecdf_pval(0.01)
P3 = plot(ecdf_pval, 0, 0.1; label="")
plot!([0, 0.1], [0, 0.1]; label="", c=:black, ls=:dot, lw=0.5)
plot!(xguide="nominal significance level α", yguide="probability of pvalue ≤ α")
plot!(xtick=0:0.01:1, ytick=0:0.01:1, xrotation=30)
title!("MannWhitneyUTest")
plot!(size=(450, 450))

plot(P1, P2, P3; size=(600, 400), layout=@layout[[a; b] c{0.7w}])

# %%
ecdf_pval = simBM(distX, distY, m, n)
@show ecdf_pval(0.05)
@show ecdf_pval(0.01)
P4 = plot(ecdf_pval, 0, 0.1; label="")
plot!([0, 0.1], [0, 0.1]; label="", c=:black, ls=:dot, lw=0.5)
plot!(xguide="nominal significance level α", yguide="probability of pvalue ≤ α")
plot!(xtick=0:0.01:1, ytick=0:0.01:1, xrotation=30)
title!("Brunner-Munzel test")# for distX m=$m vs. distY n=$n")
plot!(size=(450, 450))

plot(P1, P2, P3, P4; size=(640, 512), layout=@layout[[a b]; [c{0.65h} d{0.65h}]])

# %%
ecdf_pval = simT(distX, distY, m, n)
@show ecdf_pval(0.05)
@show ecdf_pval(0.01)

P5 = plot(ecdf_pval, 0, 0.1; label="")
plot!([0, 0.1], [0, 0.1]; label="", c=:black, ls=:dot, lw=0.5)
plot!(xguide="nominal significance level α", yguide="probability of pvalue ≤ α")
plot!(xtick=0:0.01:1, ytick=0:0.01:1, xrotation=30)
title!("Welch t-test")# for distX m=$m vs. distY n=$n")
plot!(size=(450, 450))

plot(P1, P2, P3, P5; size=(640, 512), layout=@layout[[a b]; [c{0.65h} d{0.65h}]])

# %%
p = 1/12
distX = Categorical(p, 2p, 1-6p, 2p, p)
q = 3p/2
distY = Categorical(q, 0, 1-2q, 0, q)
m, n = 200, 50
@show m, n
@show rd.((mean(distX), mean(distY)))
@show rd.((var(distX), var(distY)))
@show rd.((skewness(distX), skewness(distY)))
@show rd.((kurtosis(distX), kurtosis(distY)))
println()

ymax = maximum([probs(distX); probs(distY)]) + 0.02
P1 = bar(1:5, x -> pdf(distX, round(Int, x)); label="", title="distX, m=$m")
plot!(ylim=(-0.02, ymax))
P2 = bar(1:5, x -> pdf(distY, round(Int, x)); label="", title="distY, n=$n", c=2)
plot!(ylim=(-0.02, ymax))

ecdf_pval = simMWU(distX, distY, m, n)
@show ecdf_pval(0.05)
@show ecdf_pval(0.01)
P3 = plot(ecdf_pval, 0, 0.1; label="")
plot!([0, 0.1], [0, 0.1]; label="", c=:black, ls=:dot, lw=0.5)
plot!(xguide="nominal significance level α", yguide="probability of pvalue ≤ α")
plot!(xtick=0:0.01:1, ytick=0:0.01:1, xrotation=30)
title!("MannWhitneyUTest")
plot!(size=(450, 450))

plot(P1, P2, P3; size=(600, 400), layout=@layout[[a; b] c{0.7w}])

# %%
p = 1/12
distX = Categorical(p, 2p, 1-6p, 2p, p)
q = 3p/2
distY = Categorical(q, 0, 1-2q, 0, q)
m, n = 200, 100
@show m, n
@show rd.((mean(distX), mean(distY)))
@show rd.((var(distX), var(distY)))
@show rd.((skewness(distX), skewness(distY)))
@show rd.((kurtosis(distX), kurtosis(distY)))
println()

ymax = maximum([probs(distX); probs(distY)]) + 0.02
P1 = bar(1:5, x -> pdf(distX, round(Int, x)); label="", title="distX, m=$m")
plot!(ylim=(-0.02, ymax))
P2 = bar(1:5, x -> pdf(distY, round(Int, x)); label="", title="distY, n=$n", c=2)
plot!(ylim=(-0.02, ymax))

ecdf_pval = simMWU(distX, distY, m, n)
@show ecdf_pval(0.05)
@show ecdf_pval(0.01)
P3 = plot(ecdf_pval, 0, 0.1; label="")
plot!([0, 0.1], [0, 0.1]; label="", c=:black, ls=:dot, lw=0.5)
plot!(xguide="nominal significance level α", yguide="probability of pvalue ≤ α")
plot!(xtick=0:0.01:1, ytick=0:0.01:1, xrotation=30)
title!("MannWhitneyUTest")
plot!(size=(450, 450))

plot(P1, P2, P3; size=(600, 400), layout=@layout[[a; b] c{0.7w}])

# %%
p = 1/12
distX = Categorical(p, 2p, 1-6p, 2p, p)
q = 3p/2
distY = Categorical(q, 0, 1-2q, 0, q)
m, n = 100, 100
@show m, n
@show rd.((mean(distX), mean(distY)))
@show rd.((var(distX), var(distY)))
@show rd.((skewness(distX), skewness(distY)))
@show rd.((kurtosis(distX), kurtosis(distY)))
println()

ymax = maximum([probs(distX); probs(distY)]) + 0.02
P1 = bar(1:5, x -> pdf(distX, round(Int, x)); label="", title="distX, m=$m")
plot!(ylim=(-0.02, ymax))
P2 = bar(1:5, x -> pdf(distY, round(Int, x)); label="", title="distY, n=$n", c=2)
plot!(ylim=(-0.02, ymax))

ecdf_pval = simMWU(distX, distY, m, n)
@show ecdf_pval(0.05)
@show ecdf_pval(0.01)
P3 = plot(ecdf_pval, 0, 0.1; label="")
plot!([0, 0.1], [0, 0.1]; label="", c=:black, ls=:dot, lw=0.5)
plot!(xguide="nominal significance level α", yguide="probability of pvalue ≤ α")
plot!(xtick=0:0.01:1, ytick=0:0.01:1, xrotation=30)
title!("MannWhitneyUTest")
plot!(size=(450, 450))

plot(P1, P2, P3; size=(600, 400), layout=@layout[[a; b] c{0.7w}])

# %%
p = 1/12
distX = Categorical(p, 2p, 1-6p, 2p, p)
q = 3p/2
distY = Categorical(q, 0, 1-2q, 0, q)
m, n = 100, 200
@show m, n
@show rd.((mean(distX), mean(distY)))
@show rd.((var(distX), var(distY)))
@show rd.((skewness(distX), skewness(distY)))
@show rd.((kurtosis(distX), kurtosis(distY)))
println()

ymax = maximum([probs(distX); probs(distY)]) + 0.02
P1 = bar(1:5, x -> pdf(distX, round(Int, x)); label="", title="distX, m=$m")
plot!(ylim=(-0.02, ymax))
P2 = bar(1:5, x -> pdf(distY, round(Int, x)); label="", title="distY, n=$n", c=2)
plot!(ylim=(-0.02, ymax))

ecdf_pval = simMWU(distX, distY, m, n)
@show ecdf_pval(0.05)
@show ecdf_pval(0.01)
P3 = plot(ecdf_pval, 0, 0.1; label="")
plot!([0, 0.1], [0, 0.1]; label="", c=:black, ls=:dot, lw=0.5)
plot!(xguide="nominal significance level α", yguide="probability of pvalue ≤ α")
plot!(xtick=0:0.01:1, ytick=0:0.01:1, xrotation=30)
title!("MannWhitneyUTest")
plot!(size=(450, 450))

plot(P1, P2, P3; size=(600, 400), layout=@layout[[a; b] c{0.7w}])

# %%
