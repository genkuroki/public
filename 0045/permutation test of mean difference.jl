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

# %% tags=[]
using Printf
using Random
Random.seed!(4649373)

using Distributions
using QuadGK
using Roots
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

この函数は数値のベクトルのデータ `X`, `Y` について, 
Brunner-Munzel検定関係の統計量達を計算する. 詳細は以下の通り.

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
    pvalue = sehat > 0 ? 2ccdf(TDist(df), abs(tvalue)) : phat ≈ p ? 1.0 : 0.0
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
x ⪅ y = x < y || x ≈ y

function pvalue_permutation_test(f, X, Y, XY=[X; Y]; Nshuffles=10^4)
    m, n = length(X), length(Y)
    f_obs = f(X, Y)
    a = b = 0
    for _ in 1:Nshuffles
        shuffle!(XY)
        @views X, Y = XY[1:m], XY[m+1:m+n]
        a += f(X, Y) ⪅ f_obs
        b += f_obs ⪅ f(X, Y)
    end
    min(1, 2a/Nshuffles, 2b/Nshuffles)
end

function pvalue_permutation_diffmean_test(X, Y, XY=[X; Y]; Nshuffles=10^4)
    pvalue_permutation_test(X, Y, XY; Nshuffles) do X, Y
        mean(X) - mean(Y)
    end
end

function pvalue_permutation_student_t_test(X, Y, XY=[X; Y]; Nshuffles=10^4)
    pvalue_permutation_test(X, Y, XY; Nshuffles) do X, Y
        m, n = length(X), length(Y)
        S² = ((m-1)*var(X) + (n-1)*var(Y)) / (m+n-2)
        (mean(X) - mean(Y)) / √(S²*(1/m + 1/n))
    end
end

function pvalue_permutation_welch_t_test(X, Y, XY=[X; Y]; Nshuffles=10^4)
    pvalue_permutation_test(X, Y, XY; Nshuffles) do X, Y
        (mean(X) - mean(Y)) / √(var(X)/length(X) + var(Y)/length(Y))
    end
end

function true_pvalue(f, distx, disty, X, Y; Nsamples=10^4)
    f_obs = f(X, Y)
    Xtmp = similar(X)
    Ytmp = similar(Y)
    a = b = 0
    for _ in 1:Nsamples
        XX = rand!(distx, Xtmp)
        YY = rand!(disty, Ytmp)
        a += f(XX, YY) ⪅ f_obs
        b += f_obs ⪅ f(XX, YY)
    end
    min(1, 2a/Nsamples, 2b/Nsamples)
end

function true_pvalue_diffmean(distx, disty, X, Y; Nsamples=10^4)
    true_pvalue(distx, disty, X, Y; Nsamples) do X, Y
        mean(X) - mean(Y)
    end
end

function true_pvalue_student_t_test(distx, disty, X, Y; Nsamples=10^4)
    true_pvalue(distx, disty, X, Y; Nsamples) do X, Y
        m, n = length(X), length(Y)
        S² = ((m-1)*var(X) + (n-1)*var(Y)) / (m+n-2)
        (mean(X) - mean(Y)) / √(S²*(1/m + 1/n))
    end
end

function true_pvalue_welch_t_test(distx, disty, X, Y; Nsamples=10^4)
    true_pvalue(distx, disty, X, Y; Nsamples) do X, Y
        (mean(X) - mean(Y)) / √(var(X)/length(X) + var(Y)/length(Y))
    end
end

function pvalue_bootstrap(f, X, Y; Nsamples=10^4)
    Xtmp = similar(X)
    Ytmp = similar(Y)
    a = b = 0
    for _ in 1:Nsamples
        XX = sample!(X, Xtmp)
        YY = sample!(Y, Ytmp)
        a += f(XX, YY) ≤ 0
        b += 0 ≤ f(XX, YY)
    end
    min(1, 2a/Nsamples, 2b/Nsamples)
end

function pvalue_bootstrap_diffmean(X, Y; Nsamples=10^4)
    pvalue_bootstrap(X, Y; Nsamples) do XX, YY
        mean(XX) - mean(YY)
    end
end

# %%
X = [2, 4, 3, 4, 4, 5]
Y = [1, 7, 8, 9, 7, 8]
@show pvalue_permutation_diffmean_test(X, Y; Nshuffles=10^7)
@show pvalue_permutation_student_t_test(X, Y; Nshuffles=10^7)
@show pvalue_student_t_test(X, Y)
@show pvalue_permutation_welch_t_test(X, Y; Nshuffles=10^7)
@show pvalue_welch_t_test(X, Y)
@show pvalue_bootstrap_diffmean(X, Y; Nsamples=10^7)
@show pvalue_mann_whitney_u_test(X, Y)
@show pvalue_brunner_munzel_test(X, Y)
;

# %%
X = [2, 4, 3, 4, 4, 2, 3, 3, 4]
Y = [1, 7, 8, 9]
@show pvalue_permutation_diffmean_test(X, Y; Nshuffles=10^6)
@show pvalue_permutation_student_t_test(X, Y; Nshuffles=10^6)
@show pvalue_student_t_test(X, Y)
@show pvalue_permutation_welch_t_test(X, Y; Nshuffles=10^6)
@show pvalue_welch_t_test(X, Y)
@show pvalue_bootstrap_diffmean(X, Y; Nsamples=10^7)
@show pvalue_mann_whitney_u_test(X, Y)
@show pvalue_brunner_munzel_test(X, Y)
;

# %%
Random.seed!(4649373_4)
distx = Normal(0, 1)
disty = Normal(0, 2)
m, n = 20, 10
X = rand(distx, m)
Y = rand(disty, n)
@show mean(X) - mean(Y)
@show std(X), std(Y)
@show (mean(X) - mean(Y)) / √(var(X)/length(X) + var(Y)/length(Y))
println()

@show pvalue_permutation_diffmean_test(X, Y; Nshuffles=10^6)
@show pvalue_permutation_student_t_test(X, Y; Nshuffles=10^6)
@show pvalue_student_t_test(X, Y)
@show true_pvalue_diffmean(distx, disty, X, Y; Nsamples=10^6)
@show true_pvalue_student_t_test(distx, disty, X, Y; Nsamples=10^6)
@show pvalue_permutation_welch_t_test(X, Y; Nshuffles=10^6)
@show pvalue_welch_t_test(X, Y)
@show true_pvalue_welch_t_test(distx, disty, X, Y; Nsamples=10^6)
@show pvalue_bootstrap_diffmean(X, Y; Nsamples=10^7)
@show pvalue_mann_whitney_u_test(X, Y)
@show pvalue_brunner_munzel_test(X, Y)
;

# %%
Random.seed!(4649373_11)
distx = Normal(0, 1)
disty = Normal(0, 2)
m, n = 10, 20
X = rand(distx, m)
Y = rand(disty, n)
@show mean(X) - mean(Y)
@show std(X), std(Y)
@show (mean(X) - mean(Y)) / √(var(X)/length(X) + var(Y)/length(Y))
println()

@show pvalue_permutation_diffmean_test(X, Y; Nshuffles=10^6)
@show pvalue_permutation_student_t_test(X, Y; Nshuffles=10^6)
@show pvalue_student_t_test(X, Y)
@show true_pvalue_diffmean(distx, disty, X, Y; Nsamples=10^6)
@show true_pvalue_student_t_test(distx, disty, X, Y; Nsamples=10^6)
@show pvalue_permutation_welch_t_test(X, Y; Nshuffles=10^6)
@show pvalue_welch_t_test(X, Y)
@show true_pvalue_welch_t_test(distx, disty, X, Y; Nsamples=10^6)
@show pvalue_bootstrap_diffmean(X, Y; Nsamples=10^7)
@show pvalue_mann_whitney_u_test(X, Y)
@show pvalue_brunner_munzel_test(X, Y)
;

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

# %% tags=[]
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

# %%
distx = inversegammadist(4, 10)
disty = Normal(mean(distx), 1)
plot(1distx; label="distx")
plot!(1disty; label="disty", ls=:dash)
plot!(xlim=(-0.1, 20))
plot!(size=(300, 180)) |> display

@show mean(distx) ≈ mean(disty)
@show std(distx), std(disty)
println()
m, n = 30, 60

L = Nshuffles = 5000
@time plot_sim_diff_shuffle(; distx, m, disty, n, L, Nshuffles)

# %%
distx = inversegammadist(4, 10)
disty = Normal(mean(distx), 1)
plot(1distx; label="distx")
plot!(1disty; label="disty", ls=:dash)
plot!(xlim=(-0.1, 20))
plot!(size=(300, 180)) |> display

@show mean(distx) ≈ mean(disty)
@show std(distx), std(disty)
println()
m, n = 50, 100

L = Nshuffles = 5000
@time plot_sim_diff_shuffle(; distx, m, disty, n, L, Nshuffles)

# %%
distx = inversegammadist(4, 10)
disty = Normal(mean(distx), 1)
plot(1distx; label="distx")
plot!(1disty; label="disty", ls=:dash)
plot!(xlim=(-0.1, 20))
plot!(size=(300, 180)) |> display

@show mean(distx) ≈ mean(disty)
@show std(distx), std(disty)
println()
m, n = 100, 200

L = Nshuffles = 5000
@time plot_sim_diff_shuffle(; distx, m, disty, n, L, Nshuffles)

# %%
distx = inversegammadist(4, 10)
disty = Normal(mean(distx), 1)
plot(1distx; label="distx")
plot!(1disty; label="disty", ls=:dash)
plot!(xlim=(-0.1, 20))
plot!(size=(300, 180)) |> display

@show mean(distx) ≈ mean(disty)
@show std(distx), std(disty)
println()
m, n = 200, 400

L = Nshuffles = 5000
@time plot_sim_diff_shuffle(; distx, m, disty, n, L, Nshuffles)

# %%
distx = inversegammadist(4, 10)
disty = Normal(mean(distx), 1)
plot(1distx; label="distx")
plot!(1disty; label="disty", ls=:dash)
plot!(xlim=(-0.1, 20))
plot!(size=(300, 180)) |> display

@show mean(distx) ≈ mean(disty)
@show std(distx), std(disty)
println()
m, n = 20, 40

L = Nshuffles = 5000
@time plot_sim_diff_shuffle(; distx, m, disty, n, L, Nshuffles)

# %%
distx = inversegammadist(4, 10)
disty = Normal(mean(distx), 1)
plot(1distx; label="distx")
plot!(1disty; label="disty", ls=:dash)
plot!(xlim=(-0.1, 20))
plot!(size=(300, 180)) |> display

@show mean(distx) ≈ mean(disty)
@show std(distx), std(disty)
println()
m, n = 10, 20

L = Nshuffles = 5000
@time plot_sim_diff_shuffle(; distx, m, disty, n, L, Nshuffles)

# %%
distx = inversegammadist(4, 10)
disty = Normal(mean(distx), 1)
plot(1distx; label="distx")
plot!(1disty; label="disty", ls=:dash)
plot!(xlim=(-0.1, 20))
plot!(size=(300, 180)) |> display

@show mean(distx) ≈ mean(disty)
@show std(distx), std(disty)
println()
m, n = 5, 10

L = Nshuffles = 5000
@time plot_sim_diff_shuffle(; distx, m, disty, n, L, Nshuffles)

# %%
distx = Normal(0, 4)
disty = Normal(0, 1)
plot(1distx; label="distx")
plot!(1disty; label="disty", ls=:dash)
plot!(xlim=(-15, 15))
plot!(size=(300, 180)) |> display

@show mean(distx) ≈ mean(disty)
@show std(distx), std(disty)
println()
m, n = 10, 20

L = Nshuffles = 5000
@time plot_sim_diff_shuffle(; distx, m, disty, n, L, Nshuffles)

# %%
distx = Normal(0, 4)
disty = Normal(0, 1)
plot(1distx; label="distx")
plot!(1disty; label="disty", ls=:dash)
plot!(xlim=(-15, 15))
plot!(size=(300, 180)) |> display

@show mean(distx) ≈ mean(disty)
@show std(distx), std(disty)
println()
m, n = 5, 10

L = Nshuffles = 5000
@time plot_sim_diff_shuffle(; distx, m, disty, n, L, Nshuffles)

# %%
distx = inversegammadist(4, 2)
disty = Normal(mean(distx), 1)
plot(1distx; label="distx")
plot!(1disty; label="disty", ls=:dash)
plot!(xlim=(0, 20))
plot!(size=(300, 180)) |> display

@show mean(distx) ≈ mean(disty)
@show std(distx), std(disty)
println()
m, n = 50, 100

L = Nshuffles = 5000
@time plot_sim_diff_shuffle(; distx, m, disty, n, L, Nshuffles)

# %%
distx = inversegammadist(4, 2)
disty = Normal(mean(distx), 1)
plot(1distx; label="distx")
plot!(1disty; label="disty", ls=:dash)
plot!(xlim=(0, 20))
plot!(size=(300, 180)) |> display

@show mean(distx) ≈ mean(disty)
@show std(distx), std(disty)
println()
m, n = 20, 40

L = Nshuffles = 5000
@time plot_sim_diff_shuffle(; distx, m, disty, n, L, Nshuffles)

# %%
distx = inversegammadist(4, 2)
disty = Normal(mean(distx), 1)
plot(1distx; label="distx")
plot!(1disty; label="disty", ls=:dash)
plot!(xlim=(0, 20))
plot!(size=(300, 180)) |> display

@show mean(distx) ≈ mean(disty)
@show std(distx), std(disty)
println()
m, n = 10, 20

L = Nshuffles = 5000
@time plot_sim_diff_shuffle(; distx, m, disty, n, L, Nshuffles)

# %%
distx = inversegammadist(4, 2)
disty = Normal(mean(distx), 1)
plot(1distx; label="distx")
plot!(1disty; label="disty", ls=:dash)
plot!(xlim=(0, 20))
plot!(size=(300, 180)) |> display

@show mean(distx) ≈ mean(disty)
@show std(distx), std(disty)
println()
m, n = 5, 10

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
@show mean(distx), mean(disty)
@show std(distx), std(disty)
println()
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
distx, disty = Normal(0, 2), Normal(0, 1)
m, n = 25, 100
L = Nshuffles = 10000
@time plot_sim_diff_shuffle(; distx, m, disty, n, L, Nshuffles)

# %%
