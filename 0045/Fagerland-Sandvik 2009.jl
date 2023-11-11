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

# %% [markdown]
# # Fagerland-Sandvik (2009) の再現
#
# * 黒木 玄
# * 2023-11-10

# %% [markdown] toc=true
# <h1>目次<span class="tocSkip"></span></h1>
# <div class="toc"><ul class="toc-item"><li><span><a href="#再現用の函数達の定義" data-toc-modified-id="再現用の函数達の定義-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>再現用の函数達の定義</a></span></li><li><span><a href="#再現" data-toc-modified-id="再現-2"><span class="toc-item-num">2&nbsp;&nbsp;</span>再現</a></span><ul class="toc-item"><li><span><a href="#訂正の必要性" data-toc-modified-id="訂正の必要性-2.1"><span class="toc-item-num">2.1&nbsp;&nbsp;</span>訂正の必要性</a></span></li><li><span><a href="#訂正版1" data-toc-modified-id="訂正版1-2.2"><span class="toc-item-num">2.2&nbsp;&nbsp;</span>訂正版1</a></span></li><li><span><a href="#訂正版2-(ガンマ分布を逆ガンマ分布で置き換えた場合)" data-toc-modified-id="訂正版2-(ガンマ分布を逆ガンマ分布で置き換えた場合)-2.3"><span class="toc-item-num">2.3&nbsp;&nbsp;</span>訂正版2 (ガンマ分布を逆ガンマ分布で置き換えた場合)</a></span></li><li><span><a href="#原論文の方法をそのまま再現" data-toc-modified-id="原論文の方法をそのまま再現-2.4"><span class="toc-item-num">2.4&nbsp;&nbsp;</span>原論文の方法をそのまま再現</a></span></li></ul></li></ul></div>

# %% [markdown]
# 論文 Fagerland-Sandvik (2009) https://onlinelibrary.wiley.com/doi/10.1002/sim.3561 の再現をやってみた. 結果は
#
# * https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fsim.3561&file=sim_3561_sm_SupplMat.pdf
#
# で公開されている. 以下はその論文の説明である.

# %% [markdown]
# ![Fagerland-Sandvik%202009%20A.png](attachment:Fagerland-Sandvik%202009%20A.png)

# %% [markdown]
# ![Fagerland-Sandvik%202009%20B.png](attachment:Fagerland-Sandvik%202009%20B.png)

# %% [markdown]
# ![Fagerland-Sandvik%202009%20C.png](attachment:Fagerland-Sandvik%202009%20C.png)

# %% [markdown]
# ## 再現用の函数達の定義

# %%
using Random
Random.seed!(4649373)

using Distributions
using QuadGK
using Roots
using StatsPlots
default(fmt=:png, titlefontsize=10, size=(400, 250))

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
X = fill(1, 100)
Y = fill(0, 100)
@show pvalue_brunner_munzel_test(X, Y), rcopy(lawstat.brunner_munzel_test(X, Y))[:p_value]
brunner_munzel_test(X, Y) |> pairs

# %%
X = fill(0, 100)
Y = fill(1, 100)
@show pvalue_brunner_munzel_test(X, Y), rcopy(lawstat.brunner_munzel_test(X, Y))[:p_value]
brunner_munzel_test(X, Y) |> pairs

# %%
X = zeros(10^3)
Y = zeros(10^3)
@show pvalue_brunner_munzel_test(X, Y), rcopy(lawstat.brunner_munzel_test(X, Y))[:p_value]
brunner_munzel_test(X, Y) |> pairs

# %%
X = zeros(10^3)
Y = zeros(10^3)
Y[1] += eps()
@show pvalue_brunner_munzel_test(X, Y), rcopy(lawstat.brunner_munzel_test(X, Y))[:p_value]
brunner_munzel_test(X, Y) |> pairs

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
gammadist(σ, β) = β == 0 ? Normal(0, σ) : Gamma(4/β^2, β*σ/2)

gam = gammadist.(1:5, 2:2:10)
[std.(gam), skewness.(gam), kurtosis.(gam), shape.(gam)]

# %%
β = range(eps(), 6, 1000)
P = plot(β, β -> kurtosis(gammadist(1, β)); label="")
plot!(xguide="skewness β", yguide="", title="kurtosis of gammadist(1, β)")

β = range(0.1, 6, 1000)
Q = plot(β, β -> shape(gammadist(1, β)); label="")
plot!(xguide="skewness β", yguide="", title="shape of gammadist(1, β)")

plot(P, Q; size=(800, 250))
plot!(bottommargin=4Plots.mm)

# %%
@doc gammadist

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
c = find_zero(β -> 3 + 8/β^2 + 4/β * √(1 + 4/β^2) - 4.0, 3.0)

# %%
β = range(eps(), 7, 1000)
P = plot(β, β -> kurtosis(inversegammadist(1, β)); label="")
vline!([c]; label="$(round(c; digits=6))", ls=:dash)
plot!(xtick=0:10, ylim=(-400, 10000))
plot!(xguide="skewness β", yguide="", title="kurtosis of inversegammadist(1, β)")

β = range(0.5, 10, 1000)
Q = plot(β, β -> shape(inversegammadist(1, β)); label="")
hline!([3]; label="3", ls=:dash)
plot!(xguide="skewness β", yguide="", title="shape of inversegammadist(1, β)")

plot(P, Q; size=(800, 250))
plot!(bottommargin=4Plots.mm)

# %%
@doc inversegammadist

# %%
undefarray(T::Type, n...) = Array{T}(undef, n...)
undefarray(n...) = undefarray(Float64, n...)
ECDF(A, x) = count(≤(x), A) / length(A)

list_skewness = Tuple((β, β) for β in 0:0.5:3)
list_skewness2 = Tuple((β, β-0.5) for β in 1:0.5:3)
list_stdratio = reverse((1.0, 1.1, 1.25, 1.5, 2.0))

rd(x) = round(100x; digits=1)

function print_sim(;
        list_skewness=list_skewness,
        list_stdratio = list_stdratio,
        distfunc = gammadist,
        distfuncx = distfunc,
        distfuncy = distfunc,
        α = 0.05,
        L = 10^5,
        m = 25,
        n = 25,
        shifttype = :auto,
        correct = true,
    )
    
    nx = length(list_skewness)
    ny = length(list_stdratio)
    # `er` stands for "true alpha Error Rate".
    er_wmw = undefarray(nx, ny)
    er_bm = undefarray(nx, ny)
    er_st = undefarray(nx, ny)
    er_we = undefarray(nx, ny)
    
    for (i, (β1, β2)) in enumerate(list_skewness), (j, σ) in enumerate(list_stdratio)
        distx = distfuncx(σ, β1)
        disty = distfuncy(1, β2)
        (; pval_wmw, pval_bm, pval_st, pval_we) = sim(distx, disty, m, n; shifttype, correct, L)
        er_wmw[i, j] = ECDF(pval_wmw, α)
        er_bm[i, j] = ECDF(pval_bm, α)
        er_st[i, j] = ECDF(pval_st, α)
        er_we[i, j] = ECDF(pval_we, α)
    end
    
    println("skewness = ", list_skewness)
    println("stdratio = ", list_stdratio)
    println("distx = $distfuncx,  disty = $distfuncy, m = $m,  n = $n,  shifttype = $shifttype")
    println()
    println("Wilcoxon-Mann-Whitney:"); Base.print_matrix(stdout, rd.(er_wmw'))
    println("\n\nBrunner-Munzel:"); Base.print_matrix(stdout, rd.(er_bm'))
    println("\n\nStudent t:"); Base.print_matrix(stdout, rd.(er_st'))
    println("\n\nWelch t:"); Base.print_matrix(stdout, rd.(er_we'))
    println("\n")
end

function sim(distx, disty, m, n; shifttype=:mean, correct=true, L=10^5)
    meansh = mean(disty) - mean(distx)
    mediansh = median(disty) - median(distx)
    tiesh = tieshift(distx, disty)
    if shifttype == :mean
        distx_sh = distx + meansh
    elseif shifttype == :median
        distx_sh = distx + mediansh
    elseif shifttype == :tie || shifttype == :auto
        distx_sh = distx + tiesh
    else
        distx_sh = distx
    end
    
    pval_wmw = undefarray(L)
    pval_bm = undefarray(L)
    pval_st = undefarray(L)
    pval_we = undefarray(L)
    nth = Threads.nthreads()
    Xtmp = [undefarray(m) for _ in 1:nth]
    Ytmp = [undefarray(n) for _ in 1:nth]
    Threads.@threads for i in 1:L
        tid = Threads.threadid()
        X = rand!(distx_sh, Xtmp[tid])
        Y = rand!(disty, Ytmp[tid])
        pval_wmw[i] = pvalue_mann_whitney_u_test(X, Y; correct)
        pval_bm[i] = pvalue_brunner_munzel_test(X, Y)
        if shifttype == :auto
            @. X = X - tiesh + meansh
        end
        pval_st[i] = pvalue_student_t_test(X, Y)
        pval_we[i] = pvalue_welch_t_test(X, Y)
    end
    
    (; pval_wmw, pval_bm, pval_st, pval_we)
end

function print_sim(distx, disty, m, n; shifttype=:mean, L=10^5, α=0.05)
    (; pval_wmw, pval_bm, pval_st, pval_we) = sim(distx, disty, m, n; shifttype, L)
    er_wmw = ECDF(pval_wmw, α)
    er_bm = ECDF(pval_bm, α)
    er_st = ECDF(pval_st, α)
    er_we = ECDF(pval_we, α)
    println("distx = $distx")
    println("disty = $disty")
    println("m = $m,  n = $n,  shifttype = $shifttype")
    println()
    println("Wilcoxon-Mann-Whitney: ", rd(er_wmw), "%")
    println("Brunner-Munzel:        ", rd(er_bm), "%")
    println("Student t:             ", rd(er_st), "%")
    println("Welch t:               ", rd(er_we), "%")
    println()
end

# %%
@time print_sim(gammadist(1.1, 3), gammadist(1, 3), 1000, 1000; shifttype=:mean, L=10^4)

# %% tags=[]
@time print_sim(gammadist(1.1, 3), gammadist(1, 3), 1000, 1000; shifttype=:median, L=10^4)

# %%
@time print_sim(gammadist(1.1, 3), gammadist(1, 3), 1000, 1000; shifttype=:tie, L=10^4)

# %% [markdown]
# ## 再現

# %% [markdown]
# ### 訂正の必要性
#
# 論文 Fagerland-Sandvik (2009) は不適切な比較をしていたので, その訂正版を以下で作成した.
#
# __解説:__ どこが不適切なのか.  論文 Fagerland-Sandvik (2009) では, 2つのテスト用の仮想的な母集団分布を適当にシフトして
#
# * 母平均を等しくした場合
#
# と
#
# * 母中央値を等しくした場合
#
# を扱っている. 
#
# 「2つの母平均は等しい」はStudentの $t$ 検定とWelchの $t$ 検定(Welchの $t$ 検定は Fagerland-Sandvik (2009) でがWelch U test (modified T test)と呼ばれている)の帰無仮説なので, それらのアルファエラー率を調べるときには「母平均を等しくした場合」を見る必要がある.
#
# しかし, 「2つの母中央値は等しい」はWilcoxon-Mann-WhitneyのU検定やBrunner-Munzel検定の帰無仮説ではない.  だから, それらのアルファエラー率を調べるときには「母中央値を等しくした場合」を見てはいけない.  それらの検定が中央値に関する検定でないことは論文でも指摘されているのだが, なぜかアルファエラー率を調べるために「2つの母中央値を等しくした場合」のみを計算してしまっている.
#
# WMW検定の検定統計量Uの定義は, 無作為抽出された2つの標本 $[X_1,\ldots,X_m]$, $[Y_1,\ldots,Y_n]$ について
#
# $$
# U = \sum_{i=1}^m \sum_{j=1}^n h(X_i, Y_j), \quad
# h(x, y) = \begin{cases}
# 1 & (x < y) \\
# 1/2 & (x = y) \\
# 0 & (x > y)
# \end{cases}
# $$
#
# であり, WMW検定では $U$ と $mn/n$ の差の絶対値が十分に大きいときに帰無仮説が棄却される.
#
# BM検定の検定統計量 $\hat{p}$ は
#
# $$
# \hat{p} = \frac{U}{mn} = \frac{1}{mn}\sum_{i=1}^m \sum_{j=1}^n h(X_i, Y_j)
# $$
#
# であり, BM検定では $\hat{p}$ と $1/2$ の差の絶対値が十分に大きいときに帰無仮説が棄却される
#
# WMW検定とBM検定の検定統計量は本質的に同じであると考えられる. 
#
# WMW検定とBM検定の本質的な違いは検定統計量の分散の見積もり方の違いである.
#
# 検定統計量 $\hat{p}$ は $Y$ 側の母勝率
#
# $$
# p = P(X_1 < Y_1) + P(X_1 = Y_1)/2
# $$
#
# の不偏一致推定量になっている. 
#
# そして, 検定統計量 $\hat{p}$ の分散の推定値は, WMW検定でもBM検定でも, $m,n\to\infty$ で $0$ に近付く.
#
# だから, もしも母勝率 $p$ が $1/2$ に等しくないならば, WMW検定でもBM検定でも, $m,n\to\infty$ で帰無仮説が棄却される確率は $1$ に近付く.
#
# 母中央値が等しいことと母勝率が $1/2$ に等しいことは同値ではない(これも Fagerland-Sandvik (2009) で指摘されている).
#
# だから, 母中央値が等しいという条件の下で, WMW検定やBM検定達のアルファエラー率を計算すると, 標本サイズを大きくするほど名目有意水準の $5 \%$ からかけ離れた値になって行くことになる.
#
# Fagerland-Sandvik (2009)は実際にそういう結果を得ている.
#
# そういう結果が得られることは当然なのに Fagerland-Sandvik (2009) では以下のように述べてしまっている. これはひどくミスリーディングである!

# %% [markdown]
# ![Fagerland-Sandvik%202009%20D.png](attachment:Fagerland-Sandvik%202009%20D.png)

# %% [markdown]
# __翻訳:__
#
# >興味深い結果は, 総標本サイズが大きくなると, 名目有意水準を維持するWMW検定の能力が低下することであった.  この効果は極めて大きく, WMW検定が2標本T検定(Studentのt検定)の漸近的頑健性を共有しないことを示す.
#
# >この結果をこのセクションの最初に定義した頑健性の基準にかけると, シミュレートされた有意水準の36%が10%頑強性を持ち, 11%が20%頑強性を持ち, 53%は頑強ではない.  驚くべきことに, 標準偏差の差が10％しかないことと歪度の度合いが中程度以上であることが合わさって, ほとんどの設定で頑強でない有意水準になった.  上述したように, この非頑強性は標本サイズが大きい場合に特に深刻である.  m = n = 1000の標本サイズでの追加シミュレーション(資料に載っていない)では, 母平均が等しいという帰無仮説のもとでは99%, 母中央値が等しいという帰無仮説のもとでは40%の実質有意水準が得られた($\theta = 1.1$ (標準偏差の比が $1.1$), $\beta_A = \beta_B = 3.0$ (母歪度はどちらも $3.0$))。
#
# >もしも目的が中央値の比較であれば, たとえ標準偏差が等しい場合であっても, WMW検定の有意水準は許容できないだろう.  これは、歪度の度合いが大きいと真の有意水準が高くなることはテーブルⅢで例証されています。

# %% [markdown]
# ![Fagerland-Sandvik%202009%20E.png](attachment:Fagerland-Sandvik%202009%20E.png)

# %% [markdown]
# WMW検定の検定統計量を見れば, WMW検定を母平均の差や母中央値の差の検定として使えないことは自明なのに, わざわざこのように述べるのはおかしい.
#
# この点の強調はWMW検定やBM検定について著しく不利な印象を読者に与えることになる.
#
# フェアな比較をするためには
#
# * 母勝率を $1/2$ にした場合
#
# でWMW検定やBM検定のアルファエラー率を計算してみる必要がある.

# %% [markdown]
# __テーブルⅢの再現__

# %%
@time print_sim(; m=25, n=100, shifttype=:mean)

# %% [markdown]
# 上の結果は論文 Fagerland-Sandvik (2009) のテーブルⅢの等平均の場合を再現できている.

# %%
@time print_sim(; m=25, n=100, shifttype=:median)

# %% [markdown]
# 上の結果は論文 Fagerland-Sandvik (2009) のテーブルⅢの党中央値の場合を再現できて__いない.__
#
# それは当然である. 標準偏差と歪度が等しい2つのガンマ分布は一致し, 2つの同じ分布の標本にWilcoxon-Mann-Whitney検定を適用すると帰無仮説が棄却される確率は名目有意水準に一致する

# %%
@time print_sim(; m=25, n=100, shifttype=:median, list_skewness=list_skewness2)

# %% [markdown]
# 上の結果を見れば, もしも目的が中央値の比較であれば, たとえ標準偏差が等しい場合であっても, WMW検定やBM検定を使うべきではないことがわかる.
#
# しかし, これはすでに述べたように自明であり, この点を強調することはミスリーディングである.
#
# WMW検定やBM検定は, それらの検定統計量を見れば, 「母勝率が $1/2$ である」という帰無仮説の検定だとみなされるべきであることがわかる.

# %% [markdown]
# __フェアな比較__

# %%
@time print_sim(; m=25, n=100)

# %% [markdown]
# この結果を見ると,
#
# * 「母勝率は $1/2$ である」という帰無仮説の検定ではWMW検定よりもBM検定を使うべきである.
# * 「2つの母平均は等しい」という帰無仮説の検定ではStudentのt検定よりもWelchのt検定を使うべきである.
#
# ということがわかる.  ただし,
#
# * 例外的に「2つの母分散は等しい」という条件を2つの母集団分布が満たしていることが事前に何らかの理由で分かっている場合には例外的に, 「2つの母平均は等しい」という帰無仮説の検定ではWelchのt検定よりもStudentのt検定を使った方がよい.
#
# これは非常に例外的な場合である.  Welchのt検定の非頑強性は標本サイズが十分に大きくないことによって生じており, 下の方のセルでの計算のように標本サイズを大きくすると緩和される.
#
# Welchのt検定の頑強性の基礎は中心極限定理による標本平均の分布の正規分布近似である. 母歪度の絶対値が大きいと標本平均の分布の正規分布近似の精度を上げるためにはより大きな標本サイズが必要になる. 
#
# 中心極限定理によって標本平均の分布が十分に正規分布近似されていると期待される場合には, Studentの $t$ 検定のような脆弱な検定法は使わずに, Welchの $t$ 検定を使った方がよいだろう.

# %%
@time print_sim(; m=100, n=400, L=10^4)

# %%
distx = gammadist(2, 3)
disty = gammadist(1, 3)
@show tieshift(distx, disty)
@show mean(disty) - mean(distx)
@show median(disty) - median(distx);

# %% [markdown]
# ついでに述べておくと, 以下に引用する論文 Fagerland-Sandvik (2009) におけるWMW検定の基本に関する説明も間違っている.

# %% [markdown]
# ![Fagerland-Sandvik%202009%20F.png](attachment:Fagerland-Sandvik%202009%20F.png)

# %% [markdown]
# $U$ 統計量の分散が $mn(m+n+1)/12$ になることを示すためには, $P(X<Y)=1/2$ を仮定するだけでは足りず, $X$ と $Y$ の分布が互いに等しいという極めて強い条件を仮定する必要がある.  このことがWMW検定の脆弱性の主原因になっている.
#
# Brunner-Munzel検定では「$X$ と $Y$ の分布が互いに等しい」という仮定を大幅に緩めてある種の中心極限定理による正規分布近似だけを仮定して, 検定統計量 $\hat{p}=mnU$ の分散を推定している.

# %% [markdown]
# ### 訂正版1
#
# 以下では, WMW検定とBM検定については母集団分布の適当なシフトによって「母勝率を $1/2$ にした場合」でアルファエラー率を計算し, Studentのt検定とWelchのt検定については母集団分布の適当なシフトによって「母平均が等しい場合」でアルファエラー率を計算している(シフトタイプが自動(auto)の場合).

# %%
@time print_sim(; m=25, n=25, list_skewness=list_skewness)

# %%
@time print_sim(; m=50, n=50, list_skewness=list_skewness)

# %%
@time print_sim(; m=25, n=100, list_skewness=list_skewness)

# %%
@time print_sim(; m=100, n=25, list_skewness=list_skewness)

# %%
@time print_sim(; m=100, n=100, list_skewness=list_skewness)

# %%
@time print_sim(; m=25, n=25, list_skewness=list_skewness2)

# %%
@time print_sim(; m=50, n=50, list_skewness=list_skewness2)

# %%
@time print_sim(; m=25, n=100, list_skewness=list_skewness2)

# %%
@time print_sim(; m=100, n=25, list_skewness=list_skewness2)

# %%
@time print_sim(; m=100, n=100, list_skewness=list_skewness2)

# %% [markdown]
# 上で使われた最も歪度が大きな場合(歪度は $3$)のガンマ分布をプロットしてみよう. 

# %%
gammadist(1, 3)

# %%
@show skewness(gammadist(1, 3))
@show kurtosis(gammadist(1, 3));

# %%
plot(gammadist(1, 3), -0.1, 2.1; label="gammadist(1, 3)", ylim=(-0.1, 4.1))

# %% [markdown]
# このように, 歪度が $3$ のガンマ分布は左端で密度函数の値が無限大になるような分布になっている.  これはあまりにも特殊な状況であるようにも思われる.  そこで, 次の節ではガンマ分布を逆ガンマ分布で置き換えた場合も扱ってみよう.

# %%
inversegammadist(1, 3)

# %%
@show skewness(inversegammadist(1, 3))
@show kurtosis(inversegammadist(1, 3));

# %%
plot(inversegammadist(1, 3), -0.2, 6; label="inversegammadist(1, 3)")

# %% [markdown]
# ### 訂正版2 (ガンマ分布を逆ガンマ分布で置き換えた場合)

# %%
@time print_sim(; m=25, n=25, list_skewness=list_skewness, distfunc=inversegammadist)

# %%
@time print_sim(; m=50, n=50, list_skewness=list_skewness, distfunc=inversegammadist)

# %%
@time print_sim(; m=25, n=100, list_skewness=list_skewness, distfunc=inversegammadist)

# %%
@time print_sim(; m=100, n=25, list_skewness=list_skewness, distfunc=inversegammadist)

# %%
@time print_sim(; m=100, n=100, list_skewness=list_skewness, distfunc=inversegammadist)

# %%
@time print_sim(; m=25, n=25, list_skewness=list_skewness2, distfunc=inversegammadist)

# %%
@time print_sim(; m=50, n=50, list_skewness=list_skewness2, distfunc=inversegammadist)

# %%
@time print_sim(; m=25, n=100, list_skewness=list_skewness2, distfunc=inversegammadist)

# %%
@time print_sim(; m=100, n=25, list_skewness=list_skewness2, distfunc=inversegammadist)

# %%
@time print_sim(; m=100, n=100, list_skewness=list_skewness2, distfunc=inversegammadist)

# %% [markdown]
# ### 原論文の方法をそのまま再現

# %% [markdown]
# 以下の shifttype=:mean と shifttype=:median の場合が原論文の再現になっている.

# %%
@time print_sim(; m=25, n=25, shifttype=:mean, list_skewness=list_skewness)

# %%
@time print_sim(; m=25, n=25, shifttype=:median, list_skewness=list_skewness)

# %%
@time print_sim(; m=25, n=25, shifttype=:tie, list_skewness=list_skewness)

# %%
@time print_sim(; m=50, n=50, shifttype=:mean, list_skewness=list_skewness)

# %%
@time print_sim(; m=50, n=50, shifttype=:median, list_skewness=list_skewness)

# %%
@time print_sim(; m=50, n=50, shifttype=:tie, list_skewness=list_skewness)

# %%
@time print_sim(; m=25, n=100, shifttype=:mean, list_skewness=list_skewness)

# %%
@time print_sim(; m=25, n=100, shifttype=:median, list_skewness=list_skewness)

# %% [markdown]
# 上の結果は論文 Fagerland-Sandvik (2009) の Table 6 を再現しない.

# %% [markdown]
# ![Fagerland-Sandvik%202009%20Table%206.png](attachment:Fagerland-Sandvik%202009%20Table%206.png)

# %% [markdown]
# 論文 Fagerland-Sandvik (2009) の Table 6 (上に引用)の結果は正しくないように思われる. 
#
# このTable 6のWMWとBMの表におけるstd ratioが1.00の段の右端の値(それぞれ8.7と7.1)がおかしい.
#
# 歪度と標準偏差が等しい2つのガンマ分布は等しい.  ゆえにその場合には最初から中央値も等しい. その場合にはWMW検定の極めて強い帰無仮説の「2つの分布は等しい」が成立している場合である.  そのような場合に結果が歪度に依存するはずがない.  それにも関わらず, WMWの右下の値が8.7%と5%よりも結構大きな値になっているのはおかしい.
#
# 他のテーブルの数値も確認してみると, Table 8 にも同様の問題があり, 論文 Fagerland-Sandvik (2009) の equal median の場合の計算結果は信用できないように見える.
#
# 2つの分布の中央値を揃える計算の精度が低過ぎたのではないかと思われる.
#
# 以上の点について私が誤解しているならば詳しく教えて欲しい.

# %%
@time print_sim(; m=25, n=100, shifttype=:tie, list_skewness=list_skewness)

# %%
@time print_sim(; m=100, n=25, shifttype=:mean, list_skewness=list_skewness)

# %%
@time print_sim(; m=100, n=25, shifttype=:median, list_skewness=list_skewness)

# %%
@time print_sim(; m=100, n=25, shifttype=:tie, list_skewness=list_skewness)

# %%
@time print_sim(; m=100, n=100, shifttype=:mean, list_skewness=list_skewness)

# %%
@time print_sim(; m=100, n=100, shifttype=:median, list_skewness=list_skewness)

# %%
@time print_sim(; m=100, n=100, shifttype=:tie, list_skewness=list_skewness)

# %%
@time print_sim(; m=25, n=25, shifttype=:mean, list_skewness=list_skewness2)

# %%
@time print_sim(; m=25, n=25, shifttype=:median, list_skewness=list_skewness2)

# %%
@time print_sim(; m=25, n=25, shifttype=:tie, list_skewness=list_skewness2)

# %%
@time print_sim(; m=50, n=50, shifttype=:mean, list_skewness=list_skewness2)

# %%
@time print_sim(; m=50, n=50, shifttype=:median, list_skewness=list_skewness2)

# %%
@time print_sim(; m=50, n=50, shifttype=:tie, list_skewness=list_skewness2)

# %%
@time print_sim(; m=25, n=100, shifttype=:mean, list_skewness=list_skewness2)

# %%
@time print_sim(; m=25, n=100, shifttype=:median, list_skewness=list_skewness2)

# %%
@time print_sim(; m=25, n=100, shifttype=:tie, list_skewness=list_skewness2)

# %%
@time print_sim(; m=100, n=25, shifttype=:mean, list_skewness=list_skewness2)

# %%
@time print_sim(; m=100, n=25, shifttype=:median, list_skewness=list_skewness2)

# %%
@time print_sim(; m=100, n=25, shifttype=:tie, list_skewness=list_skewness2)

# %%
@time print_sim(; m=100, n=100, shifttype=:mean, list_skewness=list_skewness2)

# %%
@time print_sim(; m=100, n=100, shifttype=:median, list_skewness=list_skewness2)

# %%
@time print_sim(; m=100, n=100, shifttype=:tie, list_skewness=list_skewness2)

# %%
@time print_sim(gammadist(1.1, 3), gammadist(1, 3), 1000, 1000; shifttype=:mean, L=10^4)

# %%
@time print_sim(gammadist(1.1, 3), gammadist(1, 3), 1000, 1000; shifttype=:median, L=10^4)

# %%
@time print_sim(gammadist(1.1, 3), gammadist(1, 3), 1000, 1000; shifttype=:tie, L=10^4)

# %%
@time print_sim(gammadist(2, 3), gammadist(1, 3), 1000, 1000; shifttype=:tie, L=10^4)

# %%
@time print_sim(gammadist(5, 5), gammadist(1, 5), 1000, 1000; shifttype=:tie, L=10^4)

# %% [markdown]
# ## その他の計算

# %%
@time print_sim(Normal(0, 4), Normal(0, 1), 100, 100; shifttype=:mean, L=10^6)

# %%
@time print_sim(Normal(0, 4), Normal(0, 1), 100, 200; shifttype=:mean, L=10^6)

# %%
gammadist2(μ, σ) = Gamma(μ^2/σ^2, σ^2/μ)
[
    [shape(gammadist2(3, σ)) for σ in (√3, √3/2, √3/3, √3/4)],
    [scale(gammadist2(3, σ)) for σ in (√3, √3/2, √3/3, √3/4)],
    [mean(gammadist2(3, σ)) for σ in (√3, √3/2, √3/3, √3/4)],
    [var(gammadist2(3, σ)) for σ in (√3, √3/2, √3/3, √3/4)],
    [std(gammadist2(3, σ)) for σ in (√3, √3/2, √3/3, √3/4)]
]

# %%
# 等母平均の場合: この場合はWMW検定とBM検定にとって不利な状況になっている.

for σ in (√3, √3/2, √3/3, √3/4)
    println("-"^20, " 標準治療側の標準偏差 = $(round(√3; digits=2)),  試験治療側の標準偏差 = $(round(σ; digits=2))")
    print_sim(gammadist2(3, √3), gammadist2(3, σ), 100, 100; shifttype=:mean)
end

# %%
# フェアな比較: WMWとBMでは P(X<Y)+P(X=Y)/2 = 1/2 という状況で確認を行う.

for σ in (√3, √3/2, √3/3, √3/4)
    println("-"^20, " 標準治療側の標準偏差 = $(round(√3; digits=2)),  試験治療側の標準偏差 = $(round(σ; digits=2))")
    print_sim(gammadist2(3, √3), gammadist2(3, σ), 100, 100; shifttype=:auto)
end

# %%
