---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.10.3
  kernelspec:
    display_name: Julia 1.7.3
    language: julia
    name: julia-1.7
---

# Brunner-Munzel検定について

* 黒木玄
* 2022-08-05

__文献__

* E. Brunner and U. Munzel. The nonparametric Behrens-Fisher problem: Asymptotic theory and a small-sample
approximation. Biometrical Journal, 42:17–25, 2000.
\[[pdf](https://www.researchgate.net/profile/Edgar-Brunner/publication/264799502_Nonparametric_Hypotheses_and_Rank_Statistics_for_Unbalanced_Factorial_Designs/links/5756a00408ae155a87bc5c8c/Nonparametric-Hypotheses-and-Rank-Statistics-for-Unbalanced-Factorial-Designs.pdf)\]

* Karin Neubert and Edgar Brunner, A studentized permutation test for the non-parametric Behrens-Fisher problem, Computational Statistics and Data Analysis, Vol. 51, pp. 5192-5204 (2007).
https://doi.org/10.1016/j.csda.2006.05.024

* Claus P. Nowak, Markus Pauly, Edgar Brunner. The nonparametric Behrens-Fisher problem in small samples.
https://arxiv.org/abs/2208.01231

<!-- #region toc=true -->
<h1>目次<span class="tocSkip"></span></h1>
<div class="toc"><ul class="toc-item"><li><span><a href="#準備" data-toc-modified-id="準備-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>準備</a></span><ul class="toc-item"><li><span><a href="#パッケージの読み込みなど" data-toc-modified-id="パッケージの読み込みなど-1.1"><span class="toc-item-num">1.1&nbsp;&nbsp;</span>パッケージの読み込みなど</a></span></li><li><span><a href="#組み合わせの生成子" data-toc-modified-id="組み合わせの生成子-1.2"><span class="toc-item-num">1.2&nbsp;&nbsp;</span>組み合わせの生成子</a></span></li><li><span><a href="#Welchのt検定" data-toc-modified-id="Welchのt検定-1.3"><span class="toc-item-num">1.3&nbsp;&nbsp;</span>Welchのt検定</a></span></li><li><span><a href="#単峰型の函数が正の値になる場所を見つける函数" data-toc-modified-id="単峰型の函数が正の値になる場所を見つける函数-1.4"><span class="toc-item-num">1.4&nbsp;&nbsp;</span>単峰型の函数が正の値になる場所を見つける函数</a></span></li><li><span><a href="#2つの分布が「互角」になるシフトの仕方を求める函数" data-toc-modified-id="2つの分布が「互角」になるシフトの仕方を求める函数-1.5"><span class="toc-item-num">1.5&nbsp;&nbsp;</span>2つの分布が「互角」になるシフトの仕方を求める函数</a></span></li></ul></li><li><span><a href="#Brunner-Munzel検定" data-toc-modified-id="Brunner-Munzel検定-2"><span class="toc-item-num">2&nbsp;&nbsp;</span>Brunner-Munzel検定</a></span><ul class="toc-item"><li><span><a href="#Brunner-Munzel検定の実装" data-toc-modified-id="Brunner-Munzel検定の実装-2.1"><span class="toc-item-num">2.1&nbsp;&nbsp;</span>Brunner-Munzel検定の実装</a></span></li><li><span><a href="#よく使われているっぽいテストデータで正しく実装されているかを確認" data-toc-modified-id="よく使われているっぽいテストデータで正しく実装されているかを確認-2.2"><span class="toc-item-num">2.2&nbsp;&nbsp;</span>よく使われているっぽいテストデータで正しく実装されているかを確認</a></span></li><li><span><a href="#組み合わせの生成子" data-toc-modified-id="組み合わせの生成子-2.3"><span class="toc-item-num">2.3&nbsp;&nbsp;</span>組み合わせの生成子</a></span></li><li><span><a href="#Brunner-Munzel検定のpermutation版の実装" data-toc-modified-id="Brunner-Munzel検定のpermutation版の実装-2.4"><span class="toc-item-num">2.4&nbsp;&nbsp;</span>Brunner-Munzel検定のpermutation版の実装</a></span></li><li><span><a href="#permutation版が正しく実装されているかの確認" data-toc-modified-id="permutation版が正しく実装されているかの確認-2.5"><span class="toc-item-num">2.5&nbsp;&nbsp;</span>permutation版が正しく実装されているかの確認</a></span></li></ul></li><li><span><a href="#計算例" data-toc-modified-id="計算例-3"><span class="toc-item-num">3&nbsp;&nbsp;</span>計算例</a></span></li><li><span><a href="#Brunner-Munzel検定とWelchのt検定の比較" data-toc-modified-id="Brunner-Munzel検定とWelchのt検定の比較-4"><span class="toc-item-num">4&nbsp;&nbsp;</span>Brunner-Munzel検定とWelchのt検定の比較</a></span><ul class="toc-item"><li><span><a href="#第一種の過誤の確率" data-toc-modified-id="第一種の過誤の確率-4.1"><span class="toc-item-num">4.1&nbsp;&nbsp;</span>第一種の過誤の確率</a></span></li><li><span><a href="#Brunner-Munzel検定は中央値に関する検定ではないことの証拠" data-toc-modified-id="Brunner-Munzel検定は中央値に関する検定ではないことの証拠-4.2"><span class="toc-item-num">4.2&nbsp;&nbsp;</span>Brunner-Munzel検定は中央値に関する検定ではないことの証拠</a></span></li><li><span><a href="#BM検定による互角シフトの信頼区間とWelchのt検定による平均の差の信頼区間の比較" data-toc-modified-id="BM検定による互角シフトの信頼区間とWelchのt検定による平均の差の信頼区間の比較-4.3"><span class="toc-item-num">4.3&nbsp;&nbsp;</span>BM検定による互角シフトの信頼区間とWelchのt検定による平均の差の信頼区間の比較</a></span></li></ul></li><li><span><a href="#小サンプルでのpermutation版の検定とBM検定とWelchのt検定の比較" data-toc-modified-id="小サンプルでのpermutation版の検定とBM検定とWelchのt検定の比較-5"><span class="toc-item-num">5&nbsp;&nbsp;</span>小サンプルでのpermutation版の検定とBM検定とWelchのt検定の比較</a></span></li></ul></div>
<!-- #endregion -->

## 準備


### パッケージの読み込みなど

```julia
using Base.Threads
using BenchmarkTools
using Distributions
using PrettyPrinting
using QuadGK
using Random
using RCall
using Roots
using StatsBase
using StatsFuns
using StatsPlots
default(fmt=:png, size=(400, 250),
    titlefontsize=10, guidefontsize=8, tickfontsize=6)

x ⪅ y = x < y || x ≈ y
x ⪆ y = x > y || x ≈ y
safemul(x, y) = x == 0 ? x : isinf(x) ? typeof(x)(Inf) : x*y
safediv(x, y) = x == 0 ? x : isinf(y) ? zero(y) : x/y
```

### 組み合わせの生成子

```julia
"""
    nextcombination!(n, t, c = typeof(t)[min(t-1, i) for i in 1:t])

`[1,2,…,n]` からの重複無しの `t` 個の組み合わせ `c` をすべて生成したい.

`nextcombination!(n, t, c)` は配列で表現された組み合わせ `c` をその次の組み合わせに書き換えて, `c` を返す.

初期条件を `c = typeof(t)[min(t-1, i) for i in 1:t]` にすると, `binomial(n, t)` 回の `nextcombination!(n, t, c)` ですべての組み合わせが生成される.
"""
function nextcombination!(n, t, c = typeof(t)[min(t-1, i) for i in 1:t])
    t == 0 && return c
    @inbounds for i in t:-1:1
        c[i] += 1
        c[i] > (n - (t - i)) && continue
        for j in i+1:t
            c[j] = c[j-1] + 1
        end
        break
    end
    c
end

"""
    mycombinations!(n::Integer, t, c)

事前に割り当てられた組み合わせを格納する配列 `c` を使って, `[1,2,…,n]` からの重複無しの `t` 個の組み合わせのすべてを生成する生成子を返す.
"""
function mycombinations!(n::Integer, t, c)
    for i in 1:t c[i] = min(t - 1, i) end
    (nextcombination!(n, t, c) for _ in 1:binomial(n, t))
end

"""
    mycombinations!(a, t, c)

事前に割り当てられた組み合わせを格納する配列 `c` を使って, 配列 `a` からのインデックスに重複がない `t` 個の組み合わせのすべてを生成する生成子を返す.
"""
function mycombinations!(a, t, c)
    t < 0 && (t = length(a) + 1)
    (view(a, indices) for indices in mycombinations!(length(a), t, c))
end

"""
    mycombinations(x, t)

`x` が整数ならば `[1,2,…,x]` からの, `x` が配列ならば `x` からのインデックスに重複がない `t` 個の組み合わせのすべてを生成する生成子を返す.
"""
mycombinations(x, t) = mycombinations!(x, t, Vector{typeof(t)}(undef, t))
```

### Welchのt検定

```julia
function tvalue_welch(m, x̄, sx², n, ȳ, sy²; Δμ=0)
    (x̄ - ȳ - Δμ) / √(sx²/m + sy²/n)
end

function tvalue_welch(x, y; Δμ=0)
    m, x̄, sx² = length(x), mean(x), var(x)
    n, ȳ, sy² = length(y), mean(y), var(y)
    tvalue_welch(m, x̄, sx², n, ȳ, sy²; Δμ)
end

function degree_of_freedom_welch(m, sx², n, sy²)
    (sx²/m + sy²/n)^2 / ((sx²/m)^2/(m-1) + (sy²/n)^2/(n-1))
end

function degree_of_freedom_welch(x, y)
    m, sx² = length(x), var(x)
    n, sy² = length(y), var(y)
    degree_of_freedom_welch(m, sx², n, sy²)
end

function pvalue_welch(m, x̄, sx², n, ȳ, sy²; Δμ=0)
    t = tvalue_welch(m, x̄, sx², n, ȳ, sy²; Δμ)
    ν = degree_of_freedom_welch(m, sx², n, sy²)
    2ccdf(TDist(ν), abs(t))
end

function pvalue_welch(x, y; Δμ=0)
    m, x̄, sx² = length(x), mean(x), var(x)
    n, ȳ, sy² = length(y), mean(y), var(y)
    pvalue_welch(m, x̄, sx², n, ȳ, sy²; Δμ)
end

function confint_welch(m, x̄, sx², n, ȳ, sy²; α=0.05)
    ν = degree_of_freedom_welch(m, sx², n, sy²)
    c = quantile(TDist(ν), 1-α/2)
    SEhat = √(sx²/m + sy²/n)
    [x̄-ȳ-c*SEhat, x̄-ȳ+c*SEhat]
end

function confint_welch(x, y; α=0.05)
    m, x̄, sx² = length(x), mean(x), var(x)
    n, ȳ, sy² = length(y), mean(y), var(y)
    confint_welch(m, x̄, sx², n, ȳ, sy²; α)
end
```

### 単峰型の函数が正の値になる場所を見つける函数

```julia
function findpositive(f, a, b; maxsplit = 30)
    @assert f(a) < 0
    @assert f(b) < 0
    c = (a + b)/2
    f(c) > 0 && return c
    w = b - a
    for k in 2:maxsplit
        for d in range(w/2^(k+1), w/2-w/2^(k+1), step=w/2^k)
            x = c + d
            f(x) > 0 && return x 
            x = c - d
            f(x) > 0 && return x 
        end
    end
    error("k > maxplit = $maxsplit")
end

f(x) = abs(x) < 1e-4 ? 1.0 : -1.0

@time findpositive(f, -100abs(randn()), 20abs(randn()))
```

### 2つの分布が「互角」になるシフトの仕方を求める函数

```julia
"""
    prob_x_le_y(distx::UnivariateDistribution, disty::UnivariateDistribution;
        a = 0.0)

この函数は, 連続分布 `distx`, `disty` と実数 `a` について, 
`distx` と `disty` に従って生成される乱数をそれぞれ X, Y と書くとき, 
X ≤ Y + a が成立する確率を返す.
"""
function prob_x_le_y(distx::UnivariateDistribution, disty::UnivariateDistribution,
        a = 0.0)
    H(y) = cdf(distx, y) * pdf(disty, y-a)
    quadgk(H, extrema(disty + a)...)[1]
end

"""
    tieshift(distx::UnivariateDistribution, disty::UnivariateDistribution;
        p = 0.5)

この函数は, 連続分布 `distx`, `disty` と実数 `p` について, 
`distx` と `disty` に従って生成される乱数をそれぞれ X, Y と書くとき, 
X ≤ Y + a が成立する確率が `p` に等しくなるような実数 a を返す.
"""
function tieshift(distx::UnivariateDistribution, disty::UnivariateDistribution;
        p=0.5)
    find_zero(a -> prob_x_le_y(distx, disty, a) - p, 0.0)
end

@show tieshift(Normal(0, 1), Normal(2, 2))
@show tieshift(Normal(0, 1), Laplace(2, 2))
@show tieshift(Normal(0, 1), Uniform(0, 1));
```

## Brunner-Munzel検定


### Brunner-Munzel検定の実装

````julia
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
````

```julia
@doc h_brunner_munzel
```

```julia
@doc phat_brunner_munzel
```

```julia
@doc statistics_brunner_munzel
```

```julia
@doc brunner_munzel
```

```julia
X = randn(10)
Y = randn(10)
@show shiftmin = minimum(X) - maximum(Y) - 1
@show shiftmax = maximum(X) - minimum(Y) + 1
pvalue_brunner_munzel(X, Y)
```

```julia
plot(a -> pvalue_brunner_munzel(X, Y .+ a), shiftmin, shiftmax;
    label="", title="P-value function of shift")
```

### よく使われているっぽいテストデータで正しく実装されているかを確認

<!-- #region -->
https://okumuralab.org/~okumura/stat/brunner-munzel.html

```R
x = c(1,2,1,1,1,1,1,1,1,1,2,4,1,1)
y = c(3,3,4,3,1,2,3,1,1,5,4)
brunnermunzel.test(x, y)

data:  x and y
Brunner-Munzel Test Statistic = 3.1375, df = 17.683, p-value = 0.005786
95 percent confidence interval:
 0.5952169 0.9827052
sample estimates:
P(X<Y)+.5*P(X=Y) 
        0.788961 
```
<!-- #endregion -->

```julia
X = [1,2,1,1,1,1,1,1,1,1,2,4,1,1]
Y = [3,3,4,3,1,2,3,1,1,5,4]
show_plot_brunner_munzel(X, Y)
```

```julia
X = [1,2,1,1,1,1,1,1,1,1,2,4,1,1]
Y = [3,3,4,3,1,2,3,1,1,5,4]
@rput X Y
R"""
library(lawstat)
brunner.munzel.test(X, Y)
"""
```

このように Brunner-Munzel 検定は R では [lawstat](https://cran.r-project.org/package=lawstat) パッケージの [brunner.munzel.test](https://rdrr.io/cran/lawstat/man/brunner.munzel.test.html) で使える.


### 組み合わせの生成子

```julia
"""
    complementcomb!(complcomb::AbstractVector, comb::AbstractVector)

`comb` が {1,2,…,N} から重複無しに m 個を選ぶ組み合わせを表す配列であり, `comb` の中で数は小さな順に並んでいるとし, `complcomb` は長さ N - m の配列であると仮定する.

このとき, この函数は配列 `complcomb` に配列 `comb` の補集合を格納し, `complcomb` を返す.

この函数はメモリ割り当てゼロで実行される.
"""
function complementcomb!(complcomb::AbstractVector, comb::AbstractVector)
    N = length(comb) + length(complcomb)
    k = 0
    a = 0
    @inbounds for b in comb
        for i in a+1:b-1
            k += 1
            complcomb[k] = i
        end
        a = b
    end
    @inbounds for i in a+1:N
        k +=1
        complcomb[k] = i
    end
    complcomb
end

"""
    complementcomb(N, comb::AbstractVector)

`comb` が {1,2,…,N} から重複無しに m 個を選ぶ組み合わせを表す配列であり, `comb` の中で数は小さな順に並んでいると仮定する.

この函数は `comb` の補集合の配列を返す.

この函数は返り値の配列の分だけのメモリ割り当てを行う.
"""
complementcomb(N, comb::AbstractVector) =
    complementcomb!(similar(comb, N - length(comb)), comb)
```

```julia
@doc complementcomb!
```

```julia
@doc complementcomb
```

```julia
N = 10
comb = [2, 4, 5, 8]
ccomb = similar(comb, N - length(comb))
@btime complementcomb!($ccomb, $comb);
```

```julia
N, m = 5, 3
ccomb = Vector{Int}(undef, N-m)
[(copy(comb), copy(complementcomb!(ccomb, comb))) for comb in mycombinations(1:N, m)]
```

```julia
N, m = 5, 3
ccomb = Vector{Int}(undef, N-m)
[(copy(comb), complementcomb(N, comb)) for comb in mycombinations(1:N, m)]
```

### Brunner-Munzel検定のpermutation版の実装

```julia
"""
    permutation_tvalues_brunner_munzel(X, Y,
        XandY = Vector{Float64}(undef, length(X)+length(Y)),
        Tval = Vector{Float64}(undef, binomial(length(X)+length(Y), length(X))),
        Hx = similar(X, Float64),
        Hy = similar(Y, Float64)
    )

Brunner-Munzel検定のt値を `[X; Y]` から\
インデックスの重複無しに `length(X)` 個取る組み合わせと\
その補集合への分割のすべてについて計算して, `Tval` に格納して返す.
"""
function permutation_tvalues_brunner_munzel(X, Y,
        XandY = Vector{Float64}(undef, length(X)+length(Y)),
        Tval = Vector{Float64}(undef, binomial(length(X)+length(Y), length(X))),
        Hx = similar(X, Float64),
        Hy = similar(Y, Float64),
        ccomb = Vector{Int}(undef, length(Y))
    )
    m, n = length(X), length(Y)
    N = m + n
    @views XandY[1:m] .= X
    @views XandY[m+1:N] .= Y
    for (k, comb) in enumerate(mycombinations(1:N, m))
        complementcomb!(ccomb, comb)
        Tval[k] = statistics_brunner_munzel(
            view(XandY, comb), view(XandY, ccomb), Hx, Hy).tvalue
    end
    Tval
end

"""
    pvalue_brunner_munzel_perm(X, Y,
        Tval = permutation_tvalues_brunner_munzel(X, Y),
        tval = statistics_brunner_munzel(X, Y).tvalue;
        le = ⪅
    )

Brunner-Munzel検定のpermutation版のP値を返す.
"""
function pvalue_brunner_munzel_perm(X, Y,
        Tval = permutation_tvalues_brunner_munzel(X, Y),
        tval = statistics_brunner_munzel(X, Y).tvalue;
        le = ⪅
    )
    pvalue_perm = mean(T -> le(abs(tval), abs(T)), Tval)
end
```

```julia
@doc permutation_tvalues_brunner_munzel
```

```julia
@doc pvalue_brunner_munzel_perm
```

https://okumuralab.org/~okumura/stat/brunner-munzel.html

```
bm = brunner.munzel.test(x, y)$statistic
n1 = length(x)
n2 = length(y)
N = n1 + n2
xandy = c(x, y)
foo = function(X) {
  brunner.munzel.test(xandy[X], xandy[-X])$statistic
}
z = combn(1:N, n1, foo)
mean(abs(z) >= abs(bm))
```

>結果は 0.008037645 となりました。

```julia
X = [1,2,1,1,1,1,1,1,1,1,2,4,1,1]
Y = [3,3,4,3,1,2,3,1,1,5,4]
@show X Y
@show m, n = length(X), length(Y)

Tval = @time permutation_tvalues_brunner_munzel(X, Y)
@show pvalue_brunner_munzel_perm(X, Y, Tval)
stephist(Tval; norm=true, bin=101, label="", title="permutation t-values")
```

### permutation版が正しく実装されているかの確認

* https://github.com/toshi-ara/brunnermunzel/issues/14
* https://github.com/toshi-ara/brunnermunzel/files/4395032/mwe.R.zip

<!-- #region -->
__追記 2022-08-06:__ https://twitter.com/TA25140989/status/1555825941451923457 を参照せよ.

https://github.com/toshi-ara/brunnermunzel/tree/development の修正版の brunnermunzel パッケージをインストールし直した.  以下のセルの実行結果が変わるはずなので, その記録を残しておく.

以下の記録を見なくても,

* https://github.com/genkuroki/public/blob/e4faafc52721b63876b3b705f9450eade3c902f5/0034/Brunner-Munzel.ipynb

で閲覧できるが, わざわざ見に行くのも面倒なのでこのファイルにも記録を残しておく.

以前の実行結果:

```julia
@show pval_J - pval_J_le;
```
```
pval_J - pval_J_le = [0.0, 0.0, 0.007936507936507936, 0.023809523809523808, 0.023809523809523808, 0.0793650793650793, 0.0, 0.0, 0.007936507936507936, 0.0357142857142857, 0.0, 0.015873015873015928, 0.015873015873015872, 0.0, 0.015873015873015872, 0.0, 0.0, 0.0, 0.007936507936507936, 0.0, 0.023809523809523836, 0.015873015873015928, 0.011904761904761918, 0.0, 0.015873015873015928, 0.0, 0.003968253968253968, 0.0, 0.0, 0.0, 0.0, 0.007936507936507936, 0.0, 0.015873015873015928, 0.0, 0.015873015873015928, 0.0, 0.0, 0.0, 0.007936507936507908, 0.0, 0.007936507936507908, 0.0, 0.007936507936507908, 0.0, 0.0, 0.007936507936507936, 0.007936507936507936, 0.011904761904761918, 0.007936507936507936, 0.03968253968253965, 0.0, 0.007936507936507936, 0.0, 0.05555555555555547, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.015873015873015872, 0.015873015873015872, 0.0, 0.015873015873015928, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.015873015873015928, 0.007936507936507964, 0.007936507936507964, 0.0, 0.003968253968253954, 0.0, 0.0, 0.0, 0.015873015873015872, 0.0, 0.0, 0.0357142857142857, 0.015873015873015872, 0.0, 0.023809523809523725, 0.0, 0.023809523809523808, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
```
```julia
idx = @show findall(pval_J .!= pval_J_le)
length(idx)
```
```
findall(pval_J .!= pval_J_le) = [3, 4, 5, 6, 9, 10, 12, 13, 15, 19, 21, 22, 23, 25, 27, 32, 34, 36, 40, 42, 44, 47, 48, 49, 50, 51, 53, 55, 68, 69, 71, 78, 79, 80, 82, 86, 89, 90, 92, 94]

40
```
```julia
@show pval_R - pval_J_le;
```
```
pval_R - pval_J_le = [0.0, 0.0, 0.007936507936507936, 0.023809523809523808, 0.023809523809523808, 0.0793650793650793, 0.0, 0.0, 0.007936507936507936, 0.0357142857142857, 0.0, 0.015873015873015928, 0.015873015873015872, 0.0, 0.015873015873015872, 0.0, 0.0, 0.0, 0.007936507936507936, 0.0, 0.023809523809523836, 0.007936507936507908, 0.011904761904761918, 0.0, 0.015873015873015928, 0.0, 0.003968253968253968, 0.0, 0.0, 0.0, 0.0, 0.007936507936507936, 0.0, 0.00793650793650802, 0.0, 0.015873015873015928, 0.0, 0.0, 0.0, 0.007936507936507908, 0.0, 0.007936507936507908, 0.0, 0.007936507936507908, 0.0, 0.0, 0.007936507936507936, 0.007936507936507936, 0.011904761904761918, 0.007936507936507936, 0.03968253968253965, 0.0, 0.007936507936507936, 0.0, 0.05555555555555547, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.015873015873015872, 0.015873015873015872, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.007936507936507908, 0.007936507936507964, 0.007936507936507964, 0.0, 0.003968253968253954, 0.0, 0.0, 0.0, 0.015873015873015872, 0.0, 0.0, 0.0357142857142857, 0.015873015873015872, 0.0, 0.023809523809523725, 0.0, 0.023809523809523808, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
```
```julia
idx = @show findall(pval_R .!= pval_J_le)
length(idx)
```
```
findall(pval_R .!= pval_J_le) = [3, 4, 5, 6, 9, 10, 12, 13, 15, 19, 21, 22, 23, 25, 27, 32, 34, 36, 40, 42, 44, 47, 48, 49, 50, 51, 53, 55, 68, 69, 78, 79, 80, 82, 86, 89, 90, 92, 94]

39
```
```julia
@show pval_J - pval_R;
```
```
pval_J - pval_R = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00793650793650802, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.007936507936507908, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.015873015873015928, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.00793650793650802, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
```
```julia
idx = @show findall(.!(pval_J .≈ pval_R))
length(idx)
```
```
findall(.!(pval_J .≈ pval_R)) = [22, 34, 71, 78]

4
```
```julia
all(pval_J .≥ pval_R .≥ pval_J_le)
```
```
true
```

以前のコメントは以下の通り:

```
なるほど！

https://github.com/toshi-ara/brunnermunzel/issues/14

に書いてあるように, 22, 34, 71 and 78 の4つで, 値が一致していない.

〇〇以下または〇〇以上の判定を $x \approx y$ のときも true にする必要があるのだが, その部分で違いが生じているものと思われる.

現時点では https://CRAN.R-project.org/package=brunnermunzel  にアクセスすると,

>Package ‘brunnermunzel’ was removed from the CRAN repository.
>
>Formerly available versions can be obtained from the [archive](https://cran.r-project.org/src/contrib/Archive/brunnermunzel/).
>
>Archived on 2022-03-04 as check problems were not corrected in time. , LENGTH_1 checks.
>
>A summary of the most recent check results can be obtained from the [check results archive](https://cran-archive.r-project.org/web/checks/2022/2022-03-04_check_results_brunnermunzel.html).
>
>Please use the canonical form https://CRAN.R-project.org/package=brunnermunzel to link to this page.

と表示される.
```
<!-- #endregion -->

```julia
R"""
library(brunnermunzel)
set.seed(1290)
reps = 100
xx = c()
yy = c()
pval_R = numeric(reps)
for (i in seq_len(reps)){
  x = rnorm(5)
  y = rnorm(5)
  
  xx = c(xx, x)
  yy = c(yy, y)
  
  res_bm_perm <- brunnermunzel.permutation.test(x,y)
  pval_R[i] <- res_bm_perm$p.value
}
"""

@rget xx yy pval_R
XX = reshape(xx, 5, 100)
YY = reshape(yy, 5, 100)

pval_J = zeros(100)
pval_J_le = zeros(100)
for i in 1:100
    pval_J[i]    = pvalue_brunner_munzel_perm(XX[:,i], YY[:,i]; le = ⪅)
    pval_J_le[i] = pvalue_brunner_munzel_perm(XX[:,i], YY[:,i]; le = ≤)
end
```

```julia
@show pval_J - pval_J_le;
```

```julia
idx = @show findall(pval_J .!= pval_J_le)
length(idx)
```

```julia
@show pval_R - pval_J_le;
```

```julia
idx = @show findall(pval_R .!= pval_J_le)
length(idx)
```

```julia
@show pval_J - pval_R;
```

```julia
idx = @show findall(.!(pval_J .≈ pval_R))
length(idx)
```

```julia
all(pval_J .≥ pval_R .≥ pval_J_le)
```

__2022-08-06:__ やった! 値が完全に一致した! permutation版Brunner-Munzel検定について,

* https://github.com/toshi-ara/brunnermunzel/tree/development

の実装と私による実装の計算結果は以上の場合において完全に一致している.


## 計算例

```julia
m, n = 10, 10
X, Y = rand(Normal(0, 1), m), rand(Normal(0, 2), n)
@show pval_brmu = pvalue_brunner_munzel(X, Y)
@show pval_perm = pvalue_brunner_munzel_perm(X, Y);
```

```julia
Random.seed!(4)

m, n = 10, 20
X, Y = rand(Normal(0, 1), m), rand(Normal(0, 2), n)
show_plot_brunner_munzel(X, Y)
```

```julia
m, n = 10, 20
X, Y = rand(Normal(0, 1), m), rand(Normal(0, 2), n)
show_plot_brunner_munzel(X, Y)
```

```julia
m, n = 100, 50
X, Y = rand(Normal(0, 1), m), rand(Normal(0, 2), n)
show_plot_brunner_munzel(X, Y)
```

```julia
Random.seed!(4)

m, n = 10, 20
X, Y = rand(1:m, m), rand(1:n, n)
show_plot_brunner_munzel(X, Y)
```

```julia
m, n = 10, 20
X, Y = rand(1:m, m), rand(1:n, n)
show_plot_brunner_munzel(X, Y)
```

```julia
distx, disty = LogNormal(), LogNormal(1)
m, n = 10, 10
X, Y = rand(distx, m), rand(disty, n)
show_plot_brunner_munzel(X, Y)
```

```julia
distx, disty = LogNormal(), LogNormal(1)
m, n = 40, 40
X, Y = rand(distx, m), rand(disty, n)
show_plot_brunner_munzel(X, Y; xlim=(-5, 2))
```

```julia
distx, disty = LogNormal(), LogNormal(1)
m, n = 160, 160
X, Y = rand(distx, m), rand(disty, n)
show_plot_brunner_munzel(X, Y; xlim=(-3, 0))
```

## Brunner-Munzel検定とWelchのt検定の比較


### 第一種の過誤の確率

```julia
function sim_brunner_mumzel(;
        distx = Normal(0, 1), disty = Normal(0, 2), m = 10, n = 10,
        L = 10^6)
    pval_bm = Vector{Float64}(undef, L)
    tmpX = [Vector{Float64}(undef, m) for _ in 1:nthreads()]
    tmpY = [Vector{Float64}(undef, n) for _ in 1:nthreads()]
    tmpHx = [Vector{Float64}(undef, m) for _ in 1:nthreads()]
    tmpHy = [Vector{Float64}(undef, n) for _ in 1:nthreads()]
    @threads for i in 1:L
        X = rand!(distx, tmpX[threadid()])
        Y = rand!(disty, tmpY[threadid()])
        pval_bm[i] = pvalue_brunner_munzel(X, Y, tmpHx[threadid()], tmpHy[threadid()])
    end
    ecdf(pval_bm)
end

function sim_welch(;
        distx = Normal(0, 1), disty = Normal(0, 2), m = 10, n = 10,
        L = 10^6)
    pval_w = Vector{Float64}(undef, L)
    tmpX = [Vector{Float64}(undef, m) for _ in 1:nthreads()]
    tmpY = [Vector{Float64}(undef, n) for _ in 1:nthreads()]
    @threads for i in 1:L
        X = rand!(distx, tmpX[threadid()])
        Y = rand!(disty, tmpY[threadid()])
        pval_w[i] = pvalue_welch(X, Y)
    end
    ecdf(pval_w)
end

function printcompact(io, xs...)
    print(IOContext(io, :compact => true), xs...)
end

function distname(dist)
    replace(sprint(printcompact, dist), r"\{[^\}]*\}"=>"")
end

function plot_ecdf(ecdf_pval, distx, disty, m, n, a;
        testname = "", kwargs...)
    plot(p -> ecdf_pval(p), 0, 0.1; label="ecdf of P-values")
    plot!([0, 0.1], [0, 0.1]; label="", ls=:dot, c=:black)
    plot!(legend=:topleft)
    plot!(xtick=0:0.01:0.1, ytick=0:0.01:1)
    plot!(xguide="nominal significance level α", 
        yguide="probability of P-value < α")
    s = (a < 0 ? "" : "+") * string(round(a; digits=4))
    title!("$(testname)X: $(distname(distx)), m=$m\n\
        Y: $(distname(disty))$s, n=$n")
    plot!(size=(400, 450))
    plot!(; kwargs...)
end

function plot_pvals(;
        distx = Normal(0, 1), disty = Normal(0, 2), m = 10, n = 10,
        L = 10^6, a = nothing, Δμ = nothing, kwargs...)
    @show (mean(distx), std(distx))
    @show (mean(disty), std(disty))
    
    if isnothing(a)
        @show a = tieshift(distx, disty)
        @show prob_x_le_y(distx, disty + a)
    else
        @show a
        @show median(distx) - median(disty)
    end
    if isnothing(Δμ)
        @show Δμ = mean(distx) - mean(disty)
        @show mean(distx), mean(disty + Δμ)
    else
        @show Δμ
        @show mean(distx), mean(disty + Δμ)
    end
        
    ecdf_bm = @time sim_brunner_mumzel(;
        distx = distx,
        disty = disty + a,
        m, n, L, kwargs...)
    ecdf_w = @time sim_welch(;
        distx = distx,
        disty = disty + Δμ,
        m, n, L, kwargs...)
    ymax = max(ecdf_bm(0.1), ecdf_w(0.1))
    P1 = plot_ecdf(ecdf_bm, distx, disty, m, n, a;
        testname="Brunner-Munzel test\n",
        ylim=(-0.002, 1.02*ymax), kwargs...)
    P2 = plot_ecdf(ecdf_w, distx, disty, m, n, Δμ;
        testname="Welch t-test\n",
        ylim=(-0.002, 1.02*ymax), kwargs...)
    plot(P1, P2; size=(800, 450), topmargin=3.5Plots.mm)
end
```

```julia
plot_pvals(; distx = Normal(0, 1), disty = Normal(0, 1), m = 10, n = 10)
```

```julia
plot_pvals(; distx = Normal(0, 1), disty = Normal(0, 2), m = 10, n = 10)
```

```julia
plot_pvals(; distx = LogNormal(), disty = LogNormal(1), m = 10, n = 10)
```

```julia
plot_pvals(; distx = LogNormal(), disty = LogNormal(1), m = 20, n = 20)
```

```julia
plot_pvals(; distx = LogNormal(), disty = LogNormal(1), m = 40, n = 40)
```

```julia
plot_pvals(; distx = TDist(2), disty = TDist(2), m = 10, n = 10, Δμ = 0.0)
```

```julia
plot_pvals(; distx = TDist(2), disty = TDist(1.1), m = 10, n = 10, Δμ = 0.0)
```

### Brunner-Munzel検定は中央値に関する検定ではないことの証拠

```julia
distx, disty = Uniform(-1, 1), Exponential()
m, n, = 100, 100

@show distx, std(distx)
@show disty, std(disty)

@show a = tieshift(distx, disty)
ecdf_pval1 = @time sim_brunner_mumzel(;
    distx = distx, disty = disty + a, m, n)
P1 = plot_ecdf(ecdf_pval1, distx, disty, m, n, a;
    testname="case of tie shifting\n")

@show a = median(distx) - median(disty)
ecdf_pval2 = @time sim_brunner_mumzel(;
    distx = distx, disty = disty + a, m, n)
P2 = plot_ecdf(ecdf_pval2, distx, disty, m, n, a;
    testname="case of matching medians\n")

plot(P1, P2; size=(800, 450), topmargin=4Plots.mm)
```

```julia
distx, disty = Uniform(-1, 1), Exponential()
@show distx, std(distx)
@show disty, std(disty)

a = @show tieshift(distx, disty)
P1 = plot(distx, -2, 6; label="distx")
plot!(disty + a, -2, 6; label="disty + ($(round(a; digits=4)))", ls=:dash)
title!("case of tie shifting")

a = @show median(distx) - median(disty)
P2 = plot(distx, -2, 6; label="distx")
plot!(disty + a, -2, 6; label="disty + ($(round(a; digits=4)))", ls=:dash)
vline!([median(distx)]; label="median(distx)", ls=:dot, lw=1.5)
title!("case of matching medians")

plot(P1, P2; size=(800, 250))
```

```julia
distx, disty = Uniform(-1, 1), Exponential(4)
m, n, = 100, 100

@show distx, std(distx)
@show disty, std(disty)

@show a = tieshift(distx, disty)
ecdf_pval1 = @time sim_brunner_mumzel(;
    distx = distx, disty = disty + a, m, n)
P1 = plot_ecdf(ecdf_pval1, distx, disty, m, n, a;
    testname="case of tie shifting\n")

@show a = median(distx) - median(disty)
ecdf_pval2 = @time sim_brunner_mumzel(;
    distx = distx, disty = disty + a, m, n)
P2 = plot_ecdf(ecdf_pval2, distx, disty, m, n, a;
    testname="case of matching medians\n")

plot(P1, P2; size=(800, 450), topmargin=4Plots.mm)
```

```julia
distx, disty = Uniform(-1, 1), Exponential(4)
@show distx, std(distx)
@show disty, std(disty)

a = @show tieshift(distx, disty)
P1 = plot(distx, -4, 10; label="distx")
plot!(disty + a, -4, 10; label="disty + ($(round(a; digits=4)))", ls=:dash)
title!("case of tie shifting")

a = @show median(distx) - median(disty)
P2 = plot(distx, -4, 10; label="distx")
plot!(disty + a, -4, 10; label="disty + ($(round(a; digits=4)))", ls=:dash)
vline!([median(distx)]; label="median(distx)", ls=:dot, lw=1.5)
title!("case of matching medians")

plot(P1, P2; size=(800, 250))
```

```julia
distx, disty = Uniform(-1, 1), Exponential(0.5773502691896257)
m, n, = 100, 100

@show distx, std(distx)
@show disty, std(disty)

@show a = tieshift(distx, disty)
ecdf_pval1 = @time sim_brunner_mumzel(;
    distx = distx, disty = disty + a, m, n)
P1 = plot_ecdf(ecdf_pval1, distx, disty, m, n, a;
    testname="case of tie shifting\n")

@show a = median(distx) - median(disty)
ecdf_pval2 = @time sim_brunner_mumzel(;
    distx = distx, disty = disty + a, m, n)
P2 = plot_ecdf(ecdf_pval2, distx, disty, m, n, a;
    testname="case of matching medians\n")

plot(P1, P2; size=(800, 450), topmargin=4Plots.mm)
```

```julia
distx, disty = Uniform(-1, 1), Exponential(0.5773502691896257)
@show distx, std(distx)
@show disty, std(disty)

a = @show tieshift(distx, disty)
P1 = plot(distx, -2, 4; label="distx")
plot!(disty + a, -2, 4; label="disty + ($(round(a; digits=4)))", ls=:dash)
title!("case of tie shifting")

a = @show median(distx) - median(disty)
P2 = plot(distx, -2, 4; label="distx")
plot!(disty + a, -2, 4; label="disty + ($(round(a; digits=4)))", ls=:dash)
vline!([median(distx)]; label="median(distx)", ls=:dot, lw=1.5)
title!("case of matching medians")

plot(P1, P2; size=(800, 250))
```

### BM検定による互角シフトの信頼区間とWelchのt検定による平均の差の信頼区間の比較

```julia
function plot_confints(;
        distx = Normal(0, 1), disty = Normal(0, 2), m = 10, n = 10,
        L = 100, kwargs...)
    a = tieshift(distx, disty)
    Δμ = mean(distx) - mean(disty)
    BM = fill(zeros(2), 0)
    W = fill(zeros(2), 0)
    for _ in 1:L
        X = rand(distx, m)
        Y = rand(disty, n)
        push!(BM, brunner_munzel(X, Y .+ a).confint_shift)
        push!(W, confint_welch(X, Y .+ Δμ))
    end
    P = plot()
    for i in 1:L
        plot!(fill(i, 2), [first(BM[i]), last(BM[i])]; label="", c=1, lw=2)
        plot!(fill(i+0.3, 2), [first(W[i]), last(W[i])]; label="", c=2, lw=2)
    end
    title!("X: $(distname(distx)), m=$m,   Y: $(distname(disty)), n=$n")
    plot!(size=(1000, 250))
end
```

```julia
plot_confints(distx = Normal(0, 1), disty = Normal(0, 2), m = 10, n = 10)
```

```julia
plot_confints(distx = Normal(2, 1), disty = Normal(0, 2), m = 10, n = 10)
```

```julia
plot_confints(distx = LogNormal(0), disty = LogNormal(1), m = 10, n = 10)
```

```julia
plot_confints(distx = LogNormal(0), disty = LogNormal(1), m = 40, n = 40)
```

```julia
plot_confints(distx = LogNormal(0), disty = LogNormal(1), m = 160, n = 160)
```

```julia
plot_confints(distx = TDist(2), disty = TDist(2), m = 10, n = 10)
```

```julia
plot_confints(distx = TDist(2), disty = TDist(1.1), m = 10, n = 10)
```

```julia
function plot_limits(;
        distx = Normal(0, 1), disty = Normal(0, 2), m = 10, n = 10,
        L = 1000, kwargs...)
    
    @show distx, m
    @show disty, n

    a = tieshift(distx, disty)
    Δμ = mean(distx) - mean(disty)

    BM = fill(zeros(2), 0)
    W = fill(zeros(2), 0)
    for _ in 1:L
        X = rand(distx, m)
        Y = rand(disty, n)
        push!(BM, brunner_munzel(X, Y .+ a).confint_shift)
        push!(W, confint_welch(X, Y .+ Δμ))
    end

    lower = [(first(BM[i]), first(W[i])) for i in 1:L]
    upper = [(last(BM[i]), last(W[i])) for i in 1:L]

    P1 = scatter(lower; label="", msc=:auto, ms=2, ma=0.5)
    plot!(identity; label="")
    plot!(xguide="Brunner-Munzel", yguide="Welch")
    title!("lower")

    P2 = scatter(upper; label="", msc=:auto, ms=2, ma=0.5)
    plot!(identity; label="")
    plot!(xguide="Brunner-Munzel", yguide="Welch")
    title!("upper")

    plot(P1, P2; size=(640, 320))
end
```

```julia
plot_limits(distx = Normal(0, 1), disty = Normal(0, 2), m = 10, n = 10)
```

```julia
plot_limits(distx = Normal(2, 1), disty = Normal(0, 2), m = 10, n = 10)
```

```julia
plot_limits(distx = LogNormal(), disty = LogNormal(1), m = 10, n = 10)
```

```julia
plot_limits(distx = LogNormal(), disty = LogNormal(1), m = 40, n = 40)
```

```julia
@time plot_limits(distx = LogNormal(), disty = LogNormal(1), m = 160, n = 160)
```

```julia
plot_limits(distx = TDist(2), disty = TDist(2), m = 10, n = 10)
```

```julia
plot_limits(distx = TDist(2), disty = TDist(1.1), m = 10, n = 10)
```

## 小サンプルでのpermutation版の検定とBM検定とWelchのt検定の比較

```julia
function sim_brunner_mumzel_perm(;
        distx = Normal(0, 1), disty = Normal(0, 2), m = 5, n = 5,
        L = 10^2)
    pval_bm_perm = Vector{Float64}(undef, L)
    tmpX = [Vector{Float64}(undef, m) for _ in 1:nthreads()]
    tmpY = [Vector{Float64}(undef, n) for _ in 1:nthreads()]
    tmpXandY = [Vector{Float64}(undef, m+n) for _ in 1:nthreads()]
    tmpTval = [Vector{Float64}(undef, binomial(m+n, m)) for _ in 1:nthreads()]
    tmpHx = [Vector{Float64}(undef, m) for _ in 1:nthreads()]
    tmpHy = [Vector{Float64}(undef, n) for _ in 1:nthreads()]
    tmpccomb = [Vector{Int}(undef, n) for _ in 1:nthreads()]
    @threads for i in 1:L
        tid = threadid()
        X = rand!(distx, tmpX[tid])
        Y = rand!(disty, tmpY[tid])
        Tval = permutation_tvalues_brunner_munzel(X, Y,
            tmpXandY[tid], tmpTval[tid], tmpHx[tid], tmpHy[tid], tmpccomb[tid])
        tval = statistics_brunner_munzel(X, Y, tmpHx[tid], tmpHy[tid]).tvalue
        pval_bm_perm[i] = pvalue_brunner_munzel_perm(X, Y, Tval, tval)
    end
    ecdf(pval_bm_perm)
end
```

```julia
@time ecdf_bm_perm = sim_brunner_mumzel_perm(
    distx = Normal(0, 1), disty = Normal(0, 2), m = 7, n = 7, L = 10^4)
```

```julia
function plot_pvals_with_perm(;
        distx = Normal(0, 1),
        disty = Normal(0, 2),
        m = 7,
        n = 7,
        L = 10^4,
        kwargs...
    )
    a = tieshift(distx, disty)
    @time ecdf_bm_perm = sim_brunner_mumzel_perm(; distx, disty = disty + a, m, n, L)
    @time ecdf_bm = sim_brunner_mumzel(; distx, disty = disty + a, m, n, L)
    Δμ = mean(distx) - mean(disty)
    @time ecdf_w = sim_welch(; distx, disty = disty + Δμ, m, n, L)
    @show a Δμ

    plot(legend=:topleft)
    plot!(α -> ecdf_bm_perm(α), 0, 0.1; label="BM permutation")
    plot!(α -> ecdf_bm(α), 0, 0.1; label="Brunner-Munzel", ls=:dash)
    plot!(α -> ecdf_w(α), 0, 0.1; label="Welch", ls=:dashdot)
    plot!(identity; label="", c=:black, ls=:dot)
    plot!(xtick=0:0.01:0.1, ytick=0:0.01:1)
    plot!(xguide="nominal significance level α", 
        yguide="probability of P-value < α")
    a_ = string(round(a; digits=4))
    Δμ_ = string(round(Δμ; digits=4))
    title!("X: $(distname(distx)), m=$m\n\
        Y: $(distname(disty))+(a, Δμ), n=$n\n\
        a=$a_, Δμ=$Δμ_")
    plot!(size=(400, 450), titlefontsize=9)
    plot!(; kwargs...)
end
```

```julia
plot_pvals_with_perm(
    distx = Normal(0, 1), disty = Normal(0, 2), m = 5, n = 5, L = 10^4)
```

```julia
plot_pvals_with_perm(
    distx = Normal(0, 1), disty = Normal(0, 2), m = 7, n = 7, L = 10^4)
```

```julia
plot_pvals_with_perm(
    distx = Normal(0, 1), disty = Normal(0, 2), m = 5, n = 10, L = 10^4)
```

```julia
plot_pvals_with_perm(
    distx = Normal(0, 1), disty = Normal(0, 2), m = 10, n = 5, L = 10^4)
```

```julia
plot_pvals_with_perm(
    distx = Normal(0, 1), disty = Normal(0, 2), m = 10, n = 10, L = 2000)
```

```julia
plot_pvals_with_perm(
    distx = Exponential(1), disty = Exponential(2),
    m = 5, n = 5, L = 10^4)
```

```julia
plot_pvals_with_perm(
    distx = Exponential(1), disty = Exponential(2),
    m = 7, n = 7, L = 10^4)
```

```julia
plot_pvals_with_perm(
    distx = Exponential(1), disty = Exponential(2),
    m = 5, n = 10, L = 10^4)
```

```julia
plot_pvals_with_perm(
    distx = Exponential(1), disty = Exponential(2),
    m = 10, n = 5, L = 10^4)
```

```julia
plot_pvals_with_perm(
    distx = Exponential(1), disty = Exponential(2),
    m = 10, n = 10, L = 2000)
```

```julia
plot_pvals_with_perm(
    distx = Uniform(-1, 1), disty = Exponential(0.5773502691896257),
    m = 5, n = 5, L = 10^4)
```

```julia
plot_pvals_with_perm(
    distx = Uniform(-1, 1), disty = Exponential(0.5773502691896257),
    m = 7, n = 7, L = 10^4)
```

```julia
plot_pvals_with_perm(
    distx = Uniform(-1, 1), disty = Exponential(0.5773502691896257),
    m = 5, n = 10, L = 10^4)
```

```julia
plot_pvals_with_perm(
    distx = Uniform(-1, 1), disty = Exponential(0.5773502691896257),
    m = 10, n = 5, L = 10^4)
```

```julia
plot_pvals_with_perm(
    distx = Uniform(-1, 1), disty = Exponential(0.5773502691896257),
    m = 10, n = 10, L = 2000)
```

```julia
plot_pvals_with_perm(
    distx = LogNormal(), disty = LogNormal(1), m = 5, n = 5, L = 10^4)
```

```julia
plot_pvals_with_perm(
    distx = LogNormal(), disty = LogNormal(1), m = 5, n = 10, L = 10^4)
```

```julia
plot_pvals_with_perm(
    distx = LogNormal(), disty = LogNormal(1), m = 10, n = 5, L = 10^4)
```

```julia
plot_pvals_with_perm(
    distx = LogNormal(0), disty = LogNormal(1), m = 10, n = 10, L = 2000)
```

```julia

```
