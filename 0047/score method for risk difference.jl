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
#     display_name: Julia 1.10.1
#     language: julia
#     name: julia-1.10
# ---

# %% [markdown]
# # 母比率の差に関するP値と信頼区間のスコア法による構成
#
# * 黒木玄
# * 2024-03-02
# $
# \newcommand\tp{\tilde{p}}
# \newcommand\tq{\tilde{q}}
# \newcommand\op{\operatorname}
# \newcommand\pvalue{\op{pvalue}}
# \newcommand\confint{\op{confint}}
# $
#
# 母比率の差に関するスコア法の普及の歴史については
#
# * [割合の差のスコア検定](https://triadsou.hatenablog.com/entry/2024/03/01/001513), Triad sou., 2024-03-01
#
# を参照. 母比率の差のスコア検定の詳しい解説は
#
# * 佐藤俊哉, [コホート研究で用いる疫学指標のefficient scoreにもとづく信頼区間](https://www.jstage.jst.go.jp/article/jappstat1971/17/1/17_1_43/_article/-char/ja/), 1988
#
# にある. 1988年のときには個人で使えるコンピュータ環境は現在と比較するとずっと貧弱で大変だった.  佐藤俊哉さんの論文ではBASICで母比率の差のスコア法が実装されている.

# %% [markdown]
# 母比率の差のP値のスコア方による構成は以下の通り.
#
# __データの型:__ 2つの母集団からそれぞれ $m$ 人と $n$ 人を無作為抽出する. $m$ 人中 $a$ 人はある条件を満たし, 残りの $b=m-a$ 人はその条件を満たさず, $n$ 人中 $c$ 人はその条件を満たし, 残りの $d=n-c$ はその条件を満たさなかったとする.
#
# __統計モデル:__ $a, c$ は独立に成功確率 $p, q$ の二項分布に従ってランダムに生成されているというモデルを考える. 母比率のモデル化 $p,q$ の差を $|Delra=p-a$ とおく. このモデルの確率分布の確率質量函数は次の形になる:
#
# $$
# L := P(a,b,c,d|m,n,\Delta,q) = \binom{m}{a}p^a(1-p)^b\;\binom{n}{c}q^c(1-q)^d.
# $$
#
# ただし,
#
# $$
# a+b=m,\quad c+d=n,\quad p=q+\Delta).
# $$
#
# 以下ではしばらくのあいだ, $a,b,c,d$ は統計モデルに従う確率変数(モデル内確率変数)であるとする. $a,b,c,d$ を観測値を表す記号だとみなす場合には別に断ることにする.  モデル内確率変数と現実で得たデータの数値(観測値)は厳密に区別しなければいけない.

# %% [markdown]
# __スコア統計量:__ 
#
# $$
# U_1(\Delta, q) = (\log L)_\Delta = \frac{\partial}{\partial\Delta}\log L, \quad
# U_2(\Delta, q) = (\log L)_q = \frac{\partial}{\partial q}\log L
# $$
#
# とおく.  このとき, 
#
# $$
# \begin{aligned}
# &
# U_1 := U_1(\Delta, q) = \frac{a}{p}-\frac{b}{1-p},
# \\ &
# U_2 := U_2(\Delta, q) = \frac{a}{p}-\frac{b}{1-p}+\frac{c}{q}-\frac{d}{1-q}.
# \end{aligned}
# $$
#
# ここで　$p=q+\Delta$ であることに注意せよ.
#
# $a,b,c,d$ がモデル内確率変数のとき, 二項分布の期待値より, 
#
# $$
# E[a]=mp, \quad E[b]=m(1-p), \quad E[c]=nq, \quad E[d]=n(1-q)
# $$
#
# なので
#
# $$
# E[U_1] = 0, \quad E[U_2] = 0.
# $$

# %% [markdown]
# __Fisher情報量行列:__ $a,b,c,d$ がモデル内確率変数のときのスコア統計量 $U_1, U_2$ の分散共分散行列を
#
# $$
# \begin{bmatrix}
# I_{11} & I_{12} \\
# I_{21} & I_{22} \\
# \end{bmatrix} =
# \begin{bmatrix}
# E[U_1U_1] & E[U_1U_2] \\
# E[U_2U_1] & E[U_2U_2] \\
# \end{bmatrix}
# $$
#
# と書く. このとき
#
# $$
# I = 
# \begin{bmatrix}
# I_{11} & I_{12} \\
# I_{21} & I_{22} \\
# \end{bmatrix} =
# \begin{bmatrix}
# E[-(\log L)_{\Delta\Delta}] & E[-(\log L)_{\Delta q}] \\
# E[-(\log L)_{q\Delta}] & E[-(\log L)_{q q}] \\
# \end{bmatrix}.
# $$
#
# 具体的には, $m/p+m/(1-p)=m/(p(1-p))$ などより, 
#
# $$
# \begin{aligned}
# &
# I_{11} = I_{12} = I_{21} = E\left[\frac{a}{p^2}+\frac{b}{(1-p)^2}\right] =
# \frac{m}{p(1-p)},
# \\ &
# I_{22} = E\left[\frac{a}{p^2}+\frac{b}{(1-p)^2}+\frac{c}{q^2}+\frac{d}{(1-q)^2}\right] =
# \frac{m}{p(1-p)}+\frac{n}{q(1-q)}.
# \end{aligned}
# $$
#
# $I_{ij}$ を $\Delta, q$ の函数とみなすときには $I_{ij}(\Delta, q)$ と書く.

# %% [markdown]
# $a,b,c,d$ がモデル内確率変数のとき, $m,n$ が十分大きいならば二項分布に関する中心極限定理より, スコア統計量 $U_1, U_2$ の同時分布は期待値がゼロで分散共分散行列が $I$ の多変量正規分布に近似的に従う. 以下ではこの条件を仮定する.

# %% [markdown]
# __最尤法:__ $\tq$ がパラメータ $\Delta$ の値が与えられているときの尤度 $L$ を最大化することは, 
#
# $$
# \tp = \tq + \Delta
# $$
#
# とおくとき, 
#
# $$
# U_2(\Delta, \tq) = \frac{a}{\tp}-\frac{b}{1-\tp}+\frac{c}{\tq}-\frac{d}{1-\tq} = 0
# $$
#
# と同値である. $U_2(\Delta, q)$ は $q$ について単調減少函数なのでその条件で $\tq$ は一意的に決まる.
#
# __注意:__ 厳密には同値ではないが, ここでは気にしないことにする. 実装時には問題になる.
#
# $U_1(\Delta, \tq)$ が従う分布は, 条件 $U_2=0$ の下での $U_1$ の条件付き確率分布になっていることに注意せよ.

# %% [markdown]
# __条件付き確率分布:__ 一般に期待値がゼロで分散共分散行列
#
# $$
# \Sigma = 
# \begin{bmatrix}
# \sigma_y^2  & \sigma_{yx} \\
# \sigma_{xy} & \sigma_x^2 \\
# \end{bmatrix}, 
# \quad
# \sigma_{yx}=\sigma_{xy}
# $$
#
# の多変量正規分布に従う確率辺数 $Y,X$ について, 条件 $X=x$ の下での $Y$ が従う条件付き確率分布は, 期待値が $\dfrac{\sigma_{xy}}{\sigma_x^2}x$ で分散が $\dfrac{\sigma_x^2\sigma_y^2-\sigma_{xy}^2}{\sigma_x^2}$ の正規分布になる.
#
# __注意:__ この結果は線形回帰の基礎でもある.
#
# __証明__: $\Sigma$ の逆行列は
#
# $$
# \Sigma^{-1} =
# \frac{1}{\sigma_x^2\sigma_y^2 - \sigma_{xy}^2}
# \begin{bmatrix}
# \sigma_x^2 & -\sigma_{xy} \\
# -\sigma_{xy} & \sigma_y^2 \\
# \end{bmatrix}
# $$
#
# なので
#
# $$
# \begin{aligned}
# \begin{bmatrix} y & x \end{bmatrix} \Sigma^{-1} \begin{bmatrix} y \\ x \end{bmatrix} &=
# \frac{1}{\sigma_x^2\sigma_y^2 - \sigma_{xy}^2}(\sigma_x^2 y^2 - 2\sigma_{xy}xy + \sigma_y^2 x^2)
# \\ &=
# \frac{1}{\sigma_x^2\sigma_y^2 - \sigma_{xy}^2}\left(
# \sigma_x^2\left(y - \frac{\sigma_{xy}}{\sigma_x^2}x\right)^2 + \frac{\sigma_x^2\sigma_y^2-\sigma_{xy}^2}{\sigma_x^2}x^2
# \right)
# \\ &=
# \frac{\sigma_x^2}{\sigma_x^2\sigma_y^2 - \sigma_{xy}^2}\left(y - \frac{\sigma_{xy}}{\sigma_x^2}x\right)^2 + \frac{x^2}{\sigma_x^2}.
# \end{aligned}
# $$
#
# このことから, $(Y,X)$ の同時確率密度函数は
#
# $$
# \begin{aligned}
# p(y,x) &= \text{const.} \exp\left(-\frac{1}{2}\begin{bmatrix} y & x \end{bmatrix} \Sigma^{-1} \begin{bmatrix} y \\ x \end{bmatrix}\right)
# \\ &=
# \text{const.}
# \exp\left(-\frac{1}{2}\frac{\sigma_x^2}{\sigma_x^2\sigma_y^2 - \sigma_{xy}^2}\left(y - \frac{\sigma_{xy}}{\sigma_x^2}x\right)^2\right)
# \exp\left(-\frac{x^2}{2\sigma_x^2}\right)
# \end{aligned}
# $$
#
# が $X$ の分布の期待値ゼロで分散 $\sigma_x^2$ なのでその密度函数は
#
# $$
# p(x)=\text{const.}\exp\left(-\frac{x^2}{2\sigma_x^2}\right)
# $$
#
# なので, $X=x$ という条件の下での $Y$ が従う条件付き確率分布の密度函数は
#
# $$
# p(y|x) = \frac{p(y,x)}{p(x)} =
# \text{const.}
# \exp\left(-\frac{1}{2}\frac{\sigma_x^2}{\sigma_x^2\sigma_y^2 - \sigma_{xy}^2}\left(y - \frac{\sigma_{xy}}{\sigma_x^2}x\right)^2\right)
# $$
#
# となる. ゆえに, $X=x$ という条件の下での $Y$ が従う条件付き確率分布は, 期待値が $\dfrac{\sigma_{xy}}{\sigma_x^2}x$ で分散が $\dfrac{\sigma_x^2\sigma_y^2-\sigma_{xy}^2}{\sigma_x^2}$ の正規分布になる. 
#
# __証明終__

# %% [markdown]
# __スコア法のP値:__ $Y=U_1$, $X=U_2$, $\Sigma=I$, $x=0$ の場合に上の結果を適用すると, 条件 $U_2=0$ の下での $U_1$ の条件付き確率分布すなわち $U_1(\Delta, \tq)$ の分布は, 期待値 $0$ で分散が
#
# $$
# \begin{aligned}
# V = V(\Delta, q) =
# \frac{I_{11}I_{22}-I_{12}^2}{I_22} &=
# \frac{\dfrac{m}{p(1-p)}\dfrac{n}{q(1-q)}}{\dfrac{m}{p(1-p)}+\dfrac{n}{q(1-q)}} =
# \left(\dfrac{p(1-p)}{m}+\dfrac{q(1-q)}{n}\right)^{-1}
# \end{aligned}
# $$
#
# の正規分布で近似されることがわかる. ここで $p=q+\Delta$ である. さらに $m,n$ が十分に大きいならば, 大数の法則よりこの $V$ は
#
# $$
# V(\Delta, \tq) = \left(\dfrac{\tp(1-\tp)}{m}+\dfrac{\tq(1-\tq)}{n}\right)^{-1}
# $$
#
# で近似される. このとき,
#
# $$
# Z := U_1(\Delta, \tq)\sqrt{V(\Delta, \tq)^{-1}}
# $$
#
# は標準正規分布に近似的に従う. この結果を使えば, $a,b,c,d$ が観測値の場合の $Z$ の値を $z$ と書いたときの「標準正規分布に従う確率変数の絶対値が $|z|$ 以上になる確率」としてP値が定義される:
#
# $$
# \pvalue(a,b,c,d|\Delta) = 2\op{ccdf}(\op{Normal}(0,1), |z|).
# $$
#
# ここで $\op{ccdf}(\op{Normal}(0,1), |z|)$ は標準正規分布 $\op{Normal}(0,1)$ において値が $|z|$ 以上になる確率を意味し, 2倍は両側の確率を足し合わせることを意味している. ($\op{ccdf}$ は補累積分布函数(complementary cumulative distribution function)を意味している.)
#
# $a,b,c,d$ が観測値の場合の $Z$ の値 $z$ はちょうど $a,b,c,d,\Delta$ のみから計算される値になっていることに注意せよ. $m=a+b$, $n=c+d$ であり, $\tq$ の値は条件 $U_2(\Delta, \tq)=0$ によって $a,b,c,d,\Delta$ から決まるのであった.
#
# 以上のように定義されるP値を __スコア法__ のP値と呼ぶ.

# %% [markdown]
# __Pearsonのχ²検定との関係:__ 上の代わりに, 
#
# $$
# Z^2 = U_1(\Delta, \tq)^2V(\Delta, \tq)^{-1}
# $$
#
# が自由度 $1$ の $\chi^2$ 分布に近似的に従うことを用いて
#
# $$
# \pvalue(a,b,c,d|\Delta) = \op{ccdf}(\op{Chisq}(1), z^2).
# $$
#
# によって同じP値を別の方法で定義することもできる.  ここで $\op{ccdf}(\op{Chisq}(1), z^2)$ は自由度 $1$ の $\chi^2$ 分布　$\op{Chisq}(1)$ で値が $z^2$ 以上になる確率を意味している.  ($\op{ccdf}$ は補累積分布函数(complementary cumulative distribution function)を意味している.)
#
# $\Delta=0$ のとき, $Z^2$ は2×2の分割表のPearsonのχ²統計量に一致するので, このP値は2×2の分割表の独立性のPearsonのχ²検定のP値にぴったり一致する.
#
# $\Delta=0$ のとき, $Z^2$ が2×2の分割表のPearsonのχ²統計量に一致することを証明しよう.
#
# __証明:__ $\Delta=0$ のとき,
#
# $$
# 0 = U_2(\Delta=0, \tq) = \frac{a}{\tq}-\frac{b}{1-\tq}+\frac{c}{\tq}-\frac{d}{1-\tq} =
# \frac{a+c}{\tq} - \frac{b+d}{1-\tq}
# $$
#
# より $\tq = \dfrac{a+c}{a+b+c+d}$ なので,
#
# $$
# \begin{aligned}
# &
# U_1(\Delta=0, \tq) = \frac{a}{\tq} - \frac{b}{1-\tq} =
# \frac{(a+b+c+d)(ad-bc)}{(a+c)(b+d)},
# \\ &
# V(\Delta=0, \tq)^{-1} = \frac{\tq(1-\tq)}{a+b} + \frac{\tq(1-\tq)}{c+d} =
# \frac{(a+c)(b+d)}{(a+b+c+d)(a+b)(c+d)}.
# \end{aligned}
# $$
#
# ゆえに, $\Delta=0$ のとき, 
#
# $$
# Z^2 = U_1(\Delta=0, \tq)^2 V(\Delta=0, \tq)^{-1} = 
# \frac{(a+b+c+d)(ad-bc)^2}{(a+b)(c+d)(a+c)(b+d)}.
# $$
#
# これは2×2の分割表のPearsonのχ²統計量の有名な公式に一致する.
#
# __証明終__

# %% [markdown]
# __信頼区間:__ 信頼区間はP値が $\alpha$ 以上になるパラメータの値の集合として定義される.  以上の文脈では母比率の差のモデル化になっているモデルのパラメータ $\Delta$ の信頼区間が次にように定義される:
#
# $$
# \confint_\Delta(a,b,c,d|\alpha) = \{\,\Delta\in[-1,1]\mid \pvalue(a,b,c,d|\Delta) \ge \alpha\,\}.
# $$

# %% [markdown]
# __実装例:__ 以下では以上を実際に実装してみる.

# %%
using Distributions
using Roots
using StatsPlots
default(fmt=:png)

myecdf(A, x) = count(≤(x), A)/length(A)
safemul(x, y) = x == 0 ? zero(x*y) : y == 0 ? zero(x*y) : x*y
safediv(x, y) = x == 0 ? zero(x/y) : isinf(y) ? zero(x/y) : x/y

# %%
### score method for risk difference

riskdiffhat_score(a, b, c, d) = safediv(a, a+b) - safediv(c, c+d)

function loglik_rd(a, b, c, d, q, Δ=0.0)
    p = q + Δ
    safemul(a, log(p)) + safemul(b, log(1-p)) + safemul(c, log(q)) + safemul(d, log(1-q))
end

function scorestat_q_rd(a, b, c, d, q, Δ=0.0)
    p = q + Δ
    safediv(a, p) - safediv(b, 1-p) + safediv(c, q) - safediv(d, 1-q)
end

function d_scorestat_q_rd(a, b, c, d, q, Δ=0.0)
    p = q + Δ
    -safediv(a, p^2) - safediv(b, (1-p)^2) - safediv(c, q^2) - safediv(d, (1-q)^2)
end

function scorestat_Δ_rd(a, b, c, d, q, Δ=0.0)
    p = q + Δ
    safediv(a, p) - safediv(b, 1-p)
end

function estimate_q_given_Δ_rd(a, b, c, d, Δ=0.0; alg=Bisection())
    qmin, qmax = max(0.0, -Δ), min(1.0, 1.0-Δ)
    a+c==0 && return qmin
    b+d==0 && return qmax
    f(q) = scorestat_q_rd(a, b, c, d, q, Δ)
    S_qmin = f(qmin + eps())
    S_qmax = f(qmax - eps())
    S_qmin ≥ 0 && S_qmax ≥ 0 && return S_qmin < S_qmax ? qmin : qmax
    S_qmin ≤ 0 && S_qmax ≤ 0 && return S_qmin < S_qmax ? qmax : qmin
    find_zero(f, (qmin + eps(), qmax - eps()), alg)
end

function varinv_scorestat_q_rd(a, b, c, d, q, Δ=0.0)
    p = q + Δ
    safediv(p*(1-p), a+b) + safediv(q*(1-q), c+d)
end

function chisqstat_rd_score(a, b, c, d; Δ=0.0, alg=Bisection())
    Δ = clamp(Δ, -1 + eps(), 1 - eps())
    q̃ = estimate_q_given_Δ_rd(a, b, c, d, Δ; alg)
    S = scorestat_Δ_rd(a, b, c, d, q̃, Δ)
    Vinv = varinv_scorestat_q_rd(a, b, c, d, q̃, Δ)
    safemul(S^2, Vinv)
end

function pvalue_rd_score(a, b, c, d; Δ=0.0, alg=Bisection())
    χ² = chisqstat_rd_score(a, b, c, d; Δ, alg)
    ccdf(Chisq(1), χ²)
end

function confint_rd_score(a, b, c, d; α=0.05, alg=Bisection())
    χ²_α = cquantile(Chisq(1), α)
    RDhat = riskdiffhat_score(a, b, c, d)
    g(Δ) = chisqstat_rd_score(a, b, c, d; Δ, alg) - χ²_α
    L = if g(-1 + eps()) > 0
        find_zero(g, (-1 + eps(), RDhat), alg)
    else
        -1.0
    end
    U = if g(1 - eps()) > 0
        find_zero(g, (RDhat, 1 - eps()), alg)
    else
        1.0
    end
    [L, U]
end

# %%
a, b, c, d = 8, 0, 8, 0
@show Δ0 = riskdiffhat_score(a, b, c, d)
@show chisqstat_rd_score(a, b, c, d; Δ=Δ0)
chisqstat_rd_score(a, b, c, d; Δ=-0.9999)

# %%
a, b, c, d = 1, 2, 3, 4
pvalue_rd_score(1, 2, 3, 4), ccdf(Chisq(1), (a+b+c+d)*(a*d-b*c)^2/((a+b)*(c+d)*(a+c)*(b+d)))

# %%
a, b, c, d = 3, 0, 4, 0
pvalue_rd_score(a, b, c, d), ccdf(Chisq(1), safediv((a+b+c+d)*(a*d-b*c)^2, (a+b)*(c+d)*(a+c)*(b+d)))

# %%
riskdiffhat(a, b, c, d) = safediv(a, a+b) - safediv(c, c+d)

function stderr_riskdiffhat(a, b, c, d)
    m, n = a+b, c+d
    p̂, q̂ = safediv(a, m), safediv(c, n)
    √(safediv(p̂*(1-p̂), m) + safediv(q̂*(1-q̂), n))
end

function pvalue_rd_wald(a, b, c, d; Δ=0)
    RDhat = riskdiffhat(a, b, c, d)
    SEhat_riskdiffhat = stderr_riskdiffhat(a, b, c, d)
    2ccdf(Normal(0, 1), safediv(abs(RDhat - Δ), SEhat_riskdiffhat))
end

function confint_rd_wald(a, b, c, d; α=0.05)
    z = quantile(Normal(), 1-α/2)
    RDhat = riskdiffhat(a, b, c, d)
    SEhat_riskdiffhat = stderr_riskdiffhat(a, b, c, d)
    [RDhat - z*SEhat_riskdiffhat, RDhat + z*SEhat_riskdiffhat]
end

riskdiffhat_zou_donner(a, b, c, d) = safediv(a, a+b) - safediv(c, c+d)

function stderr_riskdiffhat_zou_donner(a, b, c, d; u=1)
    m, n = a+b, c+d
    p̂, q̂ = safediv(a, m), safediv(c, n)
    √(safediv(p̂*(1-p̂), m-u) + safediv(q̂*(1-q̂), n-u))
end

function pvalue_rd_zou_donner(a, b, c, d; Δ=0, u=1)
    ((a==0 && d==0) || (b==0 && c==0)) && return 1.0
    RDhat = riskdiffhat_zou_donner(a, b, c, d)
    SEhat_riskdiffhat = stderr_riskdiffhat_zou_donner(a, b, c, d; u)
    Z = safediv((1 - RDhat^2)*abs(atanh(RDhat) - atanh(Δ)), SEhat_riskdiffhat)
    2ccdf(Normal(), abs(Z))
end

function confint_rd_zou_donner(a, b, c, d; α=0.05, u=1)
    z = quantile(Normal(), 1-α/2)
    RDhat = riskdiffhat_zou_donner(a, b, c, d)
    SEhat_riskdiffhat = stderr_riskdiffhat_zou_donner(a, b, c, d; u)
    m = atanh(RDhat)
    d = safediv(z*SEhat_riskdiffhat, 1 - RDhat^2)
    [tanh(m-d), tanh(m+d)]
end

a, b, c, d  = 58, 22, 62, 38

@show riskdiffhat(a, b, c, d)
@show stderr_riskdiffhat(a, b, c, d)
@show pvalue_rd_wald(a, b, c, d)
@show pvalue_rd_zou_donner(a, b, c, d)
@show pvalue_rd_score(a, b, c, d)
@show pvalue_rd_wald(a, b, c, d; Δ=0.2)
@show pvalue_rd_zou_donner(a, b, c, d; Δ=0.2)
@show pvalue_rd_score(a, b, c, d; Δ=0.2)
@show confint_rd_wald(a, b, c, d; α=0.05)
@show confint_rd_zou_donner(a, b, c, d; α=0.05)
@show confint_rd_score(a, b, c, d; α=0.05)

plot(Δ -> pvalue_rd_wald(a, b, c, d; Δ), -0.2, 0.4; label="")
plot!(Δ -> pvalue_rd_zou_donner(a, b, c, d; Δ), -0.2, 0.4; label="", ls=:dash)
plot!(Δ -> pvalue_rd_score(a, b, c, d; Δ), -0.2, 0.4; label="", ls=:dashdot)
plot!(confint_rd_wald(a, b, c, d; α=0.05), fill(0.05, 2); label="95% CI")
plot!(confint_rd_wald(a, b, c, d; α=0.20), fill(0.20, 2); label="80% CI")
plot!(xguide="ratio difference", yguide="P-value")
plot!(xtick=-2:0.05:2, ytick=0:0.05:1)

# %%
a, b = 8, 2
c, d = 3, 7

@show riskdiffhat(a, b, c, d)
@show stderr_riskdiffhat(a, b, c, d)
@show pvalue_rd_wald(a, b, c, d; Δ=-1.0)
@show pvalue_rd_zou_donner(a, b, c, d; Δ=-1.0)
@show pvalue_rd_score(a, b, c, d; Δ=-1.0)
@show pvalue_rd_wald(a, b, c, d)
@show pvalue_rd_zou_donner(a, b, c, d)
@show pvalue_rd_score(a, b, c, d)
@show pvalue_rd_wald(a, b, c, d; Δ=0.2)
@show pvalue_rd_zou_donner(a, b, c, d; Δ=0.2)
@show pvalue_rd_score(a, b, c, d; Δ=0.2)
@show pvalue_rd_wald(a, b, c, d; Δ=1.0)
@show pvalue_rd_zou_donner(a, b, c, d; Δ=1.0)
@show pvalue_rd_score(a, b, c, d; Δ=1.0)
@show confint_rd_wald(a, b, c, d; α=0.05)
@show confint_rd_zou_donner(a, b, c, d; α=0.05)
@show confint_rd_score(a, b, c, d; α=0.05)

plot(Δ -> pvalue_rd_wald(a, b, c, d; Δ), -1, 1.0; label="Wald")
plot!(Δ -> pvalue_rd_zou_donner(a, b, c, d; Δ); label="Zou-Donner", ls=:dash)
plot!(Δ -> pvalue_rd_score(a, b, c, d; Δ); label="score", ls=:dashdot)
plot!(confint_rd_wald(a, b, c, d; α=0.05), fill(0.05, 2); label="95% CI")
plot!(confint_rd_wald(a, b, c, d; α=0.20), fill(0.20, 2); label="80% CI")
plot!(xguide="ratio difference", yguide="P-value")
plot!(xtick=-2:0.1:2, ytick=0:0.05:1)

# %%
a, b = 190, 10
c, d = 180, 20

@show riskdiffhat(a, b, c, d)
@show stderr_riskdiffhat(a, b, c, d)
@show pvalue_rd_wald(a, b, c, d; Δ=-1.0)
@show pvalue_rd_zou_donner(a, b, c, d; Δ=-1.0)
@show pvalue_rd_score(a, b, c, d; Δ=-1.0)
@show pvalue_rd_wald(a, b, c, d)
@show pvalue_rd_zou_donner(a, b, c, d)
@show pvalue_rd_score(a, b, c, d)
@show pvalue_rd_wald(a, b, c, d; Δ=0.2)
@show pvalue_rd_zou_donner(a, b, c, d; Δ=0.2)
@show pvalue_rd_score(a, b, c, d; Δ=0.2)
@show pvalue_rd_wald(a, b, c, d; Δ=1.0)
@show pvalue_rd_zou_donner(a, b, c, d; Δ=1.0)
@show pvalue_rd_score(a, b, c, d; Δ=1.0)
@show confint_rd_wald(a, b, c, d; α=0.05)
@show confint_rd_zou_donner(a, b, c, d; α=0.05)
@show confint_rd_score(a, b, c, d; α=0.05)

plot(Δ -> pvalue_rd_wald(a, b, c, d; Δ), -0.5, 1.0; label="Wald")
plot!(Δ -> pvalue_rd_zou_donner(a, b, c, d; Δ); label="Zou-Donner", ls=:dash)
plot!(Δ -> pvalue_rd_score(a, b, c, d; Δ); label="score", ls=:dashdot)
plot!(confint_rd_wald(a, b, c, d; α=0.05), fill(0.05, 2); label="95% CI")
plot!(confint_rd_wald(a, b, c, d; α=0.20), fill(0.20, 2); label="80% CI")
plot!(xguide="ratio difference", yguide="P-value")
plot!(xtick=-2:0.1:2, ytick=0:0.05:1)

# %%
a, b = 8, 0
c, d = 9, 0

@show riskdiffhat(a, b, c, d)
@show stderr_riskdiffhat(a, b, c, d)
@show pvalue_rd_wald(a, b, c, d; Δ=-1.0)
@show pvalue_rd_zou_donner(a, b, c, d; Δ=-1.0)
@show pvalue_rd_score(a, b, c, d; Δ=-1.0)
@show pvalue_rd_wald(a, b, c, d)
@show pvalue_rd_zou_donner(a, b, c, d)
@show pvalue_rd_score(a, b, c, d)
@show pvalue_rd_wald(a, b, c, d; Δ=0.2)
@show pvalue_rd_zou_donner(a, b, c, d; Δ=0.2)
@show pvalue_rd_score(a, b, c, d; Δ=0.2)
@show pvalue_rd_wald(a, b, c, d; Δ=1.0-eps())
@show pvalue_rd_zou_donner(a, b, c, d; Δ=1.0-eps())
@show pvalue_rd_score(a, b, c, d; Δ=1.0-eps())
@show pvalue_rd_wald(a, b, c, d; Δ=1.0)
@show pvalue_rd_zou_donner(a, b, c, d; Δ=1.0)
@show pvalue_rd_score(a, b, c, d; Δ=1.0)
@show confint_rd_wald(a, b, c, d; α=0.05)
@show confint_rd_zou_donner(a, b, c, d; α=0.05)
@show confint_rd_score(a, b, c, d; α=0.05)

plot(Δ -> pvalue_rd_wald(a, b, c, d; Δ), -1.0, 1.0; label="Wald")
plot!(Δ -> pvalue_rd_zou_donner(a, b, c, d; Δ); label="Zou-Donner", ls=:dash)
plot!(Δ -> pvalue_rd_score(a, b, c, d; Δ); label="score", ls=:dashdot)
plot!(confint_rd_wald(a, b, c, d; α=0.05), fill(0.05, 2); label="95% CI")
plot!(confint_rd_wald(a, b, c, d; α=0.20), fill(0.20, 2); label="80% CI")
plot!(xguide="ratio difference", yguide="P-value")
plot!(xtick=-2:0.1:2, ytick=0:0.05:1)

# %%
a, b = 190, 10
c, d = 180, 20
Δ = 0.5

@show estimate_q_given_Δ_rd(a, b, c, d, Δ)
@show scorestat_q_rd(a, b, c, d, 0.0, Δ)
@show scorestat_q_rd(a, b, c, d, 0.675, Δ)
@show scorestat_q_rd(a, b, c, d, 0.7, Δ)

find_zero(q -> (q; scorestat_q_rd(a, b, c, d, q, Δ)), (0.0, 1-Δ))

# %%
a, b = 190, 10
c, d = 180, 20
Δ = -1.0

@show scorestat_q_rd(a, b, c, d, 0.0, 0.3)
@show scorestat_q_rd(a, b, c, d, 0.0, 0.3)

qmin, qmax = max(0.0, -Δ), min(1.0, 1.0-Δ)
a+c==0 && @show qmin
b+d==0 && @show qmax
@show S_qmin = scorestat_q_rd(a, b, c, d, qmin+eps(), Δ)
@show S_qmax = scorestat_q_rd(a, b, c, d, qmax-eps(), Δ)
S_qmin ≥ 0 && S_qmax ≥ 0 && @show S_qmin < S_qmax ? qmin : qmax
S_qmin ≤ 0 && S_qmax ≤ 0 && @show S_qmin < S_qmax ? qmax : qmin
f(q) = scorestat_q_rd(a, b, c, d, q, Δ)
@show q0 = (qmin + qmax)/2
@show find_zero(f, q0, Order0())
@show find_zero(f, (qmin+eps(), qmax-eps()))


plot(q -> scorestat_q_rd(a, b, c, d, q, 0.3), 0.01, 0.69)

# %%
a, b = 10, 0
c, d = 20, 0

@show scorestat_q_rd(a, b, c, d, 0.0, 0.3)
@show scorestat_q_rd(a, b, c, d, 0.7, 0.3)

plot(q -> scorestat_q_rd(a, b, c, d, q, 0.3), 0.1, 0.69)

# %%
a, b = 0, 2
c, d = 0, 7

@show riskdiffhat(a, b, c, d)
@show stderr_riskdiffhat(a, b, c, d)
@show pvalue_rd_wald(a, b, c, d)
@show pvalue_rd_zou_donner(a, b, c, d)
@show pvalue_rd_score(a, b, c, d)
@show pvalue_rd_wald(a, b, c, d; Δ=0.2)
@show pvalue_rd_zou_donner(a, b, c, d; Δ=0.2)
@show pvalue_rd_score(a, b, c, d; Δ=0.2)

plot(Δ -> pvalue_rd_score(a, b, c, d; Δ), -1, 1; label="score", ls=:dashdot)
plot!(confint_rd_wald(a, b, c, d; α=0.05), fill(0.05, 2); label="95% CI")
plot!(confint_rd_wald(a, b, c, d; α=0.20), fill(0.20, 2); label="80% CI")
plot!(xguide="ratio difference", yguide="P-value")
plot!(xtick=-1:0.2:1, ytick=0:0.05:1)

# %%
a, b = 4, 0
c, d = 10, 1

@show riskdiffhat(a, b, c, d)
@show stderr_riskdiffhat(a, b, c, d)
@show pvalue_rd_wald(a, b, c, d)
@show pvalue_rd_zou_donner(a, b, c, d)
@show pvalue_rd_score(a, b, c, d)
@show pvalue_rd_wald(a, b, c, d; Δ=0.2)
@show pvalue_rd_zou_donner(a, b, c, d; Δ=0.2)
@show pvalue_rd_score(a, b, c, d; Δ=0.2)
@show pvalue_rd_wald(a, b, c, d; Δ=-1.0)
@show pvalue_rd_zou_donner(a, b, c, d; Δ=-1.0)
@show pvalue_rd_score(a, b, c, d; Δ=-1.0+√eps())
@show pvalue_rd_wald(a, b, c, d; Δ=1.0)
@show pvalue_rd_zou_donner(a, b, c, d; Δ=1.0)
@show pvalue_rd_score(a, b, c, d; Δ=1.0)

plot(Δ -> pvalue_rd_score(a, b, c, d; Δ), -1, 1; label="score", ls=:dashdot)
plot!(xguide="ratio difference", yguide="P-value")
plot!(xtick=-1:0.2:1, ytick=0:0.05:1)

# %%

# %%

# %%

# %%
a, b = 8, 0
c, d = 8, 1

@show scorestat_q_rd(a, b, c, d, 0.0, 0.2)
@show scorestat_q_rd(a, b, c, d, 0.8, 0.2)

plot(q -> scorestat_q_rd(a, b, c, d, q, 0.2), 0.2, 0.79)

# %%

# %%

# %%

# %%
function plot_ecdf_pvals(; m=10, n=10, p=0.8, q=0.3, Δ=p-q, L=10^5)
    a = rand(Binomial(m, p), L)
    b = m .- a
    c = rand(Binomial(n, q), L)
    d = n .- c
    pval_wald = pvalue_rd_wald.(a, b, c, d; Δ)
    pval_zd = pvalue_rd_zou_donner.(a, b, c, d; Δ)
    pval_score = pvalue_rd_score.(a, b, c, d; Δ)

    P = plot(α -> myecdf(pval_wald, α), 0, 0.1; label="Wald")
    plot!(α -> myecdf(pval_zd, α); label="Zou-Donner")
    plot!(α -> myecdf(pval_score, α); label="socre")
    plot!(identity; label="", c=:gray, ls=:dot)
    plot!(xtick=0:0.01:0.1)
    Δ == p - q && plot!(ytick=0:0.01:1)
    
    Q = plot(α -> myecdf(pval_wald, α), 0, 1; label="Wald")
    plot!(α -> myecdf(pval_zd, α); label="Zou-Donner")
    plot!(α -> myecdf(pval_score, α); label="socre")
    plot!(identity; label="", c=:gray, ls=:dot)
    plot!(xtick=0:0.1:1, ytick=0:0.1:1)

    plot(P, Q; size=(800, 420), layout=(1, 2))
    plot!(plot_title="m=$m, n=$n, p=$p, q=$q, Δ=$Δ", plot_titlefontsize=10)
end

# %%
plot_ecdf_pvals(; m=10, n=10, p=0.8, q=0.3)

# %%
plot_ecdf_pvals(; m=20, n=20, p=0.8, q=0.3)

# %%
plot_ecdf_pvals(; m=40, n=40, p=0.8, q=0.3)

# %%
plot_ecdf_pvals(; m=20, n=100, p=0.8, q=0.3)

# %%
plot_ecdf_pvals(; m=100, n=20, p=0.8, q=0.3)

# %%
plot_ecdf_pvals(; m=20, n=100, p=0.5, q=0.5)

# %%
plot_ecdf_pvals(; m=100, n=20, p=0.5, q=0.5)

# %%
plot_ecdf_pvals(; m=60, n=60, p=0.5, q=0.5)

# %%
plot_ecdf_pvals(; m=60, n=60, p=0.3, q=0.5)

# %%
a, b, c, d  = 58, 22, 62, 38
plot_ecdf_pvals(; m=a+b, n=c+d, p=a/(a+b), q=c/(c+d))

# %%
a, b, c, d  = 58, 22, 62, 38
plot_ecdf_pvals(; m=a+b, n=c+d, p=a/(a+b), q=c/(c+d), Δ=0.0)

# %%
a, b, c, d  = 58, 22, 62, 38
plot_ecdf_pvals(; m=a+b, n=c+d, p=a/(a+b), q=c/(c+d), Δ=0.2)

# %%

# %%

# %%

# %%
function polynom_scorestat_q_rd(a, b, c, d, q, Δ=0.0)
    p = q + Δ
    a*(1-p)*q*(1-q) - b*p*q*(1-q) + c*p*(1-p)*(1-q) - d*p*(1-p)*q
end

# %%
a, b, c, d = 8, 1, 8, 1
@show Δ = 0.2
@show riskdiffhat_score(a, b, c, d)
P = plot(q -> polynom_scorestat_q_rd(a, b, c, d, q, Δ), 0, 1-Δ)
Q = plot(q -> loglik_rd(a, b, c, d, q, Δ), 0, 1-Δ)
plot(P, Q; layout=(1, 2), size=(800, 250))

# %%
plot(q -> loglik_rd(a, b, c, d, q, Δ), 1-Δ-0.1, 1-Δ)

# %%

# %%
