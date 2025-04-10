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
# # P値函数のプロットの例
#
# * 黒木玄
# * 2023-12-02～2023-12-03
# $
# \newcommand\op{\operatorname}
# $
#
# データの数値は「$n$ 回中 $k$ 回表が出た」のスタイルで, データの生成のされ方に関する統計モデルとして, 表の出る確率が $p$ の二項分布を考える. すなわち $k$ に対応する確率変数 $K$ について $K\sim\op{Binomial}(n, p)$ が成立しているというモデルを考える.
#
# 中心極限定理によって, モデルの確率変数 $K$ について $\dfrac{K/n - p}{\sqrt{p(1-p)/n}}$ は近似的に標準正規分布に従う.
#
# このことを使うと, データの数値「$n$ 回中 $k$ 回」に関する仮説「表の出る確率は $p$ である」のP値 $\op{pvalue}_{\op{Wilson}}(k|n,p)$ を次のように定めることができる:
#
# $$
# \op{pvalue}_{\op{Wilson}}(k|n,p) = 2(1 - \op{cdf}(\op{Normal}(0, 1), |z|)), \quad
# z = \frac{k/n - p}{\sqrt{p(1-p)/n}}.
# $$
#
# ここで $\op{cdf}(\op{Normal}(0, 1), x)$ は標準正規分布 $\op{Normal}(0,1)$ の累積分布函数 (cumulative distribution function)である.
#
# 上の定義によって得られる函数 $p\mapsto \op{pvalue}(k|n,p)$ を __WilsonのP値函数__ と呼ぶことにする.
#
# $0\le\alpha\le 1$ であると仮定する.
#
# 上のP値の定義に対応する表の出る確率 $p$ の $100(1-\alpha)\%$ 信頼区間 $\op{confint}_{\op{Wilson}}(k|n,\alpha)$ は次のように定義される:
#
# $$
# \op{confint}_{\op{Wilson}}(k|n,\alpha) = \{\, p\in [0,1] \mid \op{pvalue}_{\op{Wilson}}(k|n,p) \ge \alpha\,\}.
# $$
#
# これを __Wilsonの信頼区間__ と呼ぶ.
#
# 以下では, 標準正規分布 $\op{Normal}(0,1)$ の分位点函数(quantile function, 累積分布函数の逆函数)を $\op{quantile}(\op{Normal}(0,1), p)$ と書くことにする.
#
# __練習問題1 (Wilsonの信頼区間の具体的な表示):__ 以上の状況の下で
#
# $$
# \begin{aligned}
# &
# \hat{p} = k/n, \quad z_{\alpha/2} = \op{quantile}(\op{Normal}(0,1), 1-\alpha/2)),
# \\ &
# a = 1+\frac{z_{\alpha/2}^2}{n}, \quad b = \hat{p}+\frac{z_{\alpha/2}^2}{2n}, \quad c = \hat{p}^2
# \end{aligned}
# $$
#
# とおき, $p$ に関する二次方程式 $ap^2-2bp+c=0$ の2つの解を小さい順に $L$, $U$ と書くとき, 
#
# $$
# \op{confint}_{\op{Wilson}}(k|n,\alpha) = [L, U]
# $$
#
# となることを示せ. (注意: $0\le L \le U \le 1$ も示す必要がある.  この問題の解答例はこのノートの終わりの方にある.)
#
# Wilsonの信頼区間については次のリンク先も参照:
#
# * [ウィルソンの信頼区間](https://ja.wikipedia.org/wiki/%E3%82%A6%E3%82%A3%E3%83%AB%E3%82%BD%E3%83%B3%E3%81%AE%E4%BF%A1%E9%A0%BC%E5%8C%BA%E9%96%93)
# * https://scholar.google.co.jp/scholar?cluster=5129299358902170657
# * https://x.com/genkuroki/status/1730818456881377514
# * https://x.com/genkuroki/status/1730917888184529149
#
# このノートではWilsonのP値函数のグラフをプロットし, そこに信頼区間や点推定値も描き込む.

# %%
using Printf
using Distributions
using StatsPlots
default(fmt=:png, titlefontsize=10)
safediv(x, y) = x == 0 ? zero(x/y) : x/y

function pvalue_binomial_wilson(k, n, p)
    phat, se = k/n, √(p * (1 - p) / n)
    z = safediv(phat - p, se)
    2ccdf(Normal(), abs(z))
end

function confint_binomia_wilson(k, n, α=0.05)
    phat, z = k/n, quantile(Normal(), 1-α/2)
    a, b, c = 1+z^2/n, phat+z^2/(2n), phat^2
    sqrtD = √(b^2 - a*c)
    [(b - sqrtD)/a, (b + sqrtD)/a]
end

r2(x) = @sprintf "%.2g" x

function plot_binomial_test(; k=19, n=50, α=0.05, p0=0.5, p1=0.6, plotest=true)
    plot(p -> pvalue_binomial_wilson(k, n, p), 0, 1; label="P-value function")
    c = 2
    p = p0
    Pval = pvalue_binomial_wilson(k, n, p)
    vline!([p]; label="", ls=:dot, c)
    scatter!([p], [Pval]; label="P-value($p) = $(r2(Pval))", msc=:auto, c)
    c = 3
    p = p1
    Pval = pvalue_binomial_wilson(k, n, p)
    vline!([p]; label="", ls=:dot, c)
    scatter!([p], [Pval]; label="P-value($p) = $(r2(Pval))", msc=:auto, c)
    plotest && begin
        c = :blue
        phat = k/n
        vline!([phat]; label="point estimate = $(r2(phat))", msc=:auto, c)
        c = :red
        hline!([α]; label="α = $(100α)%", ls=:dot, c)
        L, U = confint_binomia_wilson(k, n, α)
        plot!([L, U], fill(α, 2); label="$(100(1-α))% CI = [$(r2(L)), $(r2(U))]", lw=3, c)
    end
    plot!(xtick=0:0.05:1, ytick=0:0.05:1, xrotation=45)
    plot!(xguide="p", yguide="P-value")
    title!("P-value function for data k=$k, n=$n
        (model: normal approximation of binomial distribution)")
end

# %%
plot_binomial_test(; k=19, n=50, α=0.05, p0=0.5, p1=0.6)

# %%
plot_binomial_test(; k=19, n=50, α=0.05, p0=0.5, p1=0.6, plotest=false)

# %%
plot_binomial_test(; k=113, n=200, α=0.05, p0=0.5, p1=0.6)

# %%
plot_binomial_test(; k=113, n=200, α=0.05, p0=0.5, p1=0.6, plotest=false)

# %%
plot_binomial_test(; k=268, n=500, α=0.05, p0=0.5, p1=0.6)

# %%
plot_binomial_test(; k=268, n=500, α=0.05, p0=0.5, p1=0.6, plotest=false)

# %% [markdown]
# __練習問題1解答例:__ 分位点函数は累積分布函数の逆函数なので, $z_{\alpha/2} = \op{quantile}(\op{Normal}(0,1), 1-\alpha/2))$ は $\op{cdf}(\op{Normal}(0,1), z_{\alpha/2}) = 1 - \alpha/2$ と同値であり, その条件はさらに次と同値である:
#
# $$
# 2(1 - \op{cdf}(\op{Normal}(0, 1), z_{\alpha/2})) = \alpha.
# $$
#
# $\hat{p}=k/n$ とおき, $p\in[0,1]$ と仮定する.
#
# $2(1 - \op{cdf}(\op{Normal}(0, 1), x))$ は実数 $x$ の狭義単調減少函数なので, 
#
# $$
# \op{pvalue}_{\op{Wilson}}(k|n,p) 
# = 2\left(1 - \op{cdf}\left(\op{Normal}(0, 1), \frac{|\hat{p} - p|}{\sqrt{p(1-p)/n}} \right)\right)
# \ge \alpha
# $$
#
# と
#
# $$
# \frac{|\hat{p} - p|}{\sqrt{p(1-p)/n}} \le z_{\alpha/2}
# $$
#
# は同値であり, さらに次と同値になる:
#
# $$
# \frac{(\hat{p}-p)^2}{p(1-p)/n} \le z_{\alpha/2}^2.
# $$
#
# これは次と同値である:
#
# $$
# (\hat{p}-p)^2 \le \frac{z_{\alpha/2}^2}{n} p(1-p).
# $$
#
# さらにこれは次と同値である:
#
# $$
# \left(1 + \frac{z_{\alpha/2}^2}{n}\right)p^2 - 2\left(\hat{p} + \frac{z_{\alpha/2}^2}{2n}\right)p + \hat{p}^2 \le 0.
# $$
#
# 左辺は $p=0$ のとき $\hat{p}^2\ge 0$ となり, $p=1$ のとき $(1 - \hat{p})^2 \ge 0$ となるので, この不等式を満たす実数 $p$ 全体の集合は区間 $[0,1]$ に含まれる.  その両端の値は2次方程式
#
# $$
# ap^2-2bp+c=0, \quad
# a = 1+\frac{z_{\alpha/2}^2}{n}, \quad
# b = \hat{p}+\frac{z_{\alpha/2}^2}{2n}, \quad
# c = \hat{p}^2
# $$
#
# の解になっている.  これでWilsonの信頼区間
#
# $$
# \op{confint}_{\op{Wilson}}(k|n,\alpha) = \{\, p\in [0,1] \mid \op{pvalue}_{\op{Wilson}}(k|n,p) \ge \alpha\,\}.
# $$
#
# の両端の点が $p$ に関する2次方程式 $ap^2-2bp+c=0$ の2つの解に等しいことがわかった.

# %% [markdown]
# __定義:__ 以上の記号の下で __Waldの信頼区間__ が次のように定義される:
#
# $$
# \op{confint}_{\op{Wald}}(k|n,\alpha) = 
# \left[
# \hat{p} - z_{\alpha/2}\sqrt{\frac{\hat{p}(1-\hat{p})}{n}},\;
# \hat{p} + z_{\alpha/2}\sqrt{\frac{\hat{p}(1-\hat{p})}{n}}
# \right].
# $$
#
# Waldの信頼区間は区間 $[0, 1]$ をはみだすことがあり得ることに注意せよ. さらに, $\hat{p}=k/n$ が $0$ または $1$ の場合にWaldの信頼区間は定義されない. これらの性質はWaldの信頼区間の欠点である.
#
# __練習問題2 (Waldの信頼区間による近似):__ Waldの信頼区間はWilsonの信頼区間のどのような近似になっているか?

# %% [markdown]
# __解答例:__
#
# $$
# a = 1+\frac{z_{\alpha/2}^2}{n}, \quad
# b = \hat{p}+\frac{z_{\alpha/2}^2}{2n}, \quad
# c = \hat{p}^2
# $$
#
# より, 
#
# $$
# \begin{aligned}
# b^2-ac &=
# \hat{p}^2 + \frac{z_{\alpha/2}^2}{n} \hat{p} + \frac{z_{\alpha/2}^4}{4n^2}
# - \hat{p}^2 - \frac{z_{\alpha/2}^2}{n} \hat{p}^2
# \\ &=
# \frac{z_{\alpha/2}^2}{n} \hat{p}(1-\hat{p}) + \frac{z_{\alpha/2}^4}{4n^2}
# \\ &=
# z_{\alpha/2}^2 \frac{\hat{p}(1-\hat{p})}{n}
# \left(1 + \frac{z_{\alpha/2}^2}{4n\hat{p}(1-\hat{p})}\right).
# \end{aligned}
# $$
#
# これより, Wilsonの信頼区間の両端の値(それらは練習問題1の結果より $p$ に関する2次方程式 $ap^2-2bp+c=0$ の2つの解に等しい)は次のように表される:
#
# $$
# \frac{b\pm\sqrt{b^2-ac}}{a} =
# \frac
# {
# \hat{p}+\dfrac{z_{\alpha/2}^2}{2n}
# \pm z_{\alpha/2} \sqrt{\dfrac{\hat{p}(1-\hat{p})}{n}}
# \sqrt{1 + \dfrac{z_{\alpha/2}^2}{4n\hat{p}(1-\hat{p})}}
# }
# {1+\dfrac{z_{\alpha/2}^2}{n}}.
# $$
#
# ゆえに, $0<\hat{p}<1$ のとき, $n$ を十分大きくすると, $\dfrac{z_{\alpha/2}^2}{4n\hat{p}(1-\hat{p})}$ は十分小さくなり, そのとき $\dfrac{z_{\alpha/2}^2}{n}$, $\dfrac{z_{\alpha/2}^2}{2n}$ も十分小さくなるので, Wilsonの信頼区間の両端の値はWaldの信頼区間の両端の値
#
# $$
# \hat{p} \pm z_{\alpha/2}\sqrt{\frac{\hat{p}(1-\hat{p})}{n}}
# $$
#
# で近似される. $\hat{p}=k/n$ が $0$ または $1$ に近いとき, $\hat{p}(1-\hat{p})$ は小さな値になるので, その分だけ $n$ を大きくしないと, Waldの信頼区間によるWilsonの信頼区間の近似の精度が上がらないことに注意せよ.

# %% [markdown]
# __練習問題3 (Waldの信頼区間に対応するP値函数):__ Waldの信頼区間を与えるP値函数を構成せよ.

# %% [markdown]
# __解答例:__ 次のように __WaldのP値函数__ $\op{pvalue}_{\op{Wald}}(k|n,p)$ を定める:
#
# $$
# \op{pvalue}_{\op{Wald}}(k|n,p) = 2(1 - \op{cdf}(\op{Normal}(0, 1), |\tilde{z}|)), \quad
# \tilde{z} = \frac{\hat{p} - p}{\sqrt{\hat{p}(1-\hat{p})/n}}.
# $$
#
# ここで $\hat{p}=k/n$. WaldのP値函数の定義とWaldのP値函数の定義
#
# $$
# \op{pvalue}_{\op{Wilson}}(k|n,p) = 2(1 - \op{cdf}(\op{Normal}(0, 1), |z|)), \quad
# z = \frac{\hat{p} - p}{\sqrt{p(1-p)/n}}
# $$
#
# の違いは $z$ と $\tilde{z}$ の定義式の分母における $p$ と $\hat{p}$ の違いである.
#
# このとき, 練習問題1の解答例と同様に考えると, 
#
# $$
# \begin{aligned}
# &
# \op{pvalue}_{\op{Wald}}(k|n,p) \ge \alpha
# \\ &\iff
# |\tilde{z}| = \frac{|\hat{p} - p|}{\sqrt{\hat{p}(1-\hat{p})/n}} \le z_{\alpha/2}
# \\ &\iff
# \hat{p} - z_{\alpha/2}\sqrt{\frac{\hat{p}(1-\hat{p})}{n}}
# \le p \le
# \hat{p} + z_{\alpha/2}\sqrt{\frac{\hat{p}(1-\hat{p})}{n}}.
# \end{aligned}
# $$
#
# すなわち,
#
# $$
# \op{confint}_{\op{Wald}}(k|n,\alpha) = \{\, p \in \mathbb{R} \mid \op{pvalue}_{\op{Wald}}(k|n,p) \ge \alpha\,\}.
# $$
#
# これでWaldのP値函数がWaldの信頼区間を与えることがわかった.

# %%
