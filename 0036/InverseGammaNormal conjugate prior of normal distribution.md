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
    display_name: Julia 1.8.0
    language: julia
    name: julia-1.8
---

# 正規分布モデルの共役事前分布によるベイズ統計

* 黒木玄
* 2022-09-03
$
\newcommand\ds{\displaystyle}
\newcommand\op[1]{{\operatorname{#1}}}
\newcommand\R{{\mathbb R}}
\newcommand\var{\op{var}}
\newcommand\cov{\op{cov}}
\newcommand\ybar{{\bar y}}
\newcommand\sigmahat{{\hat\sigma}}
$

<!-- #region toc=true -->
<h1>目次<span class="tocSkip"></span></h1>
<div class="toc"><ul class="toc-item"><li><span><a href="#正規分布モデルの共役事前分布とその応用" data-toc-modified-id="正規分布モデルの共役事前分布とその応用-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>正規分布モデルの共役事前分布とその応用</a></span><ul class="toc-item"><li><span><a href="#逆ガンマ正規分布" data-toc-modified-id="逆ガンマ正規分布-1.1"><span class="toc-item-num">1.1&nbsp;&nbsp;</span>逆ガンマ正規分布</a></span></li><li><span><a href="#共役事前分布のBayes更新" data-toc-modified-id="共役事前分布のBayes更新-1.2"><span class="toc-item-num">1.2&nbsp;&nbsp;</span>共役事前分布のBayes更新</a></span></li><li><span><a href="#μの周辺事前・事後分布および事前・事後予測分布" data-toc-modified-id="μの周辺事前・事後分布および事前・事後予測分布-1.3"><span class="toc-item-num">1.3&nbsp;&nbsp;</span>μの周辺事前・事後分布および事前・事後予測分布</a></span></li><li><span><a href="#Jeffreys事前分布の場合" data-toc-modified-id="Jeffreys事前分布の場合-1.4"><span class="toc-item-num">1.4&nbsp;&nbsp;</span>Jeffreys事前分布の場合</a></span></li><li><span><a href="#平均と対数分散について一様な事前分布の場合" data-toc-modified-id="平均と対数分散について一様な事前分布の場合-1.5"><span class="toc-item-num">1.5&nbsp;&nbsp;</span>平均と対数分散について一様な事前分布の場合</a></span></li><li><span><a href="#通常の信頼区間と予測区間との比較" data-toc-modified-id="通常の信頼区間と予測区間との比較-1.6"><span class="toc-item-num">1.6&nbsp;&nbsp;</span>通常の信頼区間と予測区間との比較</a></span></li><li><span><a href="#データの数値から事前分布を決めた場合" data-toc-modified-id="データの数値から事前分布を決めた場合-1.7"><span class="toc-item-num">1.7&nbsp;&nbsp;</span>データの数値から事前分布を決めた場合</a></span></li></ul></li></ul></div>
<!-- #endregion -->

```julia
using Distributions
using StatsPlots
default(fmt=:png, size=(400, 250),
    titlefontsize=10, tickfontsize=6, guidefontsize=9,
    plot_titlefontsize=10)
using SymPy
```

```julia
# Override the Base.show definition of SymPy.jl:
# https://github.com/JuliaPy/SymPy.jl/blob/29c5bfd1d10ac53014fa7fef468bc8deccadc2fc/src/types.jl#L87-L105

@eval SymPy function Base.show(io::IO, ::MIME"text/latex", x::SymbolicObject)
    print(io, as_markdown("\\displaystyle " *
            sympy.latex(x, mode="plain", fold_short_frac=false)))
end
@eval SymPy function Base.show(io::IO, ::MIME"text/latex", x::AbstractArray{Sym})
    function toeqnarray(x::Vector{Sym})
        a = join(["\\displaystyle " *
                sympy.latex(x[i]) for i in 1:length(x)], "\\\\")
        """\\left[ \\begin{array}{r}$a\\end{array} \\right]"""
    end
    function toeqnarray(x::AbstractArray{Sym,2})
        sz = size(x)
        a = join([join("\\displaystyle " .* map(sympy.latex, x[i,:]), "&")
                for i in 1:sz[1]], "\\\\")
        "\\left[ \\begin{array}{" * repeat("r",sz[2]) * "}" * a * "\\end{array}\\right]"
    end
    print(io, as_markdown(toeqnarray(x)))
end
```

## 正規分布モデルの共役事前分布とその応用


### 逆ガンマ正規分布

平均 $\mu\in\R$, 分散 $v=\sigma^2\in\R_{>0}$ の正規分布の確率密度函数を次のように表す:

$$
p_\op{Normal}(y|\mu, v) =
\frac{1}{\sqrt{2\pi v}}\exp\left(-\frac{1}{2v}(y-\mu)^2\right)
\quad (y\in \R).
$$

分散パラメータ $\sigma^2$ を $v$ に書き直している理由は, $\sigma^2$ を1つの変数として扱いたいからである.

パラメータ $\kappa, \theta > 0$ の逆ガンマ分布の確率密度函数を次のように書くことにする:

$$
p_\op{InverseGamma}(v|\kappa,\theta) =
\frac{\theta^\kappa}{\Gamma(\kappa)}
v^{-\kappa-1}\exp\left(-\frac{\theta}{v}\right)
\quad (v > 0).
$$

$v$ がこの逆ガンマ分布に従う確率変数だとすると, 

$$
\begin{aligned}
&
\frac{1}{v} \sim
\op{Gamma}\left(\kappa,\, \frac{1}{\theta}\right) =
\frac{1}{2\theta}\op{Gamma}\left(\frac{2\kappa}{2},\, 2\right) =
\frac{1}{2\theta}\op{Chisq}(2\kappa),
\\ &
E[v] = \frac{\theta}{\kappa - 1}, \quad
\var(v) = \frac{E[v]^2}{\kappa - 2}.
\end{aligned}
$$

$A$ と $B$ が $\mu, v$ に関する定数因子の違いを除いて等しいことを $A\propto B$ と書くことにする.

逆ガンマ正規分布の密度函数を次のように定義する:

$$
\begin{aligned}
p_\op{InverseGammaNormal}(\mu,v|\mu_*, v_*, \kappa, \theta) &=
p_\op{Normal}(\mu|\mu_*, v_* v) p_\op{InverseGamma}(v|\kappa, \theta)
\\ &\propto
v^{-(\kappa+1/2)-1}
\exp\left(-\frac{1}{v}\left(\theta + \frac{1}{2v_*}(\mu-\mu_*)^2\right)\right).
\end{aligned}
$$

この逆ガンマ正規分布の密度函数に従う確率変数を $\mu,v$ と書くと,

$$
E[v] = \frac{\theta}{\kappa-1}, \quad
\var(v) = \frac{E[v]^2}{\kappa-2}, \quad
\cov(\mu, v) = 0, \quad
E[\mu] = \mu_*, \quad
\var(\mu) = v_* E[v].
$$

この逆ガンマ正規分布が正規分布の共役事前分布になっていることを次の節で確認する. 


### 共役事前分布のBayes更新

データの数値 $y_1,\ldots,y_n$ が与えられたとき, 正規分布モデルの尤度函数は

$$
\prod_{i=1}^n p_\op{Normal}(y_i|\mu,v) \propto
v^{-n/2}\exp\left(-\frac{1}{2v}\sum_{i=1}^n(y_i-\mu)^2\right)
$$

の形になる. このとき,

$$
\ybar = \frac{1}{n}\sum_{i=1}^n y_i, \quad
\sigmahat^2 = \frac{1}{n}\sum_{i=1}^n (y_i - \ybar)^2.
$$

とおくと, 

$$
\sum_{i=1}^n(y_i-\mu)^2 = n(\mu - \ybar)^2 + n\sigmahat^2
$$

なので, 尤度を最大化する $\mu, v$ は $\mu=\ybar$, $v=\sigmahat^2$ になることがわかる.

さらに, 次が成立することもわかる:

$$
\begin{aligned}
&
\prod_{i=1}^n p_\op{Normal}(y_i|\mu,v)\times
p_\op{InverseGammaNormal}(\mu,v|\mu_*, v_*, \kappa, \theta)
\\ &\propto
v^{-n/2}\exp\left(-\frac{n}{2v}\left((\mu-\ybar)^2 + \sigmahat^2\right)\right)\times
v^{-(\kappa+1/2)-1}
\exp\left(-\frac{1}{v}\left(\theta + \frac{1}{2v_*}(\mu-\mu_*)^2\right)\right)
\\ &=
v^{-(\kappa+n/2+1/2)-1}
\exp\left(-\frac{1}{v}\left(
\theta + \frac{n}{2}\left(\sigmahat^2 + \frac{(\ybar - \mu_*)^2}{1+nv_*}\right) +
\frac{1+nv_*}{2v_*}\left(\mu - \frac{\mu_*+nv_*\ybar}{1+nv_*}\right)^2
\right)\right).
\end{aligned}
$$

ゆえに共役事前分布から得られる事後分布のパラメータは次のようになる:

$$
\begin{alignedat}{2}
&
\tilde\kappa = \kappa + \frac{n}{2} =
\frac{n}{2}\left(1 + \frac{2\kappa}{n}\right), 
\\ &
\tilde\theta =
\theta + \frac{n}{2}\left(\sigmahat^2 + \frac{(\ybar - \mu_*)^2}{1+nv_*}\right) =
\frac{n\sigmahat^2}{2}\left(1 + \frac{2\theta}{n\sigmahat^2} + \frac{(\ybar - \mu_*)^2}{(1+nv_*)\sigmahat^2}\right),
\\ &
\tilde\mu_* = \frac{\mu_*+nv_*\ybar}{1+nv_*} =
\ybar\frac{1+\mu_*/(nv_*\ybar)}{1+1/(nv_*)}, 
\\ &
\tilde v_* = \frac{v_*}{1+nv_*} =
\frac{1}{n}\frac{1}{1+1/(nv_*)}.
\end{alignedat}
$$

```julia
@vars n ȳ v̂ μ v μ0 v0 κ θ
```

```julia
negloglik = n/2*log(v) + n/(2v)*((μ - ȳ)^2 + v̂)
```

```julia
neglogpri = (κ + 1//2 + 1)*log(v) + 1/v*(θ + 1/(2v0)*(μ-μ0)^2)
```

```julia
neglogpost = (κ + n/2 + 1//2 + 1)*log(v) +  1/v*(
    θ + n/2*(v̂ + 1/(1+n*v0)*(ȳ - μ0)^2) +
    (1 + n*v0)/(2v0)*(μ - (μ0 + n*v0*ȳ)/(1 + n*v0))^2)
```

```julia
simplify(negloglik + neglogpri - neglogpost)
```

### μの周辺事前・事後分布および事前・事後予測分布

確率密度函数

$$
p(\mu|\mu_*,v_*,\kappa,\theta) =
\int_{\R_{>0}} p_\op{InverseGammaNormal}(\mu,v|\mu_*,v_*,\kappa,\theta) \,dv
$$

で定義される$\mu$ の周辺事前分布は次になる:

$$
\mu \sim
\mu_* + \sqrt{\frac{\theta}{\kappa}v_*}\;\op{TDist}(2\kappa).
$$

確率密度函数

$$
p_*(y_\op{new}|\mu_*,v_*,\kappa,\theta) =
\iint_{\R\times\R_{>0}}
p_\op{Normal}(y_\op{new}|\mu,v)
p_\op{InverseGammaNormal}(\mu,v|\mu_*,v_*,\kappa,\theta)
\,d\mu\,dv
$$

で定義される $y_\op{new}$ の事前予測分布は次になる:

$$
y_\op{new} \sim
\mu_* + \sqrt{\frac{\theta}{\kappa}(1+v_*)}\;\op{TDist}(2\kappa).
$$

パラメータをBayes更新後のパラメータ

$$
\begin{alignedat}{2}
&
\tilde\kappa = \kappa + \frac{n}{2} =
\frac{n}{2}\left(1 + \frac{2\kappa}{n}\right), 
\\ &
\tilde\theta =
\theta + \frac{n}{2}\left(\sigmahat^2 + \frac{(\ybar - \mu_*)^2}{1+nv_*}\right) =
\frac{n\sigmahat^2}{2}\left(1 + \frac{2\theta}{n\sigmahat^2} + \frac{(\ybar - \mu_*)^2}{(1+nv_*)\sigmahat^2}\right),
\\ &
\tilde\mu_* = \frac{\mu_*+nv_*\ybar}{1+nv_*} =
\ybar\frac{1+\mu_*/(nv_*\ybar)}{1+1/(nv_*)}, 
\\ &
\tilde v_* = \frac{v_*}{1+nv_*} =
\frac{1}{n}\frac{1}{1+1/(nv_*)}.
\end{alignedat}
$$

に置き換えればこれは $\mu$ の周辺事後分布および事後予測分布になる.

その事後分布を使った区間推定の幅は

* $n$ が大きいほど狭くなる.
* $\kappa$ が大きいほど狭くなる.
* $\theta$ が大きいほど広くなる.
* $|\ybar - \mu_*|/\sigmahat$ が大きいほど広くなる.
* $|\ybar - \mu_*|/\sigmahat$ が大きくても, $v_*$ がさらに大きければ狭くなる.

<!-- #region -->
### Jeffreys事前分布の場合


パラメータ空間が $\{(\mu, v)=(\mu, \sigma^2)\in\R\times\R_{>0}\}$ の $2$ 次元の正規分布モデルのJeffreys事前分布 $p_\op{Jeffreys}(\mu,v)$ は

$$
p_\op{Jeffreys}(\mu,v) \propto v^{-3/2}
$$

になることが知られている. ただし, 右辺の $(\mu,v)\in\R\times\R_{>0}$ に関する積分は $\infty$ になるので, この場合のJeffreys事前分布はimproperである.

逆ガンマ正規分布の密度函数

$$
p_\op{InverseGammaNormal}(\mu,v|\mu_*, v_*, \kappa, \theta) \propto
v^{-(\kappa+1/2)-1}
\exp\left(-\frac{1}{v}\left(\theta + \frac{1}{2v_*}(\mu-\mu_*)^2\right)\right).
$$

と比較すると, Jeffreys事前分布に対応する共役事前分布のパラメータ値は形式的に次になることがわかる:

$$
\kappa \to 0, \quad
\theta \to 0, \quad
v_* \to \infty.
$$

そのとき, Bayes更新後のパラメータの公式は次のようにシンプルになる:

$$
\tilde\kappa = \frac{n}{2}, \quad
\tilde\theta = \frac{n\sigmahat^2}{2}, \quad
\tilde\mu_* = \ybar, \quad
\tilde v_* = \frac{1}{n}.
$$

さらに, 前節の公式から, $n\to\infty$ のとき, 一般のパラメータ値に関するBayes更新の結果は, $n\to\infty$ のとき漸近的にこのJeffreys事前分布の場合に一致する.

さらに, Jeffreys事前分布の場合には

$$
\frac{\tilde\theta}{\tilde\kappa} = \sigmahat^2, \quad
\tilde v_* = \frac{1}{n}, \quad
2\tilde\kappa = n.
$$

ゆえに, $\mu$ に関する周辺事後分布は

$$
\mu \sim
\ybar + \frac{\sigmahat}{\sqrt{n}}\;\op{TDist}(n)
$$

になり, 事後予測分布は次になる:

$$
y_\op{new} \sim
\ybar + \sigmahat\sqrt{1+\frac{1}{n}}\;\op{TDist}(n).
$$
<!-- #endregion -->

### 平均と対数分散について一様な事前分布の場合

平均 $\mu$ と分数の対数 $\log v = \log\sigma^2$ に関する一様な事前分布は

$$
p_\op{flat}(\mu, v) \propto v^{-1}
$$

になる. ただし, 右辺の $(\mu,v)\in\R\times\R_{>0}$ に関する積分は $\infty$ になるので, この事前分布はimproperである.

逆ガンマ正規分布の密度函数

$$
p_\op{InverseGammaNormal}(\mu,v|\mu_*, v_*, \kappa, \theta) \propto
v^{-(\kappa+1/2)-1}
\exp\left(-\frac{1}{v}\left(\theta + \frac{1}{2v_*}(\mu-\mu_*)^2\right)\right).
$$

と比較すると, 平均と対数分散について一様な事前分布に対応する共役事前分布のパラメータ値は形式的に次になることがわかる:

$$
\kappa \to -\frac{1}{2}, \quad
\theta \to 0, \quad
v_* \to \infty.
$$

このとき, Bayes更新後のパラメータの公式は次のようになる:

$$
\tilde\kappa = \frac{n-1}{2}, \quad
\tilde\theta = \frac{n\sigmahat^2}{2}, \quad
\tilde\mu_* = \ybar, \quad
\tilde v_* = \frac{1}{n}.
$$

この場合には

$$
\frac{\tilde\theta}{\tilde\kappa} = \frac{n\sigmahat^2}{n-1} = s^2, \quad
\tilde v_* = \frac{1}{n}, \quad
2\tilde\kappa = n-1.
$$

ここで, $s^2$ はデータの数値 $y_1,\ldots,y_n$ の不偏分散

$$
s^2 = \frac{1}{n-1}\sum_{i=1}^n(y_i - \ybar)^2 =
\frac{n\sigmahat^2}{n-1} > \sigmahat^2
$$

であり, $s$ はその平方根である.

ゆえに, $\mu$ に関する周辺事後分布は

$$
\mu \sim
\ybar + \frac{s}{\sqrt{n}}\;\op{TDist}(n-1)
$$

になり, $y_\op{new}$ に関する事後予測分布は次になる:

$$
y_\op{new} \sim
\ybar + s\sqrt{1+\frac{1}{n}}\;\op{TDist}(n-1).
$$

したがって, 前節の結果と比較すると, Jeffreys事前分布の事後分布と予測分布による区間推定よりもこの場合の区間推定は少し広くなる.


### 通常の信頼区間と予測区間との比較

通常の $t$ 分布を使う平均の信頼区間と次の値の予測区間の構成では以下を使う:

$$
\frac{\ybar - \mu}{s\big/\!\sqrt{n}} \sim
\op{TDist}(n-1), \quad
\frac{y_\op{new} - \ybar}{s\sqrt{1+1/n}} \sim
\op{TDist}(n-1).
$$

ここで, $s^2$ はデータの数値の不偏分散であり, $s$ はその平方根である.

したがって, 前節の結果と比較すると, 通常の信頼区間と予測区間は, 平均と対数分散に関する一様事前分布に関する事後分布と予測分布を用いた区間推定に一致する.


### データの数値から事前分布を決めた場合

$a,b>0$ であると仮定する.

データの数値から共役事前分布のパラメータを次の条件によって決めたと仮定する:

$$
E[\mu] = \mu_* = \ybar, \quad
E[v] = \frac{\theta}{\kappa-1} = \sigmahat^2, \quad
\var(\mu) = v_* E[v] = a\sigmahat^2, \quad
\var(v) = \frac{E[v]^2}{\kappa-2} = b\sigmahat^4.
$$

これは次と同値である:

$$
\mu_* = \ybar, \quad
v_* = a, \quad
\kappa = 2 + \frac{1}{b}, \quad
\theta = \sigmahat^2\left(1 + \frac{1}{b}\right).
$$

これのBayes更新の結果は以下のようになる:

$$
\begin{alignedat}{2}
&
\tilde\kappa = 2 + \frac{1}{b} + \frac{n}{2} =
\frac{n}{2}\left(1 + \frac{2(2+1/b)}{n}\right)
& &
\to 2 + \frac{n}{2},
\\ &
\tilde\theta =
\sigmahat^2\left(1 + \frac{1}{b} + \frac{n}{2}\right) + \frac{n}{2}\frac{(\ybar - \ybar)^2}{1+na} =
\frac{n\sigmahat^2}{2}\left(1 + \frac{2(1+1/b))}{n}\right)
& &
\to \sigmahat^2\left(1 + \frac{n}{2}\right),
\\ &
\tilde\mu_* = \frac{\ybar+nv_*\ybar}{1+nv_*} =
\ybar
& &
\to \ybar,
\\ &
\tilde v_* = \frac{a}{1+na} =
\frac{1}{n}\frac{1}{1+1/(na)}
& &
\to \frac{1}{n}.
\end{alignedat}
$$

以上における $\to$ は $a\to\infty$, $b\to\infty$ での極限を意味する.

以上の構成のポイントは, $\mu_* = \ybar$ となっているおかげで, $\tilde\mu_*$ も $\tilde\mu_* = \ybar$ となってバイアスが消え, さらに, $\tilde\theta$ の中の $\ds\frac{n}{2}\frac{(\ybar - \mu_*)^2}{1+na}$ の項が消えて, 区間推定の幅が無用に広くならずに済むことである.

ただし, この場合には 

$$
\frac{\tilde\theta}{\tilde\kappa} =
\sigmahat^2\frac{1 + 2(1+1/b)/n}{1 + 2(2+1/b)/n} < \sigmahat^2, \quad
v_* = \frac{1}{n}\frac{1}{1+1/(na)} < \frac{1}{n}
$$

なので, 区間推定の幅はJeffreys事前分布の場合よりも少し狭くなる.

しかし, $n$ が大きければそれらの違いは小さくなる.

```julia

```
