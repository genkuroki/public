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
#     display_name: Julia 1.8.4
#     language: julia
#     name: julia-1.8
# ---

# %% [markdown]
# # 不完全Gauss積分と不完全ガンマ函数の積
#
# * 黒木玄
# * 2022-12-31

# %% [markdown] toc=true
# <h1>目次<span class="tocSkip"></span></h1>
# <div class="toc"><ul class="toc-item"><li><span><a href="#不完全Gauss積分の2乗" data-toc-modified-id="不完全Gauss積分の2乗-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>不完全Gauss積分の2乗</a></span></li><li><span><a href="#不完全ガンマ函数の積" data-toc-modified-id="不完全ガンマ函数の積-2"><span class="toc-item-num">2&nbsp;&nbsp;</span>不完全ガンマ函数の積</a></span></li><li><span><a href="#ガンマ分布とベータプライム分布との関係" data-toc-modified-id="ガンマ分布とベータプライム分布との関係-3"><span class="toc-item-num">3&nbsp;&nbsp;</span>ガンマ分布とベータプライム分布との関係</a></span></li></ul></div>

# %% [markdown]
# ## 不完全Gauss積分の2乗
#
# 次の積分を不完全Gauss積分と呼ぶことにする:
#
# $$
# F(t) = \int_0^t e^{-x^2}\,dx.
# $$
#
# 以下, $t>0$ であると仮定する.
#
# 次を示そう:
#
# $$
# F(t)^2
# = \int_0^1 \frac{1-e^{-(1+s^2)t^2}}{1+s^2}\,ds
# = \frac{\pi}{4} - \int_0^1 \frac{e^{-(1+s^2)t^2}}{1+s^2}\,ds.
# \tag{$*$}
# $$
#
# この公式の両辺の $t\to\infty$ での極限を見ればGauss積分の値を計算できる:
#
# $$
# \left(\int_0^\infty e^{-x^2}\,dx\right)^2 = \frac{\pi}{4}, \quad
# \text{i.e.}\; \int_0^\infty e^{-x^2}\,dx = \frac{\sqrt{\pi}}{2}.
# $$
#
# ($*$)を証明しよう.
#
# $$
# \begin{aligned}
# F(t)^2 &=
# \iint_{0<x,y<t} e^{-(x^2+y^2)}\,dx\,dy
# \end{aligned}
# $$
#
# の $x,y$ に関する対称性より,
#
# $$
# \begin{aligned}
# F(t)^2
# &= 2\iint_{0<x<y<t} e^{-(x^2+y^2)}\,dx\,dy
# = 2\int_0^t \left(\int_0^y e^{-(x^2+y^2)}\,dx\right)\,dy
# \\
# &= 2\int_0^t \left(\int_0^1 e^{-(1+s^2)y^2}y\,ds\right)\,dy
# = 2\int_0^1 \left(\int_0^t e^{-(1+s^2)y^2}y\,dy\right)\,ds
# \\
# &= \int_0^1 \left[\frac{-e^{-(1+s^2)y^2}}{1+s^2}\right]_{y=0}^{y=t}\,ds
# = \int_0^1 \frac{1-e^{-(1+s^2)t^2}}{1+s^2}\,ds
# \\
# &= \int_0^1 \frac{ds}{1+s^2} - \int_0^1 \frac{e^{-(1+s^2)t^2}}{1+s^2}\,ds
# = \frac{\pi}{4} - \int_0^1 \frac{e^{-(1+s^2)t^2}}{1+s^2}\,ds.
# \end{aligned}
# $$
#
# 3番目の等号で $x = sy$ とおき, 4番目の等号で正値の被積分函数に関する積分順序の交換を使った.  最後の等号では $\int_0^a ds/(1+s^2) = \arctan a$ を使った.  これで($*$)が証明された.
#
# 以上の計算は $yx$ 平面上の直角二等辺三角形型の領域 $0<x<y<t$ 上での積分の直線 $x=sy$ の傾きを意味する変数 $s$ を使った計算になっている.

# %% [markdown]
# ## 不完全ガンマ函数の積
#
# 不完全ガンマ函数の1つは
#
# $$
# \gamma(a, t) = \int_0^t e^{-x}x^{a-1}\,dx \quad(a > 0)
# $$
#
# と定義される. ガンマ函数はこれの $t\to\infty$ での極限に等しい:
#
# $$
# \Gamma(a) = \int_0^\infty e^{-x}x^{a-1}\,dx.
# $$
#
# ベータ函数 $B(a, b)$ は次の表示を持つのであった:
#
# $$
# B(a, b) = \int_0^1 t^{a-1}(1-t)^{b-1}\,dt
# = \int_0^\infty \frac{s^{a-1}\,ds}{(1+s)^{a+b}}
# = \int_0^\infty \frac{s^{b-1}\,ds}{(1+s)^{a+b}}.
# $$
#
# 2つめの等号は $t = s/(1+s)$ とおけば得られる.
#
# __注意1:__ 不完全ガンマ函数の定義で $t$, $x$ を $t^2$, $x^2$ で置き換えると,
#
# $$
# \gamma(a, t^2) = 2\int_0^t e^{-x^2} x^{2a-1}\,dx
# $$
#
# なので
#
# $$
# \gamma(1/2, t^2) = 2\int_0^t e^{-x^2}\,dx = 2F(t)
# $$
#
# と不完全Gauss積分が特殊な場合として得られる.  前節の結果より,
#
# $$
# \gamma(1/2, t^2)^2 = 4\int_0^1 \frac{1-e^{-(1+s^2)t^2}}{1+s^2}\,ds.
# $$
#
#
# __注意2:__ ガンマ函数の積については次の公式がよく知られている:
#
# $$
# \Gamma(a)\Gamma(b) = \Gamma(a+b)B(a, b).
# $$
#
# すなわち,
#
# $$
# \Gamma(a)\Gamma(b)
# = \int_0^\infty \Gamma(a+b)\frac{s^{a-1}}{(1+s)^{a+b}}\,ds.
# $$
#
# __問題:__ 以上の注意1,2の結果の両方を含む $\gamma(a,t)\gamma(b,u)$ に関する公式を得よ.
#
# 以下ではこの問題について考える.  $a,b,t,u>0$ と仮定する.

# %% [markdown]
# 長方形型の領域 $0<x<t$, $0<y<u$ を $(x,y)=(0,0)$ と $(x,y)=(t,u)$ を通る直線によって2つの直角三角形型の領域に分割できる. その直線の上側に $(x,y)$ が含まれることは $y/x > u/t$ すなわち $x < ty/u$ と同値であり, 下側に含まれることは $y/x < u/t$ すなわち $y < ux/t$ と同値である.  ゆえに
#
# $$
# \begin{aligned}
# \gamma(a,t)\gamma(b,u)
# &= \iint_{0<x<t,\; 0<y<u} e^{-(x+y)}x^{a-1}y^{b-1}\,dx\,dy
# = A + B.
# \end{aligned}
# $$
#
# ここで,
#
# $$
# A = \iint_{0<y<u,\; 0<x<ty/u} e^{-(x+y)}x^{a-1}y^{b-1}\,dx\,dy, \quad
# B = \iint_{0<x<t,\; 0<y<ux/t} e^{-(x+y)}x^{a-1}y^{b-1}\,dx\,dy.
# $$
#
# $B$ は $A$ における $a,t,x$ と $b,u,y$ を交換した結果に等しいので, $A$ に関する公式が得られれば $B$ に関する公式も得られる.
#
# $$
# \begin{aligned}
# A
# &= \int_0^u\left(\int_0^{ty/u} e^{-(x+y)}x^{a-1}y^{b-1} \,dx\right)\,dy
# = \int_0^u\left(\int_0^{t/u} e^{-(1+s)y}(sy)^{a-1}y^{b-1} y\,ds\right)\,dy
# \\
# &= \int_0^u\left(\int_0^{t/u} e^{-(1+s)y}y^{a+b-1}s^{a-1}\,ds\right)\,dy
# = \int_0^{t/u}\left(\int_0^u e^{-(1+s)y}y^{a+b-1}s^{a-1}\,dy\right)\,ds
# \\
# &= \int_0^{t/u}\left(\int_0^{(1+s)u} e^{-z}\frac{z^{a+b-1}}{(1+s)^{a+b-1}}s^{a-1}\frac{dz}{1+s}\right)\,ds
# \\
# &= \int_0^{t/u}\left(\int_0^{(1+s)u} e^{-z}z^{a+b-1}\frac{s^{a-1}}{(1+s)^{a+b}}\,dz\right)\,ds
# \\
# &= \int_0^{t/u} \gamma(a+b, (1+s)u) \frac{s^{a-1}}{(1+s)^{a+b}}\,ds
# \end{aligned}
# $$
#
# 2つめの等号で $x=sy$ とおいた. それは前節の不完全Gauss積分の2乗の計算でも使われたテクニックである. 4つめの等号で正値被積分函数に関する積分順序の交換を使った.  5つめの等号で $y = z/(1+s)$ とおき, 最後の等号で不完全ガンマ函数の定義を使った.
#
# $B$ は $A$ における $a,t,x$ と $b,u,y$ を交換した結果に等しいので, 次の公式が得られる:
#
# $$
# \gamma(a,t)\gamma(b,u)
# = \int_0^{t/u} \gamma(a+b, (1+s)u) \frac{s^{a-1}}{(1+s)^{a+b}}\,ds
# + \int_0^{u/t} \gamma(a+b, (1+s)t) \frac{s^{b-1}}{(1+s)^{a+b}}\,ds.
# \tag{$**$}
# $$
#
# これが問題の答えである.
#
# 公式($**$)で $a=b=1/2$ とおくと,
#
# $$
# \gamma(1, v) = \int_0^v e^{-x}\,dx = 1 - e^{-v}
# $$
#
# より,
#
# $$
# \begin{aligned}
# \gamma(1/2, t)\gamma(1/2, u)
# &= \int_0^{t/u} (1 - e^{-(1+s)u}) \frac{s^{-1/2}}{1+s}\,ds
# +  \int_0^{u/t} (1 - e^{-(1+s)t}) \frac{s^{-1/2}}{1+s}\,ds
# \\
# &= \int_0^{t/u} (1 - e^{-(1+s^2)u}) \frac{s^{-1}}{1+s^2}2s\,ds
# +  \int_0^{u/t} (1 - e^{-(1+s^2)t}) \frac{s^{-1}}{1+s^2}2s\,ds
# \\
# &= 2\int_0^{t/u} \frac{1 - e^{-(1+s^2)u}}{1+s^2}\,ds
# +  2\int_0^{u/t} \frac{1 - e^{-(1+s^2)t}}{1+s^2}\,ds.
# \end{aligned}
# $$
#
# 2つめの等号で $s$ を $s^2$ で置き換えた.
#
# さらに $t, u$ の両方に $t^2$ を代入すると, 注意1の結果が得られる:
#
# $$
# \gamma(1/2, t^2)^2
# = 4\int_0^1 \frac{1 - e^{-(1+s^2)t^2}}{1+s^2}\,ds.
# $$
#
# 公式($**$)で $t\to\infty$ の極限を取ると,
#
# $$
# \Gamma(a)\gamma(b,u)
# = \int_0^\infty \gamma(a+b, (1+s)u) \frac{s^{a-1}}{(1+s)^{a+b}}\,ds.
# $$
#
# さらに $u\to\infty$ の極限を取ると, 注意2の結果が得られる:
#
# $$
# \Gamma(a)\Gamma(b)
# = \int_0^\infty \Gamma(a+b) \frac{s^{a-1}}{(1+s)^{a+b}}\,ds.
# $$

# %% [markdown]
# ## ガンマ分布とベータプライム分布との関係
#
# $X, Y$ はそれぞれがガンマ分布 $\operatorname{Gamma}(a, 1)$, $\operatorname{Gamma}(b, 1)$ に従う独立な確率変数であると仮定する. 
#
# このとき, $Z = X+Y$, $P = X/(X+Y)$ とおくと, $Z, P$ は独立な確率変数になり, $Z$ はガンマ分布 $\operatorname{Gamma}(a+b, 1$ に従い, $P$ はベータ分布 $\operatorname{Beta}(a, b)$ に従うことを示せる.
#
# さらに, $S = P/(1-P) = X/Y$ とおくと, $S$ は[ベータプライム分布](https://en.wikipedia.org/wiki/Beta_prime_distribution) $\operatorname{BetaPrime}(a,b)$ に従う確率変数になることも示される. $S$ の逆数 $S^{-1}=Y/X$ は分布 $\operatorname{BetaPrime}(b,a)$ に従う. ($
#
# これらはガンマ分布とベータ分布及びベータプライム分布の基本的な関係であり, 統計学の分野ではよく知られている結果である.
#
# 以上の結果を使うと, $t,u>0$ のとき,
#
# $$
# \begin{aligned}
# &
# P(0<X<t)\,P(0<Y<u)
# \\
# &= P(0<X<t\;\&\;0<Y<u)
# \\
# &= P\left(0<Y<u \;\&\; 0<X<\frac{t}{u}Y\right) + P\left(0<X<t \;\&\; 0<Y<\frac{u}{t}X\right)
# \\
# &= P\left(0<Y<u \;\&\; 0<S<\frac{t}{u}\right) + P\left(0<X<t \;\&\; 0<S^{-1}<\frac{u}{t}\right)
# \\
# &= P\left(0<Z<(1+S)u \;\&\; 0<S<\frac{t}{u}\right) + P\left(0<Z<(1+S^{-1})t \;\&\; 0<S^{-1}<\frac{u}{t}\right).
# \end{aligned}
# $$
#
# そして,
#
# $$
# P(0<X<t) = \frac{\gamma(a,t)}{\Gamma(a)}, \quad
# P(0<Y<u) = \frac{\gamma(b,u)}{\Gamma(b)}
# $$
#
# でかつ
#
# $$
# \begin{aligned}
# P\left(0<Z<(1+S)u \;\&\; 0<S<\frac{t}{u}\right)
# &= \int_0^{t/u} P(0<Z<(1+s)u) \frac{1}{B(a,b)}\frac{s^{a-1}}{(1+s)^{a+b}}\,ds
# \\
# &= \int_0^{t/u} \frac{\gamma(a+b, (1+s)u)}{\Gamma(a+b)} \frac{1}{B(a,b)}\frac{s^{a-1}}{(1+s)^{a+b}}\,ds
# \end{aligned}
# $$
#
# であり, これにおける $t,a$ と $u,b$ を交換した結果は $P\left(0<Z<(1+S^{-1})t \;\&\; 0<S^{-1}<\frac{u}{t}\right)$ に等しい.  ゆえに,
#
# $$
# \frac{\gamma(a,t)}{\Gamma(a)}\frac{\gamma(b,u)}{\Gamma(b)}
# = \int_0^{t/u} \frac{\gamma(a+b, (1+s)u)}{\Gamma(a+b)} \frac{1}{B(a,b)}\frac{s^{a-1}}{(1+s)^{a+b}}\,ds
# + \int_0^{u/t} \frac{\gamma(a+b, (1+s)t)}{\Gamma(a+b)} \frac{1}{B(b,a)}\frac{s^{b-1}}{(1+s)^{a+b}}\,ds.
# $$
#
# $\Gamma(a)\Gamma(b)=\Gamma(a+b)B(a,b)$, $B(b,a)=B(a,b)$ より, これは前節で得た結果($**$)と同値である.
#
# 以上のように前節で得た結果は本質的にガンマ分布とベータプライム分布の基本的な関係を与えることと同じことになっている.

# %%
