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
# # 潜在反応変数を含むダイアグラムの描き方について
#
# * 黒木玄
# * 2023-12-11
#
# __要約:__ \[黒木学 2017\] の図6-1, 図6-2, 図6-3, 図8-4, 図8-5では $Y_x$ のように表記された潜在反応変数から元の変数 $Y$ に向けて矢線が描かれている.  その部分は, シンプルな標準的ルールに従って, $Y_x$ から $Y$ に矢線を描かない別のダイアグラムに描き変えた方がよい.
#
# __参照文献__
#
# * \[黒木学 2017\] 黒木学, 構造的因果モデルの基礎, 共立出版, 2017年, iv + 309 pages.
# * \[Pearl 2009\] Judea Pearl, Causality --- Models, Resoning and Inference, Second Edition, Cambridge University Press, 2009, xix + 464 pages.
# * \[Shpister--Pearl 2007\] Ilya Shpitser and Judea Pearl, What counterfactuals can be tested, In: Proceedings of the 23rd Conference on Uncertainty in Artificial Intelligence, UAI 2007, 352--359. [https://arxiv.org/abs/1206.5294](https://arxiv.org/abs/1206.5294)

# %% [markdown] toc=true
# <h1>目次<span class="tocSkip"></span></h1>
# <div class="toc"><ul class="toc-item"><li><span><a href="#[黒木学-2017]-の図6-1,-図6-2,-図6-3,-図8-4,-図8-5およびそれらの訂正案" data-toc-modified-id="[黒木学-2017]-の図6-1,-図6-2,-図6-3,-図8-4,-図8-5およびそれらの訂正案-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>[黒木学 2017] の図6-1, 図6-2, 図6-3, 図8-4, 図8-5およびそれらの訂正案</a></span><ul class="toc-item"><li><span><a href="#図6-1とその訂正案" data-toc-modified-id="図6-1とその訂正案-1.1"><span class="toc-item-num">1.1&nbsp;&nbsp;</span>図6-1とその訂正案</a></span></li><li><span><a href="#図6-2とその訂正案" data-toc-modified-id="図6-2とその訂正案-1.2"><span class="toc-item-num">1.2&nbsp;&nbsp;</span>図6-2とその訂正案</a></span></li><li><span><a href="#図6-3とその訂正案" data-toc-modified-id="図6-3とその訂正案-1.3"><span class="toc-item-num">1.3&nbsp;&nbsp;</span>図6-3とその訂正案</a></span></li><li><span><a href="#図8-4とその訂正案" data-toc-modified-id="図8-4とその訂正案-1.4"><span class="toc-item-num">1.4&nbsp;&nbsp;</span>図8-4とその訂正案</a></span></li><li><span><a href="#図8-5とその訂正案" data-toc-modified-id="図8-5とその訂正案-1.5"><span class="toc-item-num">1.5&nbsp;&nbsp;</span>図8-5とその訂正案</a></span></li></ul></li><li><span><a href="#訂正した方が良い理由" data-toc-modified-id="訂正した方が良い理由-2"><span class="toc-item-num">2&nbsp;&nbsp;</span>訂正した方が良い理由</a></span></li><li><span><a href="#他の文献でのダイアグラムの描き方の例" data-toc-modified-id="他の文献でのダイアグラムの描き方の例-3"><span class="toc-item-num">3&nbsp;&nbsp;</span>他の文献でのダイアグラムの描き方の例</a></span></li></ul></div>

# %% [markdown]
# ## \[黒木学 2017\] の図6-1, 図6-2, 図6-3, 図8-4, 図8-5およびそれらの訂正案
#
# この節では \[黒木学 2017\] の図6-1, 図6-2, 図6-3, 図8-4, 図8-5およびその訂正案を示す.

# %% [markdown]
# 訂正の仕方は非常にシンプルである. 例えば, 構造的因果モデル内で
#
# $$
# Z = g_z(A, B, C, D, \epsilon_z), \quad
# Z_{a,b} = g_z(a, b, C_{a,b}, D_{a,b}, \epsilon_z)
# $$
#
# と表示されている変数 $Z$, $Z_{a,b}$ について, ダイアグラムを次のように描くだけである:
#
# * $A,B,C,D,\epsilon_z$ から $Z$ への矢線を描き, $Z$ に向けての矢線はその5本だけにする.
# * 特に $Z_{a,b}$ から $Z$ への矢線は描かない.
# * 介入で固定した値 $a,b$ をダイアグラムに描かない場合には, $C_{a,b}, D_{a,b}, \epsilon_z$ から $Z_{a,b}$ への矢線を描き, $Z_{a,b}$ に向けての矢線はその3本だけにする.
# * 例えば $C = C_{a,b}$ となっている場合には変数 $C_{a,b}$ を新たにダイアグラムに追加せずに $C_{a,b}$ の代わりに $C$ を用いる.
#
# このルールは標準的である.

# %% [markdown]
# ### 図6-1とその訂正案

# %% [markdown]
# <img src="https://raw.githubusercontent.com/genkuroki/public/main/0046/images/fig6-1_corrected.jpg" width=70%>

# %% [markdown]
# ### 図6-2とその訂正案

# %% [markdown]
# <img src="https://raw.githubusercontent.com/genkuroki/public/main/0046/images/fig6-2_corrected.jpg" width=70%>

# %% [markdown]
# ### 図6-3とその訂正案

# %% [markdown]
# <img src="https://raw.githubusercontent.com/genkuroki/public/main/0046/images/fig6-3_corrected.jpg" width=70%>

# %% [markdown]
# ### 図8-4とその訂正案

# %% [markdown]
# <img src="https://raw.githubusercontent.com/genkuroki/public/main/0046/images/fig8-4_corrected.jpg" width=70%>

# %% [markdown]
# ### 図8-5とその訂正案

# %% [markdown]
# <img src="https://raw.githubusercontent.com/genkuroki/public/main/0046/images/fig8-5_corrected.jpg">

# %% [markdown]
# ## 訂正した方が良い理由

# %% [markdown]
# 構造的因果モデルを表すための等式 $Z = g_z(A, B, C, D, \epsilon_z)$ では, 右辺の表示自体が重要であり, 等式を保ったまま右辺を勝手に変形してはいけない.  例えば, $B = g_b(A, \epsilon_b)$ となっているとき, 右辺にそれを代入して $Z = g_z(A, g_b(A, \epsilon_b), C, D, \epsilon_z)$ と書き直して, $Z$ は $B$ に因果的に依存しないとみなして, 因果ダイアグラムを描くときに $B$ から $Z$ への矢線を描かないようにしてはいけない.

# %% [markdown]
# $Y = g_y(X, W, \epsilon_y)$ に対応する $X=x$ と設定した場合の潜在結果変数は $Y_x=g_y(x, W_x, \epsilon_y)$ のように定義される. $W_x$ は $Y$ から先祖をさかのぼる帰納的によって同様のルールで定義される. このとき,
#
# $$
# g_y(X, W, \epsilon_y) = g_y(x, W_x, \epsilon_y) \quad\text{if $X=x$}
# $$
#
# が自明に成立する.  すなわち
#
# $$
# Y = Y_x \quad \text{if $X=x$}
# $$
#
# が自明に成立する.  このことを理由に因果ダイアグラムで $Y_x$ から $Y$ に矢線を描くことは複数通りの意味で不適切になる.
#
# __理由1:__ 因果ダイアグラムの描き方はシンプルな標準的ルールとして決まっている. 標準的ルールを勝手に変更して読者を惑わせることはやらない方がよい.
#
# __理由2:__ 因果推論の文脈で $Y$ は観察研究の結果のモデル化だと解釈され, $Y_x$ は観察兼研究の結果と無関係に $X=x$ と介入した場合の結果のモデル化だとみなされる. $Y_x$ は $Y$ の原因の1つにはなっていないので, 因果推論の文脈では通常因果を意味する矢線を $Y_x$ から $Y$ に引かない方がよい.
#
# __理由3:__ この理由は非本質的だが, 読者の参考になるかもしれないので説明しておくことにする.
#
# 例として, 因果ダイアグラム $X\to Y$ で要約される構造的因果モデル
#
# $$
# X = g_x(\epsilon_x), \quad Y = g_y(X, \epsilon_y)
# $$
#
# を考える.  ここで $\epsilon_x$, $\epsilon_y$ は独立な確率変数である.  このとき, $X=x と介入した場合の結果のモデル化 $Y_x$ は次のように定義される:
#
# $$
# Y_x = g_y(x, \epsilon_y).
# $$
#
# $X, Y, Y_x$ で構成されるモデルに対応するダイアグラムは標準的ルールに従うと次のようになる:

# %% [markdown]
# <img src="https://raw.githubusercontent.com/genkuroki/public/main/0046/images/2a.jpg" width=35%>

# %% [markdown]
# 仮に $X$ の取り得る値は $1,0$ の2通りであるとする. そのとき
#
# $$
# Y = \begin{cases}
# Y_1 & \text{if X=1,} \\
# Y_0 & \text{if X=0} \\
# \end{cases}
# $$
#
# が自明に成立している. この右辺は $Y_1, Y_0$ の両方と $X$ で決まるので, この等式を
#
# $$
# Y = h(X, Y_1, Y_0)
# $$
#
# と要約しよう. これと
#
# $$
# X = g_x(\epsilon_x), \quad
# Y_1 = g_y(1, \epsilon_y), \quad Y_0 = g_y(0, \epsilon_y),
# $$
#
# と合わせて得られる構造的因果モデルを仮に考えたとする. そのような考え方はしない方がよいが仮にそう考えたとする. そのとき, そのモデルに対応するダイアグラムは次のようになる:

# %% [markdown]
# <img src="https://raw.githubusercontent.com/genkuroki/public/main/0046/images/2b.jpg" width=50%>

# %% [markdown]
# これを単独の $x=1,0$ に関する次のダイアグラムで要約することはできない:

# %% [markdown]
# <img src="https://raw.githubusercontent.com/genkuroki/public/main/0046/images/2c.jpg" width=25%>

# %% [markdown]
# $Y$ の挙動を確定させるためには $Y_1$, $Y_0$ の両方が必要である.
#
# 以上のように単独の $Y_x$ から $Y$ への矢線を描いたダイアグラムを描くことはモデルの要約の仕方として不完全な方法になる.

# %% [markdown]
# ## 他の文献でのダイアグラムの描き方の例

# %% [markdown]
# 次の図はこの話題について最も標準的な文献だとみなされる \[Pearl 2009\] Judea Pearl, Causality --- Models, Resoning and Inference, Second Edition, 2009 のp.394のFigure 11.18 である:

# %% [markdown]
# <img src="https://raw.githubusercontent.com/genkuroki/public/main/0046/images/3a.jpg" width=80%>

# %% [markdown]
# この図は読者からの質問の中で示されたダイアグラムであり, このダイアグラムによる解釈についてPearlは次のように肯定している:
#
# \[Pearl 2009], p.394より
# >Author’s Reply (with Ilya Shpitser):<br>
# >Your generalization of the twin network to more than two worlds is correct, and so is
# your conclusion; $Y_x$ is not independent of $X$ given $Y_z$, $Z_x$, $Y$. In fact, a recent paper
# (Shpitser and Pearl 2007) formulates this graphical representation of counterfactuals in
# multiple worlds under the rubric “counterfactual graphs.”

# %% [markdown]
# そこで引用されている\[Shpitser--Pearl 2007\]には次の図がある:

# %% [markdown]
# <img src="https://raw.githubusercontent.com/genkuroki/public/main/0046/images/3c.jpg" width=60%>

# %% [markdown]
# これらのダイアグラムでも $Y_x$ から $Y$ への矢線は描かれておらず, 私による訂正案と同じシンプルな標準的なルールに従ってダイアグラムが描かれている.
