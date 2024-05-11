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
#     display_name: Julia 1.10.3
#     language: julia
#     name: julia-1.10
# ---

# %% [markdown]
# # 負の二項分布の連続時間極限としてガンマ分布が得られること
#
# * 黒木玄
# * 2024-05-11
# $
# \newcommand\op{\operatorname}
# $

# %% [markdown] toc=true
# <h1>目次<span class="tocSkip"></span></h1>
# <div class="toc"><ul class="toc-item"><li><span><a href="#解説" data-toc-modified-id="解説-1"><span class="toc-item-num">1&nbsp;&nbsp;</span>解説</a></span><ul class="toc-item"><li><span><a href="#負の二項分布の定義" data-toc-modified-id="負の二項分布の定義-1.1"><span class="toc-item-num">1.1&nbsp;&nbsp;</span>負の二項分布の定義</a></span></li><li><span><a href="#ガンマ分布の定義" data-toc-modified-id="ガンマ分布の定義-1.2"><span class="toc-item-num">1.2&nbsp;&nbsp;</span>ガンマ分布の定義</a></span></li><li><span><a href="#負の二項分布の連続時間極限がガンマ分布になること" data-toc-modified-id="負の二項分布の連続時間極限がガンマ分布になること-1.3"><span class="toc-item-num">1.3&nbsp;&nbsp;</span>負の二項分布の連続時間極限がガンマ分布になること</a></span></li><li><span><a href="#グラフを描くことによる確認方法" data-toc-modified-id="グラフを描くことによる確認方法-1.4"><span class="toc-item-num">1.4&nbsp;&nbsp;</span>グラフを描くことによる確認方法</a></span></li></ul></li><li><span><a href="#グラフによる確認" data-toc-modified-id="グラフによる確認-2"><span class="toc-item-num">2&nbsp;&nbsp;</span>グラフによる確認</a></span></li></ul></div>

# %% [markdown]
# ## 解説
#
# ### 負の二項分布の定義
#
# $r > 0$ と $0<p\le 1$ について, 負の二項分布 $\op{NegativeBinomial}(r, p)$ が次の確率質量函数によって定義される:
#
# $$
# p(m|r, p) = \binom{m + r - 1}{m} p^r (1 - p)^m \quad (m=0,1,2,\ldots).
# $$
#
# $r$ が正の整数のとき, $p(m|r, p)$ は当たりが出る確率が $p$ のルーレットを回し続けるとき, ちょうど $r$ 回の当たりが出るまでに出たはずれの回数が $m$ になる確率に等しい.
#
# $M$ が負の二項分布 $\op{NegativeBinomial}(r, p)$ に従う確率変数のとき, $M$ の期待値と分散は
#
# $$
# E[M] = \frac{r(1-p)}{p}, \quad
# \op{var}(M) = \frac{r(1-p)}{p^2}.
# $$
#
# ゆえに $L > 0$ について,
#
# $$
# E[M/L] = \frac{r(1-p)}{Lp}, \quad
# \op{var}(M/L) = \frac{r(1-p)}{L^2 p^2}.
# $$

# %% [markdown]
# ### ガンマ分布の定義
#
# $\alpha > 0$ と $\theta > 0$ について, ガンマ分布 $\op{Gamma}(\alpha, \theta)$ が次の確率密度函数によって定義される:
#
# $$
# p(t|\alpha, \theta) = \frac{1}{\theta^\alpha\Gamma(\alpha)} t^{\alpha-1} e^{-t/\theta} \quad (t > 0).
# $$
#
# $T$ がガンマ分布 $\op{Gamma}(\alpha, \theta)$ に従う確率変数のとき, $T$ の期待値と分散は
#
# $$
# E[T] = \alpha\theta, \quad
# \op{var}(T) = \alpha\theta^2.
# $$

# %% [markdown]
# ### 負の二項分布の連続時間極限がガンマ分布になること
#
# 以上を比較すると, $M\sim \op{NegativeBinomial}(r, p)$ に関する $M/L$ と $T \sim \op{Gamma}(\alpha, \theta)$ の期待値と分散は次の条件が成立するときに等しくなる:
#
# $$
# \alpha = r(1 - p), \quad \theta = \frac{1}{Lp}.
# $$
#
# これは次の条件と同値である:
#
# $$
# r = \frac{\alpha}{1 - 1/(L\theta)}, \quad p = \frac{1}{L\theta}.
# $$
#
# __定理:__ この条件の下で $L$ を十分に大きくすると, 負の二項分布の $1/L$ 倍 $\op{NegativeBinomial}(r, p)/L$ はガンマ分布 $\op{Gamma}(\alpha, \theta)$ で近似される.
#
# __証明:__ $r$ を上のようにおくと, $L\to\infty$ のとき $r\to \alpha$ となるので, $r=\alpha$ とおいて確認すれば十分である. さらに $p$ を上のようにおき, $t=m/L$ とおくと, 
#
# $$
# \binom{m+r-1}{m} =
# \frac{\Gamma(m+r)}{m!\Gamma(r)} = 
# \frac{\Gamma(m+r)}{m\Gamma(m)\Gamma(r)} = 
# \frac{1}{m B(r, m)}
# $$
#
# と
#
# $$
# \lim_{\beta\to\infty} \beta^{\alpha}B(\alpha, \beta) = \Gamma(\alpha)
# $$
#
# より, $L\to\infty$ のとき,
#
# $$
# \begin{aligned}
# p(m|r, p) \,dm &=
# \binom{Lt + \alpha - 1}{Lt} \frac{1}{L^\alpha \theta^\alpha} \left(1 - \frac{1}{L\theta} \right)^{Lt}\,L\,dt
# \\ &=
# \frac{1}{Lt\,B(\alpha, Lt)} \frac{1}{L^\alpha \theta^\alpha} \left(1 - \frac{1}{L\theta} \right)^{Lt}\,L\,dt
# \\ &=
# \frac{(Lt)^{\alpha-1}}{(Lt)^\alpha B(\alpha, Lt)} \frac{1}{L^\alpha \theta^\alpha} \left(1 - \frac{1}{L\theta} \right)^{Lt}\,L\,dt
# \\ &=
# \frac{1}{\theta^\alpha\,(Lt)^\alpha B(\alpha, Lt)} t^{\alpha-1} \left(1 - \frac{1}{L\theta} \right)^{Lt}\,dt
# \\ &\to
# \frac{1}{\theta^\alpha\Gamma(\alpha)}t^{\alpha-1}e^{-t/\theta}\,dt.
# \end{aligned}
# $$
#
# __証明終__
#
# __注意:__ $dm = L\,dt$ を使うイーカゲンな議論は以下のように正当化される. $f(t)$ を連続有界函数としたとき, $L\to\infty$ とすると, 
#
# $$
# \begin{aligned}
# \sum_{m=0}^\infty f\left(\frac{m}{L}\right)
# p(m|r, p)
# &=
# \sum_{m=0}^\infty f\left(\frac{m}{L}\right)
# \binom{m + \alpha - 1}{m} \frac{1}{L^\alpha \theta^\alpha} \left(1 - \frac{1}{L\theta} \right)^m
# \\ &=
# \sum_{m=0}^\infty f\left(\frac{m}{L}\right)
# \frac{1}{m\,B(\alpha, m)} \frac{1}{L^\alpha \theta^\alpha} \left(1 - \frac{1}{L\theta} \right)^m
# \\ &=
# \sum_{m=0}^\infty f\left(\frac{m}{L}\right)
# \frac{m^{\alpha-1}}{m^\alpha B(\alpha, m)} \frac{1}{L^\alpha \theta^\alpha} \left(1 - \frac{1}{L\theta} \right)^m
# \\ &=
# \sum_{m=0}^\infty f\left(\frac{m}{L}\right)
# \frac{1}{\theta^\alpha\,m^\alpha B(\alpha, m)} \left(\frac{m}{L}\right)^{\alpha-1} \left(1 - \frac{1}{L\theta} \right)^m\frac{1}{L}
# \\ &=
# \sum_{m=0}^\infty f\left(\frac{m}{L}\right)
# \frac{1}{\theta^\alpha\,(L(m/L))^\alpha B(\alpha, L(m/L))} \left(\frac{m}{L}\right)^{\alpha-1} \left(1 - \frac{1}{L\theta} \right)^{L(m/L)}\frac{1}{L}
# \\ &\to
# \int_0^\infty f(t)
# \frac{1}{\theta^\alpha\Gamma(\alpha)}t^{\alpha-1}e^{-t/\theta}\,dt.
# \end{aligned}
# $$
#
# 最後の行では,  $m/L$ を $t$ で $1/L$ を $dt$ で置き換え, $L\to\infty$ のとき $(Lt)^\alpha B(\alpha, Lt)\to\Gamma(\alpha)$, $(1-1/(L\theta))^{Lt}\to e^{-t/\theta}$ となることを使った.

# %% [markdown]
# ### グラフを描くことによる確認方法
#
# 前節の結果は大雑把には次のように言い直される:
#
# $$
# \alpha = r(1 - p), \quad \theta = \frac{1}{p} \quad
# \left(\Longleftrightarrow r = \frac{\alpha}{1 - 1/\theta}, \quad p = \frac{1}{\theta}\right)
# $$
#
# でかつ $p$ が十分小さいとき, 負の二項分布の $\op{NegativeBinomial}(r, p)$ のグラフはガンマ分布 $\op{Gamma}(\alpha, \theta)$ のグラフで近似される.
#
# このことをコンピュータで確認すればよい.

# %% [markdown]
# ## グラフによる確認

# %%
using Distributions
using StatsPlots
default(fmt=:png)

# %%
?NegativeBinomial

# %% [markdown]
# * μ = mean(NegativeBinomial(r, p)) = (1-p)*r/p
# * σ² = var(NegativeBinomial(r, p)) = (1-p)*r/p^2
# * α = μ^2/σ² = (1-p)*r
# * θ = σ²/μ = 1/p

# %%
mypdf(disy, x) = pdf(dist, x)
mypdf(dist::DiscreteUnivariateDistribution, x) = pdf(dist, round(x))

function plot_with_gammadist(dist; 
        μ = mean(dist), σ² = var(dist),
        α = μ^2/σ², θ = σ²/μ,
        a = -1.0, b = μ + 6√σ²,
        kwargs...
    )
    @show dist
    @show μ
    @show σ²
    @show gamma = Gamma(α, θ)
    plot(x -> mypdf(dist, x), a, b; label="dist")
    plot!(x -> pdf(gamma, x), a, b; label="gamma", ls=:dash)
    plot!(; kwargs...)
end

# %%
plot_with_gammadist(NegativeBinomial(1, 0.8); ylim=(-0.1, 5))

# %%
plot_with_gammadist(NegativeBinomial(1, 0.3))

# %%
plot_with_gammadist(NegativeBinomial(1, 0.1))

# %%
plot_with_gammadist(NegativeBinomial(1, 0.03))

# %%
plot_with_gammadist(NegativeBinomial(1, 0.01))

# %%
plot_with_gammadist(NegativeBinomial(2, 0.8); ylim=(-0.1, 5))

# %%
plot_with_gammadist(NegativeBinomial(2, 0.3))

# %%
plot_with_gammadist(NegativeBinomial(2, 0.1))

# %%
plot_with_gammadist(NegativeBinomial(2, 0.03))

# %%
plot_with_gammadist(NegativeBinomial(5, 0.8))

# %%
plot_with_gammadist(NegativeBinomial(5, 0.3))

# %%
plot_with_gammadist(NegativeBinomial(5, 0.1))

# %%
