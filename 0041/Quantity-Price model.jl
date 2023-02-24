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
#     display_name: Julia 1.9.0-beta4
#     language: julia
#     name: julia-1.9
# ---

# %% [markdown]
# $
# \newcommand\op{\operatorname}
# \newcommand\var{\op{var}}
# \newcommand\cov{\op{cov}}
# \newcommand\Normal{\op{Normal}}
# \newcommand\MvNormal{\op{MvNormal}}
# $
#
# https://escholarship.org/content/qt0j50v7hx/qt0j50v7hx_noSplash_f8d003538b36370a6a9c11bfe39ae39a.pdf で取り上げられているモデルは
#
# $$
# \begin{aligned}
# &
# Q = b_1 P + d_1 I + U_1,
# \\ &
# P = b_2 Q + d_2 W + U_2,
# \\ &
# \text{$\{I,W\}$ and $\{U_1,U_2\}$ are independent},
# \\ &
# \begin{bmatrix}
# E[U_1] \\
# E[U_2] \\
# \end{bmatrix} = \mu,
# \quad
# \begin{bmatrix}
# \var(U_1) & \cov(U_1,U_2) \\
# \cov(U_2,U_1) & \var(U_2) \\
# \end{bmatrix} = \Sigma.
# \end{aligned}
# $$
#
# この条件のうち前者の2つは次と同値である:
#
# $$
# Q = \frac{d_1 I + U_1 + b_1(d_2 W + U_2)}{1 - b_1 b_2}, \quad
# P = \frac{d_2 W + U_2 + b_2(d_1 I + U_1)}{1 - b_1 b_2}.
# $$
#
# このように生成したデータの数値からパラメータ値 $b_1,b_2,d_1,d_2,\mu,\Sigma$ を推定する方法を考える.
#
# $$
# \begin{bmatrix}
# A & C \\
# B & D \\
# \end{bmatrix}
# =
# T_{11}^{-1} T_{12},
# \quad
# T_{11}
# =
# \begin{bmatrix}
# \var(I)   & \cov(I,W) \\
# \cov(W,I) & \var(W) \\
# \end{bmatrix},
# \quad
# T_{12}
# =
# \begin{bmatrix}
# \cov(I,Q) & \cov(I,P) \\
# \cov(W,Q) & \cov(W,P) \\
# \end{bmatrix},
# $$
#
# とおくと,
#
# $$
# A = \frac{d_1}{1 - b_1 b_2}, \quad
# B = \frac{b_1 d_2}{1 - b_1 b_2}, \quad
# C = \frac{b_2 d_1}{1 - b_1 b_2}, \quad
# D = \frac{d_2}{1 - b_1 b_2}.
# $$
#
# すなわち,
#
# $$
# b_1 = B/D, \quad
# b_2 = C/A, \quad
# d_1 = (1 - b_1 b_2)A, \quad
# d_2 = (1 - b_1 b_2)D.
# $$
#
# そして,
#
# $$
# U_1 = Q - (b_1 P + d_1 I), \quad
# U_2 = P - (b_2 Q + d_2 W).
# $$
#
# もしくは,
#
# $$
# V_1 = Q - (AI + BW), \quad
# V_2 = P - (CI + DW)
# $$
#
# とおくと,
#
# $$
# U_1 = V_1 - b_1 V_2, \quad
# U_2 = V_2 - b_2 V_1, \quad
# \mu =
# \begin{bmatrix}
# E[U_1] \\
# E[U_2] \\
# \end{bmatrix},
# \quad
# \Sigma =
# \begin{bmatrix}
# \var(U_1) & \cov(U_1,U_2) \\
# \cov(U_2,U_1) & \var(U_2) \\
# \end{bmatrix}.
# $$
#
# データの数値に対して以上の手続きを適用すれば, パラメータ値 $b_1,b_2,d_1,d_2,\mu,\Sigma$ を推定できる.

# %%
using Distributions
using Random
#using StatsPlots
#default(fmt=:png, titlefontsize=8, tickfontsize=6, guidefontsize=8, legendfontsize=8)

# %%
function rand_qpiw(b1, b2, d1, d2, s1, s2)
    u1, u2, i, w = randn(), randn(), randn(), randn()
    q0 = d1*i + s1*u1
    p0 = d2*w + s2*u2
    q = (q0 + b1*p0)/(1 - b1*b2)
    p = (p0 + b2*q0)/(1 - b1*b2)
    q, p, i, w
end

function rand_qpiw!(b1, b2, d1, d2, s1, s2, qpiw::Matrix)
    Threads.@threads for i in axes(qpiw, 2)
        qpiw[:, i] .= rand_qpiw(b1, b2, d1, d2, s1, s2)
    end
    qpiw
end

function rand_qpiw!(b1, b2, d1, d2, s1, s2, n::Integer)
    qpiw = Matrix{Float64}(undef, 4, n)
    rand_qpiw!(b1, b2, d1, d2, s1, s2, qpiw)
end

# %%
b1, b2, d1, d2, s1, s2 = -1, 1, 0.4, 0.6, √0.5, √0.8
n = 10^6
qpiw = rand_qpiw!(b1, b2, d1, d2, s1, s2, n)
q, p, i, w = eachrow(qpiw);

# %%
X2 = [i w]
y2 = [q p]

β̂2 = X2 \ y2

# %%
[
       d1 b2*d1
    b1*d2    d2
]/(1 - b1*b2)

# %%
ww, xx, yy, zz = β̂2
bb1 = xx/zz
bb2 = yy/ww
dd1 = ww * (1 - bb1*bb2)
dd2 = zz * (1 - bb1*bb2)
[bb1, bb2, dd1, dd2]

# %%
u1 = q - (bb1*p + dd1*i)
u2 = p - (bb1*q + dd2*w)
fit(MvNormal, [u1 u2]')

# %%
mvn = fit(MvNormal, [i w q p]')
@show μ = mean(mvn)
Σ = cov(mvn)

# %%
Σ11 = Σ[1:2, 1:2]
Σ12 = Σ[1:2, 3:4]
Σ21 = Σ[3:4, 1:2]
Σ22 = Σ[3:4, 3:4]
β̂ = Σ11 \ Σ12

# %%
X = [i w]
q̂, p̂ = eachcol(X*β̂)
v1, v2 = q - q̂, p - p̂
u1, u2 = v1 - bb1*v2, v2 - bb2*v1
mvn = fit(MvNormal, [u1 u2]')
@show mean(mvn)
cov(mvn)

# %%
ww, xx, yy, zz = β̂
bb1 = xx/zz
bb2 = yy/ww
dd1 = ww * (1 - bb1*bb2)
dd2 = zz * (1 - bb1*bb2)
[bb1, bb2, dd1, dd2]

# %%
u1 = q - (bb1*p + dd1*i)
u2 = p - (bb2*q + dd2*w)
mvn = fit(MvNormal, [u1'; u2'])
@show mean(mvn)
cov(mvn)

# %%
