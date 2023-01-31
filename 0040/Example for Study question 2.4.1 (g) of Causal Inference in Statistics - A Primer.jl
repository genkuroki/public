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
#     display_name: Julia 1.9.0-beta3
#     language: julia
#     name: julia-1.9
# ---

# %%
using Distributions
using StatsPlots
using Printf

# %% [markdown]
# Model:
#
# $$
# \begin{alignedat}{2}
# &
# Z_1 \sim \operatorname{Normal}(0, 1), & &
# \\ &
# Z_2 \sim \operatorname{Normal}(0, 1),
# \\ &
# Z_3 = Z_1 + Z_2 + \varepsilon_1, & & \quad \varepsilon_1 \sim \operatorname{Normal}(0, e),
# \\ &
# X = Z_1 + Z_3 + \varepsilon_2, & & \quad \varepsilon_2 \sim \operatorname{Normal}(0, e),
# \\ &
# W = X + \varepsilon_3, & & \quad \varepsilon_3 \sim \operatorname{Normal}(0, e),
# \\ &
# Y = W + Z_2 + Z_3 + \varepsilon_4, & & \quad \varepsilon_4 \sim \operatorname{Normal}(0, e) \qquad (e = 0.2)
# \end{alignedat}
# $$

# %%
function rand_fig2_9(; e=0.2)
    z1 = randn()
    z2 = randn()
    z3 = z1 + z2 + e*randn()
    x = z1 + z3 + e*randn()
    w = x + e*randn()
    y = w + z2 + z3 + e*randn()
    [w, x, y, z1, z2, z3]
end

function rand_fig2_9(n; e=0.2)
    [rand_fig2_9(; e) for _ in 1:n] |> stack
end

# %%
N = 10^4
data = rand_fig2_9(N)

# %%
scatter(data[6,:], data[5,:]; ms=1, msc=:auto, ma=0.5, label="(z3, z2)")

# %%
scatter(data[1,:], data[5,:]; ms=1, msc=:auto, ma=0.5, label="(w, z2)")

# %%
anim = @animate for z3 in -4:0.2:4
    ks = @.(abs(data[6,:] - z3) < 0.1)
    scatter(data[1,ks], data[5,ks]; ms=3, msc=:auto, ma=0.5, label="")
    plot!(xlim=extrema(data[1,:]), ylim=extrema(data[5,:]))
    title!("(w, z2) given z3=$(@sprintf("%5.2f", z3))")
end
gif(anim, "w_and_z2_give_z3.gif")

# %%
n = N
W, X, Y, Z1, Z2, Z3 = (r[1:n] for r in eachrow(data))

A = [ones(n) Z3]
@show a = A \ Z2

B = [ones(n) Z3 W]
@show b = B \ Z2

C = [ones(n) Z3 X]
@show c = C \ Z2

D = [ones(n) Z3 Z1]
@show d = D \ Z2;

# %%
n = 100
W, X, Y, Z1, Z2, Z3 = (r[1:n] for r in eachrow(data))

A = [ones(n) Z3]
@show a = A \ Z2

B = [ones(n) Z3 W]
@show b = B \ Z2

C = [ones(n) Z3 X]
@show c = C \ Z2

D = [ones(n) Z3 Z1]
@show d = D \ Z2;

# %%
res_a = [z2 - (a[1] + a[2]*z3)          for (w,x,y,z1,z2,z3) in eachcol(data)]
res_b = [z2 - (b[1] + b[2]*z3 + b[3]*w) for (w,x,y,z1,z2,z3) in eachcol(data)]
res_c = [z2 - (c[1] + c[2]*z3 + c[3]*x) for (w,x,y,z1,z2,z3) in eachcol(data)]
res_d = [z2 - (d[1] + d[2]*z3 + d[3]*z1) for (w,x,y,z1,z2,z3) in eachcol(data)]
@show std(res_a)
@show std(res_b)
@show std(res_c)
@show std(res_d)
stephist(res_a; label="z2 ~ 1 + z3", norm=true)
stephist!(res_b; label="z2 ~ 1 + z3 + w", norm=true, ls=:dash)
stephist!(res_c; label="z2 ~ 1 + z3 + x", norm=true, ls=:dashdot)
stephist!(res_d; label="z2 ~ 1 + z3 + z1", norm=true, ls=:dashdotdot)
title!("distribution of residual errors")

# %% [markdown]
# このように, 上で定義したモデル
#
# $$
# \begin{alignedat}{2}
# &
# Z_1 \sim \operatorname{Normal}(0, 1), & &
# \\ &
# Z_2 \sim \operatorname{Normal}(0, 1),
# \\ &
# Z_3 = Z_1 + Z_2 + \varepsilon_1, & & \quad \varepsilon_1 \sim \operatorname{Normal}(0, e),
# \\ &
# X = Z_1 + Z_3 + \varepsilon_2, & & \quad \varepsilon_2 \sim \operatorname{Normal}(0, e),
# \\ &
# W = X + \varepsilon_3, & & \quad \varepsilon_3 \sim \operatorname{Normal}(0, e),
# \\ &
# Y = W + Z_2 + Z_3 + \varepsilon_4, & & \quad \varepsilon_4 \sim \operatorname{Normal}(0, e) \qquad (e = 0.2)
# \end{alignedat}
# $$
#
# によって生成されたデータについて, $Z_2$ の値を $Z_3$ の値から線形回帰で予測するとき, 説明変数(独立変数)に $W$, $X$, $Z_1$ の中の1つを追加すると予測誤差が順次小さくなる.
#
# 上のモデルは次のように描ける:
#
# $$
# \begin{matrix}
# Z_1        &          &     &          & Z_2 \\
# \downarrow & \searrow &     & \swarrow & \downarrow \\
# \downarrow &          & Z_3 &          & \downarrow \\
# \downarrow & \swarrow &     & \searrow & \downarrow \\
# X          & \to      & W   & \to      & Y   \\
# \end{matrix}
# $$
#
# $Z_3$ で条件付けると, $Z_2$ と $W$ のあいだの以下のパスは閉じられる:
#
# $$
# \begin{matrix}
#            &          &       &          & Z_2 \\
#            &          &       & \swarrow &            \\
#            &          & (Z_3) &          &            \\
#            & \swarrow &       &          &            \\
# X          & \to      & W     &          & \\
# \end{matrix}
# $$
#
# しかし, $Z_3$ で条件付けると, $Z_2$ と $W$ のあいだの以下のパスは開かれる:
#
# $$
# \begin{matrix}
# Z_1        &          &       &          & Z_2 \\
# \downarrow & \searrow &       & \swarrow &            \\
# \downarrow &          & (Z_3) &          &            \\
# \downarrow &          &       &          &            \\
# X          & \to      & W     &          & \\
# \end{matrix}
# $$
#
# ゆえに $Z_3$ で条件付けると, $Z_2$ と $W$ は独立ではなくなる.
#
# 上の数値例では, $Z_2$ を説明変数として使う回帰において($Z_2$ で条件付けるとき),  $W$ も説明変数(独立変数)に追加すると予測精度が向上している. 
#
# さらに予測精度の向上は $W$ よりも $X$ を追加した場合の方が大きく, さらに $Z_1$ を追加した場合の方が大きいことも確認されている.

# %%
