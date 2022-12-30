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
# https://twitter.com/kimu3_slime/status/1608678692817080320?s=20&t=DN3hlMCpGLP6Q9HeHs8gaQ
#
# による以下のグラフに関する指摘について:
#
# https://gendai.media/articles/-/63904?page=3
#
# ![2022-12-30.png](attachment:63251af2-2ec0-43af-92a0-5cf04ab6868f.png)

# %% [markdown]
# このグラフはどうも以下に引用する問題の解答例のグラフっぽい。
#
# http://tomo-kumagai.eco.coocan.jp/2016_math_text_kenlo.pdf (pp.20-21)
#
# ![2022-12-30 (1).png](attachment:a70fc865-0e23-4bee-abf5-350ffc427416.png)

# %%
# See http://tomo-kumagai.eco.coocan.jp/2016_math_text_kenlo.pdf, pp.20-21

using OrdinaryDiffEq
using StaticArrays
using Plots
default(fmt=:png, titlefontsize=12)

function LotkaVolterra(u, param, t)
    S, W = u
    (; α₁, α₁, β₁, β₂) = param
    dS =  α₁*S - β₁*S*W
    dW = -α₁*W + β₂*S*W
    SVector(dS, dW)
end

param = (α₁=1.2, α₂=1.1, β₁=0.6, β₂=0.7)
S0, W0 = 2.0, 1.0
u0 = SVector(S0, W0)
tspan = (0, 25)
prob = ODEProblem(LotkaVolterra, u0, tspan, param)

# %%
dt = 0.1
sol = solve(prob, Euler(); dt)
plot(sol; label=["S(t)" "W(t)"], ls=[:solid :dash], title="Euler(), dt=$dt")

# %% [markdown]
# このグラフと件の次のグラフが一致しているように見える.
#
# ![2022-12-30.png](attachment:e42a6b10-6db8-4232-9716-abbbcb2e2d0a.png)

# %% [markdown]
# 以下は時間の刻み幅を小さくした場合.

# %%
dt = 0.01
sol = solve(prob, Euler(); dt)
plot(sol; label=["S(t)" "W(t)"], ls=[:solid :dash], title="Euler(), dt=$dt")

# %%
dt = 0.001
sol = solve(prob, Euler(); dt)
plot(sol; label=["S(t)" "W(t)"], ls=[:solid :dash], title="Euler(), dt=$dt")

# %%
sol = solve(prob, Vern7())
plot(sol; label=["S(t)" "W(t)"], ls=[:solid :dash], title="Vern7()")

# %% [markdown]
# 以下は古いバージョン

# %% tags=[]
using OrdinaryDiffEq
using StaticArrays
using Plots
default(fmt=:png, titlefontsize=12)

function LotkaVolterra_old(u, param, t)
    S, W = u
    p, q, r, s = param
    dS =  p*S - q*W*S
    dW = -r*W + s*W*S
    SVector(dS, dW)
end

param = [1.4, 0.7, 1.0, 1.0]
u0 = SVector(2.0, 1.0)
tspan = (0, 25)
prob_old = ODEProblem(LotkaVolterra_old, u0, tspan, param)

# %%
sol = solve(prob_old, Euler(); dt=0.1)
plot(sol; legend=false, ls=[:solid :dash], title="Euler(), dt=0.05")

# %%
sol = solve(prob_old, Euler(); dt=0.05)
plot(sol; legend=false, ls=[:solid :dash], title="Euler(), dt=0.05")

# %%
sol = solve(prob_old, Euler(); dt=0.01)
plot(sol; legend=false, ls=[:solid :dash], title="Euler(), dt=0.05")

# %%
sol = solve(prob_old, Vern7())
plot(sol; legend=false, ls=[:solid :dash], title="Vern7()")

# %%
