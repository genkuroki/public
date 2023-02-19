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
using LinearAlgebra
dot2(x) = dot(x, x)
using Distributions
using StatsPlots
default(fmt=:png, titlefontsize=10, tickfontsize=6, guidefontsize=7)
using DataFrames
using RCall
using Random

# %% [markdown]
# データを生成するモデル:
# $$
# \begin{aligned}
# & X_1 = U_1, \quad & & U_1 \sim \operatorname{Uniform}(-\sqrt{3}, \sqrt{3}) \\
# & X_2 = X_1 + U_2, \quad & & U_2 \sim 10\operatorname{Uniform}(-\sqrt{3}, \sqrt{3}), \\
# & Y = 1 + X_1 + X_1^2 + X_2 + U_3, \quad & & U_3 \sim \operatorname{Normal}(0, 1).
# \end{aligned}
# $$
#
# 推定用のモデル1:
# $$
# Y = \beta_0 + \beta_1 X_1 + \beta_2 X_2 + \sigma U, \quad U \sim \operatorname{Normal}(0, 1).
# $$
#
# 推定用のモデル2
# $$
# Y = \beta_0 + \beta_1 X_1 + \beta_2 X_1^2 + \beta_3 X_2 + \sigma U, \quad U \sim \operatorname{Normal}(0, 1).
# $$

# %%
Random.seed!(4649373)

n = 200
x1 = rand(Uniform(-√3, √3), n)
x2 = x1 + 10rand(Uniform(-√3, √3), n)
y = @. 1 + x1 + x1^2 + x2 + randn()

anim = @animate for t in 0:3:359
    scatter(x1, x2, y; label="", title="data", msc=:auto, alpha=0.5, ms=4, camera=(t+30, 20))
    plot!(xguide="x1", yguide="x2", zguide="y")
    plot!(size=(500, 500))
end
gif(anim, "x1x2y.gif")

# %%
println("*** Y = β₀ + β₁X₁ + β₂X₂ + σU で推定 ***\n")
X = [ones(n) x1 x2]
@show β̂0, β̂1, β̂2 = X \ y
ŷ = @. β̂0 + β̂1*x1 + β̂2*x2
@show σ̂² = dot2(y - ŷ)/n
println()

X1 = [ones(n) x1]
@show α̂0, α̂1 = X1 \ y
ŷ1 = @. α̂0 + α̂1*x1
@show dot2(y - ŷ1)/n
println()

X2 = [ones(n) x2]
@show γ̂0, γ̂1 = X2 \ y
ŷ2 = @. γ̂0 + γ̂1*x2
@show dot2(y - ŷ2)/n

P1 = scatter(x1, y; label="(x1, y)", msc=:auto, alpha=0.5, ms=3)
P2 = scatter(x2, y; label="(x2, y)", msc=:auto, alpha=0.5, ms=3)

Q1 = scatter(x1, y - ŷ1; label="(x1, y - ŷ1)", msc=:auto, alpha=0.5, ms=3)
hline!([0]; label="", ls=:dot, c=:gray)
Q2 = scatter(x2, y - ŷ2; label="(x2, y - ŷ2)", msc=:auto, alpha=0.5, ms=3)
hline!([0]; label="", ls=:dot, c=:gray)

R1 = scatter(x1, y - ŷ; label="(x1, y - ŷ)", msc=:auto, alpha=0.5, ms=3)
hline!([0]; label="", ls=:dot, c=:gray)
R2 = scatter(x2, y - ŷ; label="(x2, y - ŷ)", msc=:auto, alpha=0.5, ms=3)
hline!([0]; label="", ls=:dot, c=:gray)

plot(P1, P2, Q1, Q2, R1, R2; size=(600, 800), layout=(3, 2))
plot!(legend=:outertop) |> display

anim = @animate for t in 0:3:359
    scatter(x1, x2, y - ŷ; label="", msc=:auto, alpha=0.5, ms=4, camera=(t+30, 20))
    plot!(xguide="x1", yguide="x2", zguide="y - ŷ")
    title!("residual error of model1")
    plot!(size=(500, 500))
end
gif(anim, "reserr1.gif")

# %%
println("*** Y = β₀ + β₁X₁ + β₂X₁² + β₃X₂ + σU で推定 ***\n")

X = [ones(n) x1 @.(x1^2) x2]
@show β̂0, β̂1, β̂2, β̂3 = X \ y
ŷ = @. β̂0 + β̂1*x1 + β̂2*x1^2 + β̂3*x2
@show ŝ² = dot2(y - ŷ)/(n-4)
println()

X1 = [ones(n) x1 @.(x1^2)]
@show α̂0, α̂1, α̂2 = X1 \ y
ŷ1 = @. α̂0 + α̂1*x1 + α̂2*x1^2
@show dot2(y - ŷ1)/n
println()

X2 = [ones(n) x2]
@show γ̂0, γ̂1 = X2 \ y
ŷ2 = @. γ̂0 + γ̂1*x2
@show dot2(y - ŷ2)/n

P1 = scatter(x1, y; label="(x1, y)", msc=:auto, alpha=0.5, ms=3)
P2 = scatter(x2, y; label="(x2, y)", msc=:auto, alpha=0.5, ms=3)

Q1 = scatter(x1, y - ŷ1; label="(x1, y - ŷ1)", msc=:auto, alpha=0.5, ms=3)
hline!([0]; label="", ls=:dot, c=:gray)
Q2 = scatter(x2, y - ŷ2; label="(x2, y - ŷ2)", msc=:auto, alpha=0.5, ms=3)
hline!([0]; label="", ls=:dot, c=:gray)

R1 = scatter(x1, y - ŷ; label="(x1, y - ŷ)", msc=:auto, alpha=0.5, ms=3)
hline!([0]; label="", ls=:dot, c=:gray)
R2 = scatter(x2, y - ŷ; label="(x2, y - ŷ)", msc=:auto, alpha=0.5, ms=3)
hline!([0]; label="", ls=:dot, c=:gray)

plot(P1, P2, Q1, Q2, R1, R2; size=(600, 800), layout=(3, 2))
plot!(legend=:outertop) |> display

anim = @animate for t in 0:3:359
    scatter(x1, x2, y - ŷ; label="", msc=:auto, alpha=0.5, ms=4, camera=(t+30, 20))
    plot!(xguide="x1", yguide="x2", zguide="y - ŷ")
    title!("residual error of model2")
    plot!(size=(500, 500))
end
gif(anim, "reserr2.gif")

# %%
R"""
library(dplyr)
library(ggplot2)
"""

# %%
df = DataFrame(y=y, x1=x1, x12=@.(x1^2), x2=x2)
@rput df

# %%
R"""
model1 = lm(y ~ x1 + x2, data = df)
"""

# %%
R"""
car::crPlots(model1)
""";

# %%
R"""
model2 = lm(y ~ x1 + x12 + x2, data = df)
"""

# %%
R"""
car::crPlots(model2)
""";

# %% [markdown]
# https://twitter.com/NorimitsuNishi1/status/1626514396112650240

# %%
R"""
set.seed(2349)
n <- 40
d <- tibble(x1 = runif(n, 3, 8), x2 = runif(n, 1, 5),
    epsilon = rnorm(n, 0, 0.3)) |>
    mutate(y = 10 + (4*x1) + (1 * x2^2) + epsilon)
"""
@rget d
y = @. 10 + (4*d.x1) + (1 * d.x2^2) + d.epsilon
@show y == d.y
d.x22 = @. d.x2^2
@rput d

# %%
R"""
ggplot(d) + geom_point(aes(x = x1, y = y)) +
geom_smooth(aes(x=x1, y=y), se=FALSE)
"""

# %%
R"""
ggplot(d) + geom_point(aes(x = x2, y = y)) +
geom_smooth(aes(x=x2, y=y), se=FALSE)
"""

# %% tags=[]
R"""
ggplot(d) + geom_point(aes(x = x22, y = y)) +
geom_smooth(aes(x=x22, y=y), se=FALSE)
"""

# %%

R"""
model = lm(y ~ x1, data = d)
car::crPlots(model)
"""

# %%
R"""
model = lm(y ~ x2, data = d)
car::crPlots(model)
"""

# %%
R"""
model = lm(y ~ x22, data = d)
car::crPlots(model)
"""

# %%
R"""
model = lm(y ~ x1 + x2, data = d)
car::crPlots(model)
"""

# %%
R"""
model = lm(y ~ x1 + x22, data = d)
car::crPlots(model)
"""

# %% [markdown]
# データを生成するモデル:
#
# $$
# y = 10 + 4x_1 + x_2^2 + 0.3u, \quad u \sim \operatorname{Norma}(0, 1). 
# $$
#
# モデル1:
#
# $$
# y = b_0 + b_1 x_1 + b_2 x_2 + \sigma u.
# $$
#
# モデル2:
#
# $$
# y = b_0 + b_1 x_1 + b_2 x_2^2 + \sigma u.
# $$
#
# モデル2:
#
# $$
# y = b_0 + b_1 x_1 + b_2 x_2 + b_3 x_2^2 + \sigma u.
# $$

# %%
R"""
set.seed(2349)
n <- 40
d <- tibble(x1 = runif(n, 3, 8), x2 = runif(n, 1, 5),
    epsilon = rnorm(n, 0, 0.3)) |>
    mutate(y = 10 + (4*x1) + (1 * x2^2) + epsilon)
"""

n = size(d, 1)
(; x1, x2, y) = d

anim = @animate for t in 0:3:359
    scatter(x1, x2, y; label="", title="data d", msc=:auto, alpha=0.5, ms=4, camera=(t+30, 45))
    plot!(xguide="x1", yguide="x2", zguide="y")
    plot!(size=(500, 500))
end
gif(anim, "d_x1x2y.gif")

# %%
n = size(d, 1)
(; x1, x2, y) = d

P1 = scatter(x1, y; label="(x1, y)", msc=:auto, alpha=0.8, ms=3)
P2 = scatter(x2, y; label="(x2, y)", msc=:auto, alpha=0.8, ms=3)

plot(P1, P2; size=(640, 320), legend=:outertop)

# %%
# y ~ b0 + b1*x1 + b2*x2 + u で推定

n = size(d, 1)
(; x1, x2, y) = d
X1 = [ones(n) x1 x2]

@show b1 = X1 \ y

ŷ = @. b1[1] + b1[2]*x1 + b1[3]*x2

s1 = √(dot2(y - ŷ)/(n-3))
@show s1

R1 = scatter(x1, y - ŷ; label="(x1, y - ŷ)", msc=:auto, alpha=0.8, ms=3)
hline!([0]; label="", ls=:dot, c=:gray)
R2 = scatter(x2, y - ŷ; label="(x2, y - ŷ)", msc=:auto, alpha=0.8, ms=3)
hline!([0]; label="", ls=:dot, c=:gray)

plot(R1, R2; size=(640, 320), legend=:outertop)

# %%
# y ~ b0 + b1*x1 + b2*x2^2 + u で推定

n = size(d, 1)
(; x1, x2, y) = d
X2 = [ones(n) x1 @.(x2^2)]

b_true = Float64[10, 4, 1]
@show b_true
@show b2 = X2 \ y

ŷ = @. b2[1] + b2[2]*x1 + b2[3]*x2^2

σ_true = 0.3
s2 = √(dot2(y - ŷ)/(n-3))
@show σ_true
@show s2

R1 = scatter(x1, y - ŷ; label="(x1, y - ŷ)", msc=:auto, alpha=0.8, ms=3)
hline!([0]; label="", ls=:dot, c=:gray)
R2 = scatter(x2, y - ŷ; label="(x2, y - ŷ)", msc=:auto, alpha=0.8, ms=3)
hline!([0]; label="", ls=:dot, c=:gray)

plot(R1, R2; size=(640, 320), legend=:outertop)

# %%
# y ~ b0 + b1*x1 + b2*x2 + b3*x2^2 + u で推定

n = size(d, 1)
(; x1, x2, y) = d
X3 = [ones(n) x1 x2 @.(x2^2)]

b_true = Float64[10, 4, 0, 1]
@show b_true
@show b3 = X3 \ y

ŷ = @. b3[1] + b3[2]*x1 + b3[3]*x2 + b3[4]*x2^2

σ_true = 0.3
s3 = √(dot2(y - ŷ)/(n-4))
@show σ_true
@show s3

R1 = scatter(x1, y - ŷ; label="(x1, y - ŷ)", msc=:auto, alpha=0.8, ms=3)
hline!([0]; label="", ls=:dot, c=:gray)
R2 = scatter(x2, y - ŷ; label="(x2, y - ŷ)", msc=:auto, alpha=0.8, ms=3)
hline!([0]; label="", ls=:dot, c=:gray)

plot(R1, R2; size=(640, 320), legend=:outertop)

# %%
