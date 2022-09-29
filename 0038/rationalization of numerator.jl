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
#     display_name: Julia 1.8.1
#     language: julia
#     name: julia-1.8
# ---

# %%
using Plots

f(k) = √(k+1) - √(k-1)
g(k) = 2/(√(k+1) + √(k-1))
h(k) = 1/√k

k = (10^9-1000):10^9
plot(k, f; label="√(k+1) - √(k-1)", alpha=0.5, lw=0.7)
plot!(k, g; label="2/(√(k+1) + √(k-1))", ls=:dash)
plot!(k, h; label="1/√k", ls=:dashdot)
plot!(rightmargin=5Plots.mm)

# %%
k = 10^9
@eval @show f($k) g($k) h($k);

# %%
k = big(10^9)
@eval @show f($k) g($k) h($k);

# %%
using Plots

F(a, b, c) = (-b + √(b^2 - 4a*c))/(2a)
G(a, b, c) = (2c)/(-b - √(b^2 - 4a*c))

b, c = 1, 1
a = range(-3e-8, 3e-8, 2001)
plot(a, a -> F(a, b, c); label="(-b + √(b^2 - 4a*c))/(2a)", lw=0.7)
plot!(a, a -> G(a, b, c); label="(2c)/(-b - √(b^2 - 4a*c))", ls=:dash)
plot!(rightmargin=4Plots.mm)

# %%
n = 10^9
@eval @show F(1/$n, 1, 1) G(1/$n, 1, 1);

# %%
n = big(10^9)
@eval @show F(1/$n, 1, 1) G(1/$n, 1, 1);

# %%
