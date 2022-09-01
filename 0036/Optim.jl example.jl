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
#     display_name: Julia 1.8.0
#     language: julia
#     name: julia-1.8
# ---

# %% [markdown]
# https://twitter.com/dannchu/status/1563401132902588417
#
# ![FbJQ9oxacAAPvN2.png](attachment:120d12de-b14e-4749-a52f-09160a479bd4.png)

# %%
using Optim

function F(x; w = (0.39, 0.30, 0.21, 0.10))
    p, q, r = x
    a, b, c, d = w
    (p + q + r - 1)^2 + (p^2 + 2p*r - a)^2 + (r^2 - b)^2 + (q^2 + 2q*r - c)^2 + (2p*q - d)^2
end

o = optimize(F, fill(1/3, 3), LBFGS())

# %%
o.minimizer

# %%
sum(o.minimizer)

# %%
using Optim

function G(x; w = (0.39, 0.30, 0.21, 0.10))
    p, q  = x
    r = 1 - p - q
    a, b, c, d = w
    (p^2 + 2p*r - a)^2 + (r^2 - b)^2 + (q^2 + 2q*r - c)^2 + (2p*q - d)^2
end

o = optimize(G, fill(1/3, 2), LBFGS())

# %%
sol = [o.minimizer; 1 - sum(o.minimizer)]

# %%
using Optim
using Plots

function G(x; w = (0.39, 0.30, 0.21, 0.10))
    p, q  = x
    r = 1 - p - q
    a, b, c, d = w
    (p^2 + 2p*r - a)^2 + (r^2 - b)^2 + (q^2 + 2q*r - c)^2 + (2p*q - d)^2
end
o = optimize(G, fill(1/3, 2), LBFGS())
@show minimum(o)
sol = [o.minimizer; 1 - sum(o.minimizer)]
@show sol

p = range(0, 1, 200)
q = range(0, 1, 200)
heatmap(p, q, (p, q) -> p + q ≤ 1 ? √G((p, q)) : NaN; clim=(0, 1))
scatter!([sol[1]], [sol[2]]; label="(p, q)", c=:cyan)
plot!(xlim=(0, 1), ylim=(0, 1), xtick=0:0.1:1, ytick=0:0.1:1)
plot!(xguide="p", yguide="q", size=(480, 400))

# %%
using Optim

function H(x; w = (0.39, 0.30))
    p, q, c, d  = x
    r = 1 - p - q
    a, b = w
    (p^2 + 2p*r - a)^2 + (r^2 - b)^2 + (q^2 + 2q*r - c)^2 + (2p*q - d)^2
end
o = optimize(H, fill(1/4, 4), LBFGS())
@show minimum(o)
@show o.minimizer;

# %%
