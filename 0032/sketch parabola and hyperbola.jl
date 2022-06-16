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
#     display_name: Julia 1.7.3
#     language: julia
#     name: julia-1.7
# ---

# %%
using Plots
default(fmt=:png)
using StatsFuns

X(x, y) = x/y
Y(x, y) = 1 - 1/y

# %%
x = range(-10, 10, 3000) .^5
y = @. (x+2)^2 + 2
xs = range(-10, 10, 10) .^5
ys = logit.(range(0.5, 1, 10))

plot(; framestyle=:none)
for a in -10:10
    plot!(X.(a, ys), Y.(a, ys); label="", c=2, lw=0.3)
end
for b in 3:3000
    plot!(X.(xs, b), Y.(xs, b); label="", c=2, lw=0.3)
end
plot!(X.(x,y), Y.(x, y); label="", c=1)
plot!(xlim=(-1.2, 0.2), ylim = (0.4, 1.2))
title!("y = (x + 2)² + 2")

# %%
x = range(-10, 10, 3000) .^5
y = @. √((x+2)^2 + 2)
xs = range(-10, 10, 10) .^5
ys = logit.(range(0.5, 1, 10))

plot(; framestyle=:none)
for a in -10:10
    plot!(X.(a, ys), Y.(a, ys); label="", c=2, lw=0.3)
end
for b in 3:3000
    plot!(X.(xs, b), Y.(xs, b); label="", c=2, lw=0.3)
end
plot!(X.(x,y), Y.(x, y); label="", c=1)
plot!(xlim=(-2.2, 2.2), ylim = (0, 2))
title!("y = √((x + 2)² + 2)")

# %%
x = range(-10, 10, 1000) .^5
y = @. √((x+2)^2 + 2)
xs = range(-10, 10, 10) .^5
ys = logit.(range(0.5, 1, 10))

plot(; framestyle=:none)
for a in -10:10
    plot!(X.(a, ys), Y.(a, ys); label="", c=2, lw=0.3)
end
for b in 3:3000
    plot!(X.(xs, b), Y.(xs, b); label="", c=2, lw=0.3)
end
plot!(X.(x,y), Y.(x, y); label="", c=1)
plot!(-X.(x,y), 2 .- Y.(x, y); label="", c=1, ls=:dot)
plot!(xlim=(-2.2, 2.2), ylim = (0, 2))
title!("y = √((x + 2)² + 2)")

# %%
