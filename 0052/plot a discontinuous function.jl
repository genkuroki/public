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
#     display_name: Julia 1.11.2
#     language: julia
#     name: julia-1.11
# ---

# %% [markdown]
# https://x.com/dannchu/status/1880435686261813733

# %%
using Plots
default(fmt=:png)

plot(x->(2+x)/(2-x)+sqrt(x),ylim=(-10,10),xlim=(-1,10))
plot!(x->x+1)

# %%
using Plots
default(fmt=:png)
clampnan(x; lo=-Inf, hi=Inf) = (x < lo || hi < x) ? oftype(x, NaN) : x

plot(x->clampnan((2+x)/(2-x)+sqrt(x); lo=-20, hi=20), 0, 10)
plot!(x->x+1, -1, 10)
plot!(xtick=-1:10)

# %%
f(x) = mod(x, 1) == 0.5 ? oftype(x, NaN) : round(x)
P = plot(range(-3, 3, step=0.01), round; label="round(x)")
Q = plot(range(-3, 3, step=0.01), f; label="mod(x, 1) == 0.5 ? oftype(x, NaN) : round(x)", c=2)
plot(P, Q; size=(600, 700), layout=(2, 1))

# %%
