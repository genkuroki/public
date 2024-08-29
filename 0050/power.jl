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
#     display_name: Julia 1.10.4
#     language: julia
#     name: julia-1.10
# ---

# %%
using Plots
default(fmt=:png, label="", guidefontsize=14)

f(a; alpha=0.05, power=0.8) = power*(1-a) / (alpha*a + power*(1-a))

# %%
plot(power -> f(0.5; power), 0, 1)
plot!(xguide="power")
plot!(xtick=0:0.1:1, ytick=0:0.1:1)
hline!([0.95, 0.9, 0.8], ls=:dot)
title!("a = 0.5,  alpha = 0.05")

# %%
plot(power -> f(0.7; power), 0, 1)
plot!(xguide="power")
plot!(xtick=0:0.1:1, ytick=0:0.1:1)
hline!([0.95, 0.9, 0.8], ls=:dot)
title!("a = 0.7,  alpha = 0.05")

# %%
f(0.7; power=0.5)

# %%
plot(a -> f(a), 0, 1)
plot!(xguide="a")
plot!(xtick=0:0.1:1, ytick=0:0.1:1)
hline!([0.95, 0.9, 0.8], ls=:dot)
title!("alpha = 0.05,  power = 0.80")

# %%
