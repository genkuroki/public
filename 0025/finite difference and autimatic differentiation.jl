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
#     display_name: Julia 1.7.0
#     language: julia
#     name: julia-1.7
# ---

# %% [markdown]
# https://twitter.com/mkashi/status/1469153893179596801

# %%
using Plots

f(x) = x^3 - 3x^2 + 3x - 3
df(x) = 3x^2 - 6x + 3

x = 2.1
dx = @. 2.0 ^ (-60:4)

finitediff = @. (f(x + dx) - f(x)) / dx
abserr_finirediff = @. abs(finitediff - df(x))

plot(dx, abserr_finirediff; xscale=:log10, yscale=:log10, label="")
plot!(; xtick=@.(1e1^(-20:2:2)), ytick=@.(1e1^(-20:2:2)))
plot!(; xlabel="dx", ylabel="absolute error")

# %%
df(x)

# %%
using ForwardDiff
ForwardDiff.derivative(f, x)

# %%
@show finitediff;

# %%
using FiniteDifferences
@show [central_fdm(k, 1)(f, x) for k in 2:10];

# %%
