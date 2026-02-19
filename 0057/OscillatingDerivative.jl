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
#     display_name: Julia
#     language: julia
#     name: julia
# ---

# %% [markdown]
# https://colab.research.google.com/github/genkuroki/public/blob/main/0057/OscillatingDerivative.ipynb

# %%
#haskey(ENV, "COLAB_GPU") && (import Pkg; Pkg.add("ForwardDiff"))
#using ForwardDiff
using Plots
default(fmt=:png, legend=false)

f(x) = if x == 0
    zero(float(x))
else
    exp(-1/x^2) * (2 + sin(1/x^4))
end

logf(x) = if x == 0
    zero(float(x))
else
    log(2 + sin(1/x^4)) - 1/x^2
end

# dlogf(x) = if x == 0
#     zero(float(x))
# else
#     ForwardDiff.derivative(logf, x)
# end

dlogf(x) = if x == 0
    zero(float(x))
else
    -4/x^5 * cos(1/x^4) / (2 + sin(1/x^4)) + 2/x^3
end

# %%
xs = range(-0.6, 0.6, 10001)
P1 = plot(xs, f; title=raw"$f(x) = \mathrm{if}\ x \ne 0 \ \mathrm{then}\ \exp(-1/x^2)(2 + \sin(1/x^4)) \ \mathrm{else}\ 0$")
P2 = plot(xs, x -> abs(x)^5 * dlogf(x); title=raw"$|x|^5(\log f(x))'$")
plot(P1, P2; size=(1200, 400))

# %%
