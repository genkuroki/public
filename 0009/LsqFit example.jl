# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:hydrogen
#     text_representation:
#       extension: .jl
#       format_name: hydrogen
#       format_version: '1.3'
#       jupytext_version: 1.11.2
#   kernelspec:
#     display_name: Julia 1.8.0-DEV
#     language: julia
#     name: julia-1.8
# ---

# %% [markdown]
# https://twitter.com/kado_judo0312/status/1412968953450622977

# %%
using Plots
using LsqFit

# %%
model(t, p) = p[1] * exp.(-p[2] * t)
tdata = range(0, 10, length=21)

noise = [0.075383981064656, -0.8126158128658965, 1.5758753705225542, 1.4656391027715079, -0.5883715498465719, -0.680892588706559, 1.2067535098436946, -1.9517689516680197, 0.2532440785082103, 0.7943300540482112, -1.0897646372738818, -2.176470865841504, 0.9492433626731639, 1.3235403906519267, 0.13798106753050526, 1.0267333913533003, 0.03511342272679487, -1.5026374968513874, 0.5459381122556283, 0.6502966880952431, 1.455355259859152];

ydata = model(tdata, [1.0, 2.0]) .* exp.(noise);

# %%
scatter(tdata, ydata)

# %%
scatter(tdata, ydata; yscale=:log10)

# %%
p0 = [0.5, 0.5]
fit = curve_fit(model, tdata, ydata, p0)
param = fit.param
yFit(t) = model(t, param)

# %%
scatter(tdata, ydata)
plot!(tdata , yFit)

# %%
scatter(tdata, ydata; yscale=:log10)
plot!(tdata , yFit)

# %%
logmodel(t, a) = a[1] .- a[2] * t
logfit = curve_fit(logmodel, tdata, log.(ydata), [0.0, 0.0])
@show logfit.param
logyFit(t) = logmodel(t, logfit.param)

# %%
scatter(tdata, ydata; yscale=:log10)
plot!(tdata , exp∘logyFit)

# %%
scatter(tdata, ydata)
plot!(tdata , exp∘logyFit)

# %%
