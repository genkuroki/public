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
#     display_name: Julia 1.6.2
#     language: julia
#     name: julia-1.6
# ---

# %% [markdown]
# # Examples of marginal plots
#
# Extended version of https://discourse.julialang.org/t/multiple-plots-some-refinements/64800/5

# %%
using Plots
using Distributions
using KernelDensity
using LinearAlgebra

# true distribution of sample
dist_true = MvNormal([1.0, 2.0], [2.0 -1.0; -1.0 4.0])
mu_x, mu_y = mean(dist_true)
s_x, s_y = .√diag(cov(dist_true))
distx_true, disty_true = Normal(mu_x, s_x), Normal(mu_y, s_y)

# generate test sample
n = 2^10
sample = rand(dist_true, n)

# kernel density estimation
X, Y = sample[1, :], sample[2, :]
k, kx, ky = kde.(((X, Y), X, Y))
ik, ikx, iky = InterpKDE.((k, kx, ky));

# %% [markdown]
# ## Sample marginal plots

# %%
layout = @layout [
    a             _
    b{0.8w, 0.8h} c
]

xlim, ylim = extrema.((k.x, k.y))
legend = colorbar = false
a = plot(kx.x, kx.density; xlim, legend)
c = plot(ky.density, ky.x; ylim, legend, xrotation=90)
b = contour(k.x, k.y, k.density; xlim, ylim, legend)

plot(a, b, c; layout, link=:both, size=(500, 500))

# %%
layout = @layout [
    a             _
    b{0.8w, 0.8h} c
]

xlim, ylim = extrema.((k.x, k.y))
legend = colorbar = false
a = plot(kx.x, kx.density; xlim, legend)
c = plot(ky.density, ky.x; ylim, legend, xrotation=90)
b = contourf(k.x, k.y, k.density; xlim, ylim, legend)

plot(a, b, c; layout, link=:both, size=(500, 500))

# %%
layout = @layout [
    a             _
    b{0.8w, 0.8h} c
]

xlim, ylim = extrema.((k.x, k.y))
legend = colorbar = false
a = plot(kx.x, kx.density; xlim, legend)
c = plot(ky.density, ky.x; ylim, legend, xrotation=90)
b = heatmap(k.x, k.y, k.density; xlim, ylim, legend)

plot(a, b, c; layout, link=:both, size=(500, 500))

# %%
layout = @layout [
    a             _
    b{0.8w, 0.8h} c
]

xlim, ylim = extrema.((k.x, k.y))
legend = colorbar = false
a = plot(kx.x, kx.density; xlim, legend)
c = plot(ky.density, ky.x; ylim, legend, xrotation=90)
b = scatter(X, Y; xlim, ylim, legend, marker_z=pdf.(Ref(ik), X, Y), alpha=0.7, msw=0)

plot(a, b, c; layout, link=:both, size=(500, 500))

# %%
layout = @layout [
    a             _
    b{0.8w, 0.8h} c
]

xlim, ylim = extrema.((k.x, k.y))
legend = colorbar = false
a = plot(kx.x, kx.density; xlim, legend)
c = plot(ky.density, ky.x; ylim, legend, xrotation=90)
b = scatter(X, Y; xlim, ylim, legend, marker_z=pdf.(Ref(ik), X, Y), alpha=0.7, msw=0, color=:rainbow)

plot(a, b, c; layout, link=:both, size=(500, 500))

# %% [markdown]
# ## True distribution marginal plots

# %%
layout = @layout [
    a             _
    b{0.8w, 0.8h} c
]

xlim, ylim = extrema.((k.x, k.y))
legend = colorbar = false
a = plot(kx.x, pdf.(distx_true, kx.x); xlim, label="true")
c = plot(pdf.(disty_true, ky.x), ky.x; ylim, xrotation=90, label="true")
b = contour(k.x, k.y, (x, y)->pdf(dist_true, [x, y]); xlim, ylim, legend)

plot(a, b, c; layout, link=:both, size=(500, 500))

# %%
layout = @layout [
    a             _
    b{0.8w, 0.8h} c
]

xlim, ylim = extrema.((k.x, k.y))
legend = colorbar = false
a = plot(kx.x, pdf.(distx_true, kx.x); xlim, label="true")
c = plot(pdf.(disty_true, ky.x), ky.x; ylim, xrotation=90, label="true")
b = contourf(k.x, k.y, (x, y)->pdf(dist_true, [x, y]); xlim, ylim, legend)

plot(a, b, c; layout, link=:both, size=(500, 500))

# %%
layout = @layout [
    a             _
    b{0.8w, 0.8h} c
]

xlim, ylim = extrema.((k.x, k.y))
legend = colorbar = false
a = plot(kx.x, pdf.(distx_true, kx.x); xlim, label="true")
c = plot(pdf.(disty_true, ky.x), ky.x; ylim, xrotation=90, label="true")
b = heatmap(k.x, k.y, (x, y)->pdf(dist_true, [x, y]); xlim, ylim, legend)

plot(a, b, c; layout, link=:both, size=(500, 500))

# %% [markdown]
# ## Other examples

# %%
layout = @layout [
    a             _
    b{0.8w, 0.8h} c
]

xlim, ylim = extrema.((k.x, k.y))
legend = colorbar = false

a = plot(kx.x, kx.density; xlim, label="X")
plot!(kx.x, pdf.(distx_true, kx.x); label="true", ls=:dash)

c = plot(ky.density, ky.x; ylim, xrotation=90, label="Y")
plot!(pdf.(disty_true, ky.x), ky.x; label="true", ls=:dash)

b = contour(k.x, k.y, k.density; xlim, ylim, legend)

plot(a, b, c; layout, link=:both, size=(500, 500))

# %%
layout = @layout [
    a             _
    b{0.8w, 0.8h} c
]

xlim, ylim = extrema.((k.x, k.y))
legend = colorbar = false

a = plot(kx.x, kx.density; xlim, label="X")
plot!(kx.x, pdf.(distx_true, kx.x); label="true", ls=:dash)

c = plot(ky.density, ky.x; ylim, xrotation=90, label="Y")
plot!(pdf.(disty_true, ky.x), ky.x; label="true", ls=:dash)

b = contourf(k.x, k.y, k.density; xlim, ylim, legend)

plot(a, b, c; layout, link=:both, size=(500, 500))

# %%
layout = @layout [
    a             _
    b{0.8w, 0.8h} c
]

xlim, ylim = extrema.((k.x, k.y))
legend = colorbar = false

a = plot(kx.x, kx.density; xlim, label="X")
plot!(kx.x, pdf.(distx_true, kx.x); label="true", ls=:dash)

c = plot(ky.density, ky.x; ylim, xrotation=90, label="Y")
plot!(pdf.(disty_true, ky.x), ky.x; label="true", ls=:dash)

b = heatmap(k.x, k.y, k.density; xlim, ylim, legend)

plot(a, b, c; layout, link=:both, size=(500, 500))

# %%
layout = @layout [
    a             _
    b{0.8w, 0.8h} c
]

xlim, ylim = extrema.((k.x, k.y))
legend = colorbar = false

a = plot(kx.x, kx.density; xlim, label="X")
plot!(kx.x, pdf.(distx_true, kx.x); label="true", ls=:dash)

c = plot(ky.density, ky.x; ylim, xrotation=90, label="Y")
plot!(pdf.(disty_true, ky.x), ky.x; label="true", ls=:dash)

b = scatter(X, Y; xlim, ylim, colorbar, marker_z=pdf.(Ref(ik), X, Y), alpha=0.7, msw=0, label="(X, Y)")
contour!(k.x, k.y, (x, y)->pdf(dist_true, [x, y]); xlim, ylim, legend, ls=:dash, color=2)

plot(a, b, c; layout, link=:both, size=(500, 500))

# %%
layout = @layout [
    a             _
    b{0.8w, 0.8h} c
]

xlim, ylim = extrema.((k.x, k.y))
legend = colorbar = false

a = plot(kx.x, kx.density; xlim, label="X")
plot!(kx.x, pdf.(distx_true, kx.x); label="true", ls=:dash)

c = plot(ky.density, ky.x; ylim, xrotation=90, label="Y")
plot!(pdf.(disty_true, ky.x), ky.x; label="true", ls=:dash)

b = scatter(X, Y; xlim, ylim, colorbar, marker_z=pdf.(Ref(ik), X, Y), alpha=0.7, msw=0, color=:rainbow, label="(X, Y)")
contour!(k.x, k.y, (x, y)->pdf(dist_true, [x, y]); xlim, ylim, legend, ls=:dash, color=2)

plot(a, b, c; layout, link=:both, size=(500, 500))

# %%
