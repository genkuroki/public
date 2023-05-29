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
#     display_name: Julia 1.9.0
#     language: julia
#     name: julia-1.9
# ---

# %%
using Plots
default(fmt=:png)

f(t) = t^2 - 1
g(t) = t*f(t)

ts = range(-√2, √2, 1000)
plot(f.(ts), g.(ts); label="")
plot!(xlim=(-2.1, 2.1), ylim=extrema(g.(ts)))

# %%
using Plots
default(fmt=:png)

xs, ys = -2.1:0.01:2.1, -1.4:0.01:1.4
contour(xs, ys, (x, y)->y^2-x^3-x^2; levels=[0.0], c=1, colorbar=false)

# %%
using Plots
default(fmt=:png)

f(t) = t^2 - 1
g(t) = t*f(t)

times = range(-√2, √2, 101)
times = [times; reverse(times)]
t = times[30]
anim = @animate for t in times
    ts = range(-√2, √2, 1000)
    P1 = plot(f.(ts), g.(ts); label="")
    plot!(xlim=(-2.1, 2.1), ylim=extrema(g.(ts)))
    plot!(xguide="x", yguide="y")
    plot!(xtick=-2:0.5:2, ytick=-2:0.5:2)
    scatter!([f(t)], [g(t)]; label="")
    title!("y² = x³ + x²   (x = t² - 1, y = tx)")
    P2 = plot(f.(ts), ts; label="")
    plot!(xlim=(-2.1, 2.1), ylim=extrema(ts))
    plot!(xguide="x", yguide="t")
    plot!(xtick=-2:0.5:2, ytick=-2:0.5:2)
    scatter!([f(t)], [t]; label="")
    title!("t² = x + 1   (t = y/x)")
    plot(P1, P2; size=(800, 300))
    plot!(titlefontsize=10)
    plot!(bottommargin=3Plots.mm)
end

gif(anim, "degellres.gif")

# %%
