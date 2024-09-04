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
#     display_name: Julia 1.10.5
#     language: julia
#     name: julia-1.10
# ---

# %%
using ImplicitPlots
using Plots
default(fmt=:png, legendfontsize=12)

f(x, y) = x*y
g(x, y) = x^2 + 5x*y + 4y^2 - 9
xlim = (0, 4)
ylim = (0, 2)

plot()
for (t, lc, ls) in zip(0.5:0.5:1.5, 1:3, (:dash, :dashdot, :dashdotdot))
    implicit_plot!((x, y) -> f(x, y) - t; xlim, ylim, label="xy=$t", lc, ls)
end
implicit_plot!(g; xlim, ylim, label="x²+5xy+4y²=9", lc=:red, ls=:solid)
plot!(xguide="x", yguide="y")
plot!(size=(600, 300))

# %%
using ImplicitPlots
using Plots
default(fmt=:png, legendfontsize=12)

f(x, y) = x*y
g(x, y) = x^2 + 5x*y + 4y^2 - 9
xlim = (-6, 6)
ylim = (-4, 4)

plot()
for (t, lc, ls) in zip(0.5:0.5:1.5, 1:3, (:dash, :dashdot, :dashdotdot))
    implicit_plot!((x, y) -> f(x, y) - t; xlim, ylim, label="xy=$t", lc, ls)
end
implicit_plot!(g; xlim, ylim, label="x²+5xy+4y²=9", lc=:red, ls=:solid)
plot!(xguide="x", yguide="y")
plot!(size=(600, 400))

# %%
using ImplicitPlots
using Plots
default(fmt=:png, legendfontsize=12)

f(x, y) = x*y
g(x, y) = x^2 + 5x*y + 4y^2 - 9
xlim = (0, 4)
ylim = (0, 2)

@gif for t in [fill(0.5, 20); 0.5:0.01:1.0; fill(1.0, 30); 1.0:0.01:1.5; fill(1.5, 20)]
    lw = t == 1 ? 2 : 1
    plot()
    implicit_plot!((x, y) -> f(x, y) - t; xlim, ylim, label="xy=$t", lc=:blue, lw)
    implicit_plot!(g; xlim, ylim, label="x²+5xy+4y²=9", lc=:red, ls=:solid)
    plot!(xguide="x", yguide="y")
    plot!(size=(600, 300))
end

# %%
using ImplicitPlots
using Plots
default(fmt=:png, legendfontsize=12)

f(x, y) = x*y
g(x, y) = x^2 + 5x*y + 4y^2 - 9
xlim = (-6, 6)
ylim = (-4, 4)

@gif for t in [fill(0.5, 20); 0.5:0.01:1.0; fill(1.0, 30); 1.0:0.01:1.5; fill(1.5, 20)]
    lw = t == 1 ? 2 : 1
    plot()
    implicit_plot!((x, y) -> f(x, y) - t; xlim, ylim, label="xy=$t", lc=:blue, lw)
    implicit_plot!(g; xlim, ylim, label="x²+5xy+4y²=9", lc=:red, ls=:solid)
    plot!(xguide="x", yguide="y")
    plot!(size=(600, 400))
end

# %%
