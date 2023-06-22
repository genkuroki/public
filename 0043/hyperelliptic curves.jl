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
#     display_name: Julia 1.9.1
#     language: julia
#     name: julia-1.9
# ---

# %%
using Plots
default(fmt=:png)
implplot(args...; kwargs...) = contour(args...; c=1, levels=[0], colorbar=false, kwargs...)

# %%
f(x) = x^2*(x+1)*(x+2)
F(x, y) = y^2 - f(x)
x = range(-5, 2, 1000)
y = range(-3, 3, 1000)
implplot(x, y, F)

# %%
f(x) = x^2*(x+1)*(x+2)*(x+3)
F(x, y) = y^2 - f(x)
x = range(-5, 2, 1000)
y = range(-3, 3, 1000)
implplot(x, y, F)

# %%
f(x) = x^2*(x+1)^2*(x+2)*(x+3)
F(x, y) = y^2 - f(x)
x = range(-5, 2, 1000)
y = range(-3, 3, 1000)
implplot(x, y, F)

# %%
F(x) = x + 2
G(x) = x * (x + 1)
F(x, y) = y^2 - F(x)*G(x)^2
x = range(-3, 2, 1000)
y = range(-1, 1, 1000)
implplot(x, y, F)

# %%
F(x) = x + 1
G(x) = x * (x + 1)
F(x, y) = y^2 - F(x)*G(x)^2
x = range(-3, 2, 1000)
y = range(-1, 1, 1000)
implplot(x, y, F)

# %%
t = range(-1.1, 1.1, 1000)
f(t) = t^2 - 1
G(x) = x * (x + 1)
g(t) = t * G(f(t))
plot(f.(t), g.(t))

# %%
ts = range(-√2.5, √2.5, 1000)
f(t) = t^2 - 2
G(x) = x * (x + 1)
g(t) = t * G(f(t))
tspan = range(extrema(ts)..., 100)
tspan = [tspan; reverse(tspan)]
anim = @animate for t in tspan
    P = plot(f.(ts), g.(ts); label="")
    scatter!([f(t)], [g(t)]; label="")
    plot!(xguide="x", yguide="y")
    title!("y²=x²(x+1)²(x+2), x=t²-2, y=tx(x+1)", titlefontsize=12)
    Q = plot(f.(ts), ts; label="")
    scatter!([f(t)], [t]; label="")
    title!("t²=x+2, t=y/(x(x+1))", titlefontsize=12)
    plot!(xguide="x", yguide="t")
    plot(P, Q; size=(720, 300))
    plot!(bottommargin=4Plots.mm)
end
gif(anim, "y²=x²(x+1)²(x+2).gif")

# %%
