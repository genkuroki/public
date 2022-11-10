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
#     display_name: Julia 1.8.2
#     language: julia
#     name: julia-1.8
# ---

# %%
using Plots

# %%
a = 2
f(a, x) = x^2 - a
df(a, x0, x) = 2x0*(x - x0) + f(a, x0)
g(a, x0) = (x0 + a/x0)/2
x0 = a
x1 = g(a, x0)
x2 = g(a, x1)
x3 = g(a, x2)
x4 = g(a, x3)
@show a, x0, x1, x2, x3, x4, √a
plot(legend=:topleft)
plot!(x -> f(a, x), 0.9√a, 1.1a; label="y = x² - $a")
plot!([x0, x0], [0, f(a, x0)]; label="", c=2)
plot!(x -> df(a, x0, x), x1, x0; label="", c=2)
plot!([x1, x1], [0, f(a, x1)]; label="", c=2)
plot!(x -> df(a, x1, x), x2, x1; label="", c=2)
plot!([x2, x2], [0, f(a, x2)]; label="", c=2)
plot!(x -> df(a, x2, x), x3, x2; label="", c=2)
hline!([0]; label="", c=:black)
plot!(size=(400, 400))

# %%
a = 3
f(a, x) = x^2 - a
df(a, x0, x) = 2x0*(x - x0) + f(a, x0)
g(a, x0) = (x0 + a/x0)/2
x0 = a
x1 = g(a, x0)
x2 = g(a, x1)
x3 = g(a, x2)
x4 = g(a, x3)
@show a, x0, x1, x2, x3, x4, √a
plot(legend=:topleft)
plot!(x -> f(a, x), 0.9√a, 1.1a; label="y = x² - $a")
plot!([x0, x0], [0, f(a, x0)]; label="", c=2)
plot!(x -> df(a, x0, x), x1, x0; label="", c=2)
plot!([x1, x1], [0, f(a, x1)]; label="", c=2)
plot!(x -> df(a, x1, x), x2, x1; label="", c=2)
plot!([x2, x2], [0, f(a, x2)]; label="", c=2)
plot!(x -> df(a, x2, x), x3, x2; label="", c=2)
hline!([0]; label="", c=:black)
plot!(size=(400, 400))

# %%
a = 4
f(a, x) = x^2 - a
df(a, x0, x) = 2x0*(x - x0) + f(a, x0)
g(a, x0) = (x0 + a/x0)/2
x0 = a
x1 = g(a, x0)
x2 = g(a, x1)
x3 = g(a, x2)
x4 = g(a, x3)
@show a, x0, x1, x2, x3, x4, √a
plot(legend=:topleft)
plot!(x -> f(a, x), 0.9√a, 1.1a; label="y = x² - $a")
plot!([x0, x0], [0, f(a, x0)]; label="", c=2)
plot!(x -> df(a, x0, x), x1, x0; label="", c=2)
plot!([x1, x1], [0, f(a, x1)]; label="", c=2)
plot!(x -> df(a, x1, x), x2, x1; label="", c=2)
plot!([x2, x2], [0, f(a, x2)]; label="", c=2)
plot!(x -> df(a, x2, x), x3, x2; label="", c=2)
hline!([0]; label="", c=:black)
plot!(size=(400, 400))

# %%
a = 5
f(a, x) = x^2 - a
df(a, x0, x) = 2x0*(x - x0) + f(a, x0)
g(a, x0) = (x0 + a/x0)/2
x0 = a
x1 = g(a, x0)
x2 = g(a, x1)
x3 = g(a, x2)
x4 = g(a, x3)
@show a, x0, x1, x2, x3, x4, √a
plot(legend=:topleft)
plot!(x -> f(a, x), 0.9√a, 1.1a; label="y = x² - $a")
plot!([x0, x0], [0, f(a, x0)]; label="", c=2)
plot!(x -> df(a, x0, x), x1, x0; label="", c=2)
plot!([x1, x1], [0, f(a, x1)]; label="", c=2)
plot!(x -> df(a, x1, x), x2, x1; label="", c=2)
plot!([x2, x2], [0, f(a, x2)]; label="", c=2)
plot!(x -> df(a, x2, x), x3, x2; label="", c=2)
hline!([0]; label="", c=:black)
plot!(size=(400, 400))

# %%
a = 2
F(a, x) = x^3 - a
dF(a, x0, x) = 3x0^2*(x - x0) + F(a, x0)
G(a, x0) = (2x0 + a/x0^2)/3
x0 = a
x1 = G(a, x0)
x2 = G(a, x1)
x3 = G(a, x2)
x4 = G(a, x3)
@show a, x0, x1, x2, x3, x4, a^(1/3)
plot(legend=:topleft)
plot!(x -> F(a, x), 0.9a^(1/3), 1.1a; label="y = x³ - $a")
plot!([x0, x0], [0, F(a, x0)]; label="", c=2)
plot!(x -> dF(a, x0, x), x1, x0; label="", c=2)
plot!([x1, x1], [0, F(a, x1)]; label="", c=2)
plot!(x -> dF(a, x1, x), x2, x1; label="", c=2)
plot!([x2, x2], [0, F(a, x2)]; label="", c=2)
plot!(x -> dF(a, x2, x), x3, x2; label="", c=2)
hline!([0]; label="", c=:black)
plot!(size=(400, 400))

# %%
a = 3
F(a, x) = x^3 - a
dF(a, x0, x) = 3x0^2*(x - x0) + F(a, x0)
G(a, x0) = (2x0 + a/x0^2)/3
x0 = a
x1 = G(a, x0)
x2 = G(a, x1)
x3 = G(a, x2)
x4 = G(a, x3)
@show a, x0, x1, x2, x3, x4, a^(1/3)
plot(legend=:topleft)
plot!(x -> F(a, x), 0.9a^(1/3), 1.1a; label="y = x³ - $a")
plot!([x0, x0], [0, F(a, x0)]; label="", c=2)
plot!(x -> dF(a, x0, x), x1, x0; label="", c=2)
plot!([x1, x1], [0, F(a, x1)]; label="", c=2)
plot!(x -> dF(a, x1, x), x2, x1; label="", c=2)
plot!([x2, x2], [0, F(a, x2)]; label="", c=2)
plot!(x -> dF(a, x2, x), x3, x2; label="", c=2)
hline!([0]; label="", c=:black)
plot!(size=(400, 400))

# %%
a = 1.5^3
F(a, x) = x^3 - a
dF(a, x0, x) = 3x0^2*(x - x0) + F(a, x0)
G(a, x0) = (2x0 + a/x0^2)/3
x0 = a
x1 = G(a, x0)
x2 = G(a, x1)
x3 = G(a, x2)
x4 = G(a, x3)
@show a, x0, x1, x2, x3, x4, a^(1/3)
plot(legend=:topleft)
plot!(x -> F(a, x), 0.9a^(1/3), 1.1a; label="y = x³ - $a")
plot!([x0, x0], [0, F(a, x0)]; label="", c=2)
plot!(x -> dF(a, x0, x), x1, x0; label="", c=2)
plot!([x1, x1], [0, F(a, x1)]; label="", c=2)
plot!(x -> dF(a, x1, x), x2, x1; label="", c=2)
plot!([x2, x2], [0, F(a, x2)]; label="", c=2)
plot!(x -> dF(a, x2, x), x3, x2; label="", c=2)
hline!([0]; label="", c=:black)
plot!(size=(400, 400))

# %%
