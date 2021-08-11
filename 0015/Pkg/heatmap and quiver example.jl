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

# %%
using ForwardDiff
using Plots

f(x, y) = exp(-(x^2 + y^2))
u(x, y) = -ForwardDiff.derivative(Base.Fix2(f, y), x) # = f_x(x, y)
v(x, y) = -ForwardDiff.derivative(Base.Fix1(f, x), y) # = f_y(x, y)

xlim = (-1.7, 1.7)
ylim = (-1.7, 1.7)
xs = range(xlim...; length=200)
ys = range(ylim...; length=200)

c = 0.5
x = range(-1.5, 1.5; length=11)
y = range(-1.5, 1.5; length=11)
X, Y = reim(complex.(x', y)) # meshgrid
U, V = c*u.(x', y), c*v.(x', y)

heatmap(xs, ys, f)
quiver!(vec(X-U/2), vec(Y-V/2); quiver=(vec(U), vec(V)), color=:cyan)
plot!(; xlim, ylim, size=(450, 400))

# %%
using ForwardDiff
using Plots

f(x, y) = y * exp(-(x^2 + y^2))
u(x, y) = -ForwardDiff.derivative(Base.Fix2(f, y), x) # = f_x(x, y)
v(x, y) = -ForwardDiff.derivative(Base.Fix1(f, x), y) # = f_y(x, y)

xlim = (-2.2, 2.2)
ylim = (-2.2, 2.2)
xs = range(xlim...; length=200)
ys = range(ylim...; length=200)

c = 0.5
x = range(-2.0, 2.0; length=13)
y = range(-2.0, 2.0; length=13)
X, Y = reim(complex.(x', y)) # meshgrid
U, V = c*u.(x', y), c*v.(x', y)

heatmap(xs, ys, f)
quiver!(vec(X-U/2), vec(Y-V/2); quiver=(vec(U), vec(V)), color=:cyan)
plot!(; xlim, ylim, size=(450, 400))

# %%
