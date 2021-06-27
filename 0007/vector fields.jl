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
VERSION

# %%
using Plots
using Zygote

# %%
meshgrid(x, y) = reim(complex.(x', y))

x = 1:4
y = 10(1:3)
X, Y = meshgrid(x, y)
display(X)
display(Y)

# %%
function vf!(x, y, f; scale=1, kwargs...)
    X, Y = meshgrid(x, y)
    u(x, y) = scale * f(x, y)[1]
    v(x, y) = scale * f(x, y)[2]
    U = u.(X, Y)
    V = v.(X, Y)
    X -= U/2
    Y -= V/2
    quiver!(vec(X), vec(Y); quiver = (vec(U), vec(V)), kwargs...)
end

vf(x, y, f; kwargs...) = (plot(); vf!(x, y, f; kwargs...))

# %%
g(x, y) = (-y, x)
x = y = range(-2, 2, length=21)[2:2:end]
vf(x, y, g; scale=0.2, size=(400, 400))
plot!(xlim=(-2, 2), ylim=(-2, 2))

# %%
f(x, y) = x^3 - 3x + y^2
df(x, y) = gradient(f, x, y)
xs = ys = range(-2, 2, length=101)
heatmap(xs, ys, f; color=:rainbow)
x = y = range(-2, 2, length=21)[2:2:end]
vf!(x, y, df; scale=0.05, color=:white)
plot!(xlim=extrema(xs), ylim=extrema(ys))

# %%
