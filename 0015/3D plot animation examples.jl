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
# # 3D plot animation examples
#
# * Gen Kuroki
# * 2021-08-07～

# %% [markdown]
# ## Generate test data of simple arrays
#
# The test data `X`, `Y`, and `Z` are 1-dimensional arrays (vectors) of length $10^4+1$.

# %%
# The following code is copied from https://diffeq.sciml.ai/stable/basics/plot/#Example

using DifferentialEquations, Plots
function lorenz(du,u,p,t)
 du[1] = p[1]*(u[2]-u[1])
 du[2] = u[1]*(p[2]-u[3]) - u[2]
 du[3] = u[1]*u[2] - p[3]*u[3]
end

u0 = [1., 5., 10.]
tspan = (0., 100.)
p = (10.0,28.0,8/3)
prob = ODEProblem(lorenz, u0, tspan,p)
sol = solve(prob);

# %%
# Create data of simple arrays for test plotting

t = range(sol.prob.tspan...; length=10^4+1)
X, Y, Z = ((t -> sol(t)[i]).(t) for i in 1:3)
xlim, ylim, zlim = extrema.((X, Y, Z));

# %% [markdown]
# ## 3D plotting with gr() backend

# %%
gr(fmt = :png)
plot(X, Y, Z; label="", lw=0.5)

# %%
gr(fmt = :png)
anim = @animate for i in 1:100:length(X)
    @views plot(X[1:i], Y[1:i], Z[1:i]; label="", lw=0.5)
    plot!(; xlim, ylim, zlim)
end
gif(anim, "lorenz1.gif")

# %%
gr(fmt = :png)
x, y, z = X, Y, Z
A = plot(x, y, z; label="", lw=0.5)
B = plot(x, y; label="", lw=0.2, title="x, y", titlefontsize=8)
C = plot(x, z; label="", lw=0.2, title="x, z", titlefontsize=8)
D = plot(y, z; label="", lw=0.2, title="y, z", titlefontsize=8)
layout = @layout [
    a{0.7h}
    [b c d]
]
plot(A, B, C, D; layout, size=(640, 640))

# %%
gr(fmt = :png)
anim = @animate for i in 1:100:length(X)
    @views x, y, z = X[1:i], Y[1:i], Z[1:i]
    A = plot(x, y, z; label="", lw=0.5, xlim, ylim, zlim)
    B = plot(x, y; label="", lw=0.2, title="x, y", titlefontsize=8, xlim=xlim, ylim=ylim)
    C = plot(x, z; label="", lw=0.2, title="x, y", titlefontsize=8, xlim=xlim, ylim=zlim)
    D = plot(y, z; label="", lw=0.2, title="x, y", titlefontsize=8, xlim=ylim, ylim=zlim)
    layout = @layout [
        a{0.7h}
        [b c d]
    ]
    plot(A, B, C, D; layout, size=(640, 640))
end
gif(anim, "lorenz2.gif")

# %% [markdown]
# ## Rotated 3D plotting with pyplot() backend

# %%
pyplot(fmt = :png)
@show default(:camera)
plot(X, Y, Z; label="", lw=0.3, camera=(60, 30), size=(540, 500))

# %%
pyplot(fmt = :png)
anim = @animate for (k, i) in enumerate(1:50:length(X))
    @views x, y, z = X[1:i], Y[1:i], Z[1:i]
    plot(x, y, z; label="", lw=0.3, camera=(360k/100, 30), size=(540, 500))
    plot!(; xlim, ylim, zlim)
end
PyPlot.clf()
gif(anim, "lorenz3.gif")

# %%
pyplot(fmt = :png)
x, y, z = X, Y, Z
A = plot(x, y, z; label="", lw=0.3, xlim, ylim, zlim, camera=(360*30/100, 30))
B = plot(x, y; label="", lw=0.15, title="x, y", titlefontsize=8)
C = plot(x, z; label="", lw=0.15, title="x, z", titlefontsize=8)
D = plot(y, z; label="", lw=0.15, title="y, z", titlefontsize=8)
layout = @layout [a{0.75w} [b; c; d]]
plot(A, B, C, D; layout, size=(720, 500))

# %%
pyplot(fmt = :png)
anim = @animate for (k, i) in enumerate(1:50:length(X))
    @views x, y, z = X[1:i], Y[1:i], Z[1:i]
    A = plot(x, y, z; label="", lw=0.3, xlim, ylim, zlim, camera=(360k/100, 30))
    B = plot(x, y; label="", lw=0.15, title="x, y", titlefontsize=8, xlim=xlim, ylim=ylim)
    C = plot(x, z; label="", lw=0.15, title="x, z", titlefontsize=8, xlim=xlim, ylim=zlim)
    D = plot(y, z; label="", lw=0.15, title="y, z", titlefontsize=8, xlim=ylim, ylim=zlim)
    layout = @layout [a{0.75w} [b; c; d]]
    plot(A, B, C, D; layout, size=(720, 500))
end
PyPlot.clf()
gif(anim, "lorenz4.gif")

# %%
pyplot(fmt = :png)
anim = @animate for (k, i) in enumerate(1:50:length(X))
    @views x, y, z = X[1:i], Y[1:i], Z[1:i]
    A = plot(x, y, z; label="", lw=0.3, xlim, ylim, zlim, camera=(360k/100, 30))
    B = plot(x, y; label="", lw=0.15, title="x, y", titlefontsize=6, xlim=xlim, ylim=ylim)
    C = plot(x, z; label="", lw=0.15, title="x, z", titlefontsize=6, xlim=xlim, ylim=zlim)
    D = plot(y, z; label="", lw=0.15, title="y, z", titlefontsize=6, xlim=ylim, ylim=zlim)
    layout = @layout [a{0.75w} [b; c; d]]
    plot(A, B, C, D; layout, size=(360, 250), tickfontsize=5)
end
PyPlot.clf()
gif(anim, "lorenz5.gif")

# %% [markdown]
# ## Layout examples
#
# For details, see https://docs.juliaplots.org/latest/layouts/

# %%
using Plots
gr(fmt = :auto)

# %%
A = plot(rand(10); label="A", color=1)
B = plot(rand(10); label="B", color=2)
C = plot(rand(10); label="C", color=3)
D = plot(rand(10); label="D", color=4)
E = plot(rand(10); label="E", color=5)
F = plot(rand(10); label="F", color=6)

# Simple grid layout
layout = @layout [
    a b c
    d e f
]

plot(A, B, C, D, E, F; layout)

# %%
A = plot(rand(10); label="A", color=1)
B = plot(rand(10); label="B", color=2)
C = plot(rand(10); label="C", color=3)
D = plot(rand(10); label="D", color=4)

# Use the `_` character to ignore plots
layout = @layout [
    a b _
    c _ d
]

plot(A, B, C, D; layout)

# %%
A = plot(rand(10); label="A", color=1)
B = plot(rand(10); label="B", color=2)
C = plot(rand(10); label="C", color=3)
D = plot(rand(10); label="D", color=4)
E = plot(rand(10); label="E", color=5)
F = plot(rand(10); label="F", color=6)
G = plot(rand(10); label="G", color=7)
H = plot(rand(10); label="H", color=8)
K = plot(rand(10); label="K", color=9)

# Simple size control
layout = @layout [
    a{0.2h, 0.5w} b{0.33w} c
    d{0.7h} e f
    g h k
]

plot(A, B, C, D, E, F, G, H, K; layout)

# %%
A = plot(rand(10); label="A", color=1)
B = plot(rand(10); label="B", color=2)
C = plot(rand(10); label="C", color=3)
D = plot(rand(10); label="D", color=4)

layout = @layout [
    a{0.7h} 
    [b c d]
]

plot(A, B, C, D; layout)

# %%
A = plot(rand(10); label="A", color=1)
B = plot(rand(10); label="B", color=2)
C = plot(rand(10); label="C", color=3)
D = plot(rand(10); label="D", color=4)

layout = @layout [
    a{0.75w} [b
              c
              d]
]

plot(A, B, C, D; layout)

# %%
