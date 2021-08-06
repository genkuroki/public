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
# Create data for test plotting

t = range(sol.prob.tspan...; length=10^4+1)
X, Y, Z = ((t -> sol(t)[i]).(t) for i in 1:3)
xlim, ylim, zlim = extrema.((X, Y, Z));

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
