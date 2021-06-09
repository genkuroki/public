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
#     display_name: Julia 1.7.0-DEV
#     language: julia
#     name: julia-1.7
# ---

# %%
@show VERSION
@time using DifferentialEquations
@time using Plots
@time using ProgressMeter

# %% [markdown]
# Wave equaiton: $\displaystyle
# \frac{\partial^2 u}{\partial t^2} = 
# \frac{\partial^2 u}{\partial x^2} + 
# \frac{\partial^2 u}{\partial y^2}
# $

# %%
# the wave equation on the square with free boundary condition

#f!(du, v, u, p, t) = @. du = v

function g!(dv, v, u, p, t)
    h⁻² = p[1]
    @. @views dv[2:end-1,2:end-1] = h⁻²*(
        u[1:end-2,2:end-1] + u[3:end,2:end-1] + u[2:end-1,1:end-2] + u[2:end-1,3:end] - 4u[2:end-1,2:end-1])
    @. @views dv[1  ,2:end-1] = h⁻²*(u[2,    2:end-1] + u[1,  1:end-2] + u[1,  3:end] - 3u[1,  2:end-1])
    @. @views dv[end,2:end-1] = h⁻²*(u[end-1,2:end-1] + u[end,1:end-2] + u[end,3:end] - 3u[end,2:end-1])
    @. @views dv[2:end-1,1  ] = h⁻²*(u[2:end-1,2    ] + u[1:end-2,1  ] + u[3:end,1  ] - 3u[2:end-1,1  ])
    @. @views dv[2:end-1,end] = h⁻²*(u[2:end-1,end-1] + u[1:end-2,end] + u[3:end,end] - 3u[2:end-1,end])
    dv[1  ,1  ] = h⁻²*(u[1  ,2    ] + u[2    ,1  ] - 2u[1  ,1  ]) 
    dv[end,1  ] = h⁻²*(u[end,2    ] + u[end-1,1  ] - 2u[end,1  ]) 
    dv[1  ,end] = h⁻²*(u[1  ,end-1] + u[2    ,end] - 2u[1  ,end]) 
    dv[end,end] = h⁻²*(u[end,end-1] + u[end-1,end] - 2u[end,end]) 
end

x = y = range(-1, 1; length=201)
h = step(x)
p = (1/h^2,)

f(x, y) = exp(-(x^2+y^2)/2)
u0 = f.(x', y)
v0 = zero(u0)
tspan = (0.0, 10.0)

#@time prob = DynamicalODEProblem(g!, f!, v0, u0, tspan, p)
@time prob = SecondOrderODEProblem(g!, v0, u0, tspan, p)
sol = nothing; GC.gc()
@time sol = solve(prob; saveat=0.02)
sol = nothing; GC.gc()
@time sol = solve(prob; sageat=0.02)
sol = nothing; GC.gc()
@time sol = solve(prob; saveat=0.02);

# %%
x = y = range(-1, 1; length=201)
h = step(x)
p = (1/h^2,)

f(x, y) = exp(-(x^2+y^2)/2)
u0 = f.(x', y)
v0 = zero(u0)
tspan = (0.0, 20.0)

#@time prob = DynamicalODEProblem(g!, f!, v0, u0, tspan, p)
@time prob = SecondOrderODEProblem(g!, v0, u0, tspan, p)
sol = nothing; GC.gc()
@time sol = solve(prob; saveat=0.02);

# %%
umin = minimum(minimum.(u.x[2] for u in sol.u))
umax = maximum(maximum.(u.x[2] for u in sol.u))
umin, umax

# %%
L = 201
ts = range(tspan...; length=L) |> r -> [fill(r[1], 20); r; fill(r[end], 20)]
prog = Progress(length(ts), 0)
@gif for t in ts
    heatmap(x, y, sol(t).x[2]'; clim=(umin, umax), colorbar=false, size=(300, 316), c=:CMRmap)
    title!("t = $(round(t, digits=1))"; titlefontsize=12, aspectratio=1)
    next!(prog)
end

# %%
L = 201
ts = range(tspan...; length=L) |> r -> [fill(r[1], 20); r; fill(r[end], 20)]
prog = Progress(length(ts), 0)
@gif for t in ts
    surface(x, y, sol(t).x[2]'; zlim=(umin, umax), colorbar=false, c=:CMRmap, camera=(30, 60))
    title!("t = $(round(t, digits=1))"; titlefontsize=12)
    next!(prog)
end

# %%
L = 201
ts = range(tspan...; length=L) |> r -> [fill(r[1], 20); r; fill(r[end], 20)]
prog = Progress(length(ts), 0)
@gif for t in ts
    wireframe(x[begin:4:end], y[begin:4:end], 
        sol(t).x[2][begin:4:end, begin:4:end]'; 
        zlim=(umin, umax), colorbar=false, c=:CMRmap, camera=(30, 60))
    title!("t = $(round(t, digits=1))"; titlefontsize=12)
    next!(prog)
end

# %%
g(x, y) = exp(-x^2/2-y^2/4)
u0 = g.(x', y)
v0 = zero(u0)
tspan = (0.0, 20.0)

@time prob = SecondOrderODEProblem(g!, v0, u0, tspan, p)
sol = nothing; GC.gc()
@time sol = solve(prob; saveat=0.02)

umax = maximum(maximum.(u.x[2] for u in sol.u))
umin = minimum(minimum.(u.x[2] for u in sol.u))

L = 201
ts = range(tspan...; length=L) |> r -> [fill(r[1], 20); r; fill(r[end], 20)]
prog = Progress(length(ts), 0)
@gif for t in ts
    surface(x, y, sol(t).x[2]'; zlim=(umin, umax), colorbar=false, c=:CMRmap, camera=(30, 60))
    title!("t = $(round(t, digits=1))"; titlefontsize=12)
    next!(prog)
end

# %%
L = 201
ts = range(tspan...; length=L) |> r -> [fill(r[1], 20); r; fill(r[end], 20)]
prog = Progress(length(ts), 0)
@gif for t in ts
    wireframe(x[begin:4:end], y[begin:4:end], 
        sol(t).x[2][begin:4:end, begin:4:end]'; 
        zlim=(umin, umax), colorbar=false, c=:CMRmap, camera=(30, 60))
    title!("t = $(round(t, digits=1))"; titlefontsize=12)
    next!(prog)
end

# %%
