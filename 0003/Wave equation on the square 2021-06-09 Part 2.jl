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
VERSION

# %% [markdown]
# Wave equaiton: $\displaystyle
# \frac{\partial^2 u}{\partial t^2} = 
# \frac{\partial^2 u}{\partial x^2} + 
# \frac{\partial^2 u}{\partial y^2}
# $

# %%
@time using DifferentialEquations
@time using Plots
@time using ProgressMeter
# Please do pkg> add https://github.com/bicycle1885/Fmt.jl
@time using Fmt

# %%
# The wave equation on the square with free boundary condition

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

x = y = range(-1, 1; length=501)
h = step(x)
p = (1/h^2,)

F(a, b, c, x, y) = exp(-((x-a)^2+(y-b)^2)/(2c^2))
f(x, y) = F(0.7, 0, 0.01, x, y) + F(-0.4, 0.6, 0.01, x, y) + F(-0.4, -0.6, 0.01, x, y)
u0 = f.(x', y)
v0 = zero(u0)
tspan = (0.0, 2.0)

#prob = DynamicalODEProblem(g!, f!, v0, u0, tspan, p)
prob = SecondOrderODEProblem(g!, v0, u0, tspan, p)
sol = nothing; GC.gc()
@time sol = solve(prob; saveat=0.01)
sol = nothing; GC.gc()
@time sol = solve(prob; saveat=0.01)

umin = minimum(minimum.(u.x[2] for u in sol.u))
umax = maximum(maximum.(u.x[2] for u in sol.u))
umin, umax

# %%
L = 201
ts = range(tspan...; length=L) |> r -> [fill(r[1], 20); r; fill(r[end], 20)]
prog = Progress(length(ts), 0)
@gif for t in ts
    z = sol(t).x[2]'
    @views heatmap(
        x, y, z;
        #x[begin:2:end], y[begin:2:end], z[begin:2:end, begin:2:end];
        clim=(-0.03, 0.03), colorbar=false, c=:bwr)
    title!(f"t = {$t:4.2f}"; titlefontsize=12, aspectratio=1)
    plot!(size=(420, 440), xlim=extrema(x), ylim=extrema(y))
    next!(prog)
end

# %%
function plot_sol(sol, x, y, t; color=:bwr)
    z = sol(t).x[2]'
    surface(
        x, y, z;
        #x[begin:2:end], y[begin:2:end], z[begin:2:end, begin:2:end];
        zlim=(-0.03, 0.03), colorbar=false, color, camera=(30, 85))
    title!(f"t = {$t:4.2f}"; titlefontsize=12)
    plot!(bg=:lightgrey, size=(640, 500))
end

plot_sol(sol, x, y, 0.6)

# %%
L = 201
ts = range(tspan...; length=L) |> r -> [fill(r[1], 20); r; fill(r[end], 20)]
prog = Progress(length(ts), 0)
@gif for t in ts
    plot_sol(sol, x, y, t)
    next!(prog)
end

# %%
function g_dirichlet!(dv, v, u, p, t)
    h⁻² = p[1]
    @. @views dv[2:end-1,2:end-1] = h⁻²*(
        u[1:end-2,2:end-1] + u[3:end,2:end-1] + u[2:end-1,1:end-2] + u[2:end-1,3:end] - 4u[2:end-1,2:end-1])
    @. @views dv[1  ,2:end-1] = h⁻²*(u[2,    2:end-1] + u[1,  1:end-2] + u[1,  3:end] - 4u[1,  2:end-1])
    @. @views dv[end,2:end-1] = h⁻²*(u[end-1,2:end-1] + u[end,1:end-2] + u[end,3:end] - 4u[end,2:end-1])
    @. @views dv[2:end-1,1  ] = h⁻²*(u[2:end-1,2    ] + u[1:end-2,1  ] + u[3:end,1  ] - 4u[2:end-1,1  ])
    @. @views dv[2:end-1,end] = h⁻²*(u[2:end-1,end-1] + u[1:end-2,end] + u[3:end,end] - 4u[2:end-1,end])
    dv[1  ,1  ] = h⁻²*(u[1  ,2    ] + u[2    ,1  ] - 4u[1  ,1  ]) 
    dv[end,1  ] = h⁻²*(u[end,2    ] + u[end-1,1  ] - 4u[end,1  ]) 
    dv[1  ,end] = h⁻²*(u[1  ,end-1] + u[2    ,end] - 4u[1  ,end]) 
    dv[end,end] = h⁻²*(u[end,end-1] + u[end-1,end] - 4u[end,end]) 
end

x = y = range(-1, 1; length=501)
h = step(x)
p = (1/h^2,)

F(a, b, c, x, y) = exp(-((x-a)^2+(y-b)^2)/(2c^2))
f(x, y) = F(0.7, 0, 0.01, x, y) + F(-0.4, 0.6, 0.01, x, y) + F(-0.4, -0.6, 0.01, x, y)
u0 = f.(x', y)
v0 = zero(u0)
tspan = (0.0, 2.0)

prob = SecondOrderODEProblem(g_dirichlet!, v0, u0, tspan, p)
sol = nothing; GC.gc()
@time sol = solve(prob; saveat=0.01)

umin = minimum(minimum.(u.x[2] for u in sol.u))
umax = maximum(maximum.(u.x[2] for u in sol.u))
umin, umax

# %%
L = 201
ts = range(tspan...; length=L) |> r -> [fill(r[1], 20); r; fill(r[end], 20)]
prog = Progress(length(ts), 0)
@gif for t in ts
    heatmap(x, y, sol(t).x[2]';
        clim=(-0.03, 0.03), colorbar=false, c=:bwr)
    title!(f"t = {$t:4.2f}"; titlefontsize=12, aspectratio=1)
    plot!(size=(420, 440), xlim=extrema(x), ylim=extrema(y))
    next!(prog)
end

# %%
plot_sol(sol, x, y, 0.6)

# %%
L = 201
ts = range(tspan...; length=L) |> r -> [fill(r[1], 20); r; fill(r[end], 20)]
prog = Progress(length(ts), 0)
@gif for t in ts
    plot_sol(sol, x, y, t)
    next!(prog)
end

# %%
