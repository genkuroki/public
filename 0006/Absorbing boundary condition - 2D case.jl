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
#     display_name: Julia 1.8.0-DEV
#     language: julia
#     name: julia-1.8
# ---

# %%
@show VERSION
@time using DifferentialEquations
@time using Plots
@time using ProgressMeter

# %%
# Descriptions of discretized wave equations with various boundary conditions

function f2_dirichlet!(dv, v, u, p, t)
    f2!(dv, v, u, p, t)
    dirichlet_boundary_condition!(dv, v, u, p, t)
end

function f2_neumann!(dv, v, u, p, t)
    f2!(dv, v, u, p, t)
    neumann_boundary_condition!(dv, v, u, p, t)
end

function f2_absorbing!(dv, v, u, p, t)
    f2!(dv, v, u, p, t)
    absorbing_boundary_condition!(dv, v, u, p, t)
end

function f2!(dv, v, u, p, t)
    h = p[1]; h⁻² = 1/h^2
    @. @views dv[2:end-1,2:end-1] = h⁻²*(
        u[1:end-2,2:end-1] + u[3:end,2:end-1] + u[2:end-1,1:end-2] + u[2:end-1,3:end] - 4u[2:end-1,2:end-1]) # interior
    @. @views dv[1      ,2:end-1] = h⁻²*(u[2      ,2:end-1] + u[1      ,1:end-2] + u[1    ,3:end] - 3u[1      ,2:end-1]) # left edge
    @. @views dv[end    ,2:end-1] = h⁻²*(u[end-1  ,2:end-1] + u[end    ,1:end-2] + u[end  ,3:end] - 3u[end    ,2:end-1]) # right edge
    @. @views dv[2:end-1,1      ] = h⁻²*(u[2:end-1,2      ] + u[1:end-2,1      ] + u[3:end,1    ] - 3u[2:end-1,1      ]) # bottom edge
    @. @views dv[2:end-1,end    ] = h⁻²*(u[2:end-1,end-1  ] + u[1:end-2,end    ] + u[3:end,end  ] - 3u[2:end-1,end    ]) # top edge
    dv[1  ,1  ] = h⁻²*(u[1  ,2    ] + u[2    ,1  ] - 2u[1  ,1  ]) # bottom-left corner
    dv[end,1  ] = h⁻²*(u[end,2    ] + u[end-1,1  ] - 2u[end,1  ]) # bottom-right corner
    dv[1  ,end] = h⁻²*(u[1  ,end-1] + u[2    ,end] - 2u[1  ,end]) # top-left corner
    dv[end,end] = h⁻²*(u[end,end-1] + u[end-1,end] - 2u[end,end]) # top-right corner
end

function dirichlet_boundary_condition!(dv, v, u, p, t)
    h = p[1]; h⁻² = 1/h^2
    @. @views dv[1      ,2:end-1] += h⁻²*(-u[1      ,2:end-1]) # left edge
    @. @views dv[end    ,2:end-1] += h⁻²*(-u[end    ,2:end-1]) # right edge
    @. @views dv[2:end-1,1      ] += h⁻²*(-u[2:end-1,1      ]) # bottom edge
    @. @views dv[2:end-1,end    ] += h⁻²*(-u[2:end-1,end    ]) # top edge
    dv[1  ,1  ] += h⁻²*(-2u[1  ,1  ]) # bottom-left corner
    dv[end,1  ] += h⁻²*(-2u[end,1  ]) # bottom-right corner
    dv[1  ,end] += h⁻²*(-2u[1  ,end]) # top-left corner
    dv[end,end] += h⁻²*(-2u[end,end]) # top-right corner
end

function neumann_boundary_condition!(dv, v, u, p, t)
    h = p[1]; h⁻² = 1/h^2
    @. @views dv[1      ,2:end-1] += h⁻²*(u[2      ,2:end-1] - u[1      ,2:end-1]) # left edge
    @. @views dv[end    ,2:end-1] += h⁻²*(u[end-1  ,2:end-1] - u[end    ,2:end-1]) # right edge
    @. @views dv[2:end-1,1      ] += h⁻²*(u[2:end-1,2      ] - u[2:end-1,1      ]) # bottom edge
    @. @views dv[2:end-1,end    ] += h⁻²*(u[2:end-1,end-1  ] - u[2:end-1,end    ]) # top edge
    dv[1  ,1  ] += h⁻²*(u[1  ,2    ] + u[2    ,1  ] - 2u[1  ,1  ]) # bottom-left corner
    dv[end,1  ] += h⁻²*(u[end,2    ] + u[end-1,1  ] - 2u[end,1  ]) # bottom-right corner
    dv[1  ,end] += h⁻²*(u[1  ,end-1] + u[2    ,end] - 2u[1  ,end]) # top-left corner
    dv[end,end] += h⁻²*(u[end,end-1] + u[end-1,end] - 2u[end,end]) # top-right corner
end

function absorbing_boundary_condition!(dv, v, u, p, t)
    h = p[1]; h⁻¹ = 1/h
    neumann_boundary_condition!(dv, v, u, p, t)
    @. @views dv[1      ,2:end-1] -= h⁻¹*2v[1      ,2:end-1] # left edge
    @. @views dv[end    ,2:end-1] -= h⁻¹*2v[end    ,2:end-1] # right edge
    @. @views dv[2:end-1,1      ] -= h⁻¹*2v[2:end-1,1      ] # bottom edge
    @. @views dv[2:end-1,end    ] -= h⁻¹*2v[2:end-1,end    ] # top edge
    dv[1  ,1  ] -= h⁻¹*4v[1  ,1  ] # bottom-left corner
    dv[end,1  ] -= h⁻¹*4v[end,1  ] # bottom-right corner
    dv[1  ,end] -= h⁻¹*4v[1  ,end] # top-left corner
    dv[end,end] -= h⁻¹*4v[end,end] # top-right corner
end

# %%
function gifanim(x, y, sol; 
        umin = minimum(minimum.(u.x[2] for u in sol.u)),
        umax = maximum(maximum.(u.x[2] for u in sol.u)),
        clim = (umin, umax),
        color = :bwr, 
        bg=:lightgrey, 
        size = (300, 316)
    )
    @show umin, umax
    @show clim
    sleep(0.1)
    
    ts = range(sol.prob.tspan...; length=length(sol.t)) |> r -> [fill(r[1], 20); r; fill(r[end], 20)]
    prog = Progress(length(ts), 0)
    @gif for t in ts
        heatmap(x, y, sol(t).x[2]'; colorbar=false, clim, size, color, bg)
        title!("t = $(round(t, digits=1))"; titlefontsize=12, aspectratio=1)
        plot!(; xlim=extrema(x), ylim=extrema(y))
        next!(prog)
    end
end

function gifanim3d(x, y, sol; 
        umin = minimum(minimum.(u.x[2] for u in sol.u)),
        umax = maximum(maximum.(u.x[2] for u in sol.u)),
        clim = (umin, umax),
        color = :bwr,
        bg=:lightgrey, 
        camera = (30, 80), 
        size = (600, 500)
    )
    @show umin, umax
    @show clim
    sleep(0.1)
    
    ts = range(sol.prob.tspan...; length=length(sol.t)) |> r -> [fill(r[1], 20); r; fill(r[end], 20)]
    prog = Progress(length(ts), 0)
    @gif for t in ts
        surface(x, y, sol(t).x[2]'; colorbar=false, clim, size, color, bg, camera)
        title!("t = $(round(t, digits=1))"; titlefontsize=12)
        plot!(; xlim=extrema(x), ylim=extrema(y), zlim=clim)
        next!(prog)
    end
end

# %%
x = y = 0:0.1:π
h = step(x)
p = (; h)

f(x, y) = exp(-10(x - 2)^2 - 10(y - 2)^2)
u0 = f.(x', y)
v0 = zero(u0)
tspan = (0.0, 6.0)

prob = SecondOrderODEProblem(f2_dirichlet!, v0, u0, tspan, p)
@time sol = solve(prob; saveat=0.05)

gifanim(x, y, sol; clim=(-0.3, 0.3))

# %%
gifanim3d(x, y, sol; clim=(-0.3, 0.3))

# %%
x = y = 0:0.1:π
h = step(x)
p = (; h)

f(x, y) = exp(-10(x - 2)^2 - 10(y - 2)^2)
u0 = f.(x', y)
v0 = zero(u0)
tspan = (0.0, 6.0)

prob = SecondOrderODEProblem(f2_neumann!, v0, u0, tspan, p)
@time sol = solve(prob; saveat=0.05)

gifanim(x, y, sol; clim=(-0.3, 0.3))

# %%
gifanim3d(x, y, sol; clim=(-0.3, 0.3))

# %%
x = y = 0:0.1:π
h = step(x)
p = (; h)

f(x, y) = exp(-10(x - 2)^2 - 10(y - 2)^2)
u0 = f.(x', y)
v0 = zero(u0)
tspan = (0.0, 6.0)

prob = SecondOrderODEProblem(f2_absorbing!, v0, u0, tspan, p)
@time sol = solve(prob; saveat=0.05)

gifanim(x, y, sol; clim=(-0.3, 0.3))

# %%
gifanim3d(x, y, sol; clim=(-0.3, 0.3))

# %%
x = y = range(0, π; length=201)
h = step(x)
p = (; h)

f(x, y) = exp(-10(x - 2)^2 - 10(y - 2)^2)
u0 = f.(x', y)
v0 = zero(u0)
tspan = (0.0, 6.0)

prob = SecondOrderODEProblem(f2_dirichlet!, v0, u0, tspan, p)
@time sol = solve(prob; saveat=0.05)

gifanim(x, y, sol; clim=(-0.3, 0.3))

# %%
gifanim3d(x, y, sol; clim=(-0.3, 0.3))

# %%
x = y = range(0, π; length=201)
h = step(x)
p = (; h)

f(x, y) = exp(-10(x - 2)^2 - 10(y - 2)^2)
u0 = f.(x', y)
v0 = zero(u0)
tspan = (0.0, 6.0)

prob = SecondOrderODEProblem(f2_neumann!, v0, u0, tspan, p)
@time sol = solve(prob; saveat=0.05)

gifanim(x, y, sol; clim=(-0.3, 0.3))

# %%
gifanim3d(x, y, sol; clim=(-0.3, 0.3))

# %%
x = y = range(0, π; length=201)
h = step(x)
p = (; h)

f(x, y) = exp(-10(x - 2)^2 - 10(y - 2)^2)
u0 = f.(x', y)
v0 = zero(u0)
tspan = (0.0, 6.0)

prob = SecondOrderODEProblem(f2_absorbing!, v0, u0, tspan, p)
@time sol = solve(prob; saveat=0.05)

gifanim(x, y, sol; clim=(-0.3, 0.3))

# %%
gifanim3d(x, y, sol; clim=(-0.3, 0.3))

# %%
