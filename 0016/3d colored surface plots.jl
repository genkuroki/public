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

# %%
using Plots, StaticArrays, Parameters
default(colorbar=false)

rotate(u, v, w) = let c = cosd, s = sind
    SMatrix{3, 3}(
        c(u)*c(v), c(u)*s(v)*s(w)-s(u)*c(w), c(u)*s(v)*c(w)+s(u)*s(w),
        s(u)*c(v), s(u)*s(v)*s(w)+c(u)*c(w), c(u)*s(v)*s(w)-c(u)*s(w),
            -s(v),      c(v)*s(w),                c(v)*c(w))'
end

"""Assume `surffunc(u, v)` returns `(x, y, z, h)`."""
function plot_surface(u, v, surffunc, param; R = rotate(0, 0, 0), kwargs...)
    xyzh = surffunc.(u', v, Ref(param))
    xyz = (((x, y, z, h),) -> R * SVector(x, y, z)).(xyzh)
    x, y, z = ((a -> a[i]).(xyz) for i in 1:3)
    h = (((x, y, z, h),) -> h).(xyzh)
    surface(x, y, z; fill_z = h, kwargs...)
end

# %%
pyplot(fmt=:png)

function torusfunc(u, v, param)
    @unpack a, b, d = param
    x = (b + a * cos(u)) * cos(v)
    y = (b + a * cos(u)) * sin(v)
    z =      a * sin(u)
    h = d' * SVector(x, y, z)
    x, y, z, h
end

torusparam = (a = 5, b = 10, d = SVector(1, 0, 1))
n = 50
u_torus = v_torus = range(0, 2π; length = n + 1)
plot_surface(u_torus, v_torus, torusfunc, torusparam; size=(500, 400), 
    lims=(-20, 20), color=:CMRmap, R=rotate(30, 30, 30))

# %%
pyplot(fmt=:png)
@time anim = @animate for t in range(0, 360; length=41)[1:end-1]
    plot_surface(u_torus, v_torus, torusfunc, torusparam; size=(400, 300),
        lims=(-20, 20), color=:CMRmap, R=rotate(30, 30, 30),
        camera=(t, 20), ticks=false)
end
PyPlot.clf()
gif(anim, "torus.gif")

# %%
pyplot(fmt=:png)

function spherefunc(u, v, param)
    @unpack r, d, f = param
    x = r * cos(u) * cos(v)
    y = r * cos(u) * sin(v)
    z = r * sin(u)
    h = f(d' * SVector(x, y, z))
    x, y, z, h
end

sphereparam = (
    r = 20,
    d = SVector(0, 0, 1),
    f = h -> (1 + abs(h - 1)) * sin(h)
)
n = 50
u_sphere = range(-π/2, π/2; length = n + 1)
v_sphere = range(0, 2π; length = 2n + 1)
plot_surface(u_sphere, v_sphere, spherefunc, sphereparam; size=(500, 400), 
    lims=(-20, 20), color=:gist_earth, R=rotate(30, 30, 30),
    camera=(30, 20))

# %%
pyplot(fmt=:png)
@time anim = @animate for t in range(0, 360; length=41)[1:end-1]
    plot_surface(u_sphere, v_sphere, spherefunc, sphereparam; size=(400, 300),
        lims=(-20, 20), color=:gist_earth, R=rotate(30, 30, 30),
        camera=(t, 20), ticks=false)
end
PyPlot.clf()
gif(anim, "sphere.gif")

# %%
plotly(fmt=:auto)
plot_surface(u_torus, v_torus, torusfunc, torusparam; size=(500, 500), 
    lims=(-20, 20), color=:CMRmap, R=rotate(30, 30, 30))

# %%
plotly(fmt=:auto)
plot_surface(u_sphere, v_sphere, spherefunc, sphereparam; size=(500, 500), 
    lims=(-20, 20), color=:gist_earth, R=rotate(30, 30, 30))

# %%
using Pkg
println("Julia v", VERSION)
Pkg.status("IJulia")
Pkg.status("PyCall"; mode = PKGMODE_MANIFEST)
Pkg.status("Plots")
Pkg.status("PyPlot")
Pkg.status("Plotly")
Pkg.status("StaticArrays")
Pkg.status("Parameters")

# %%
