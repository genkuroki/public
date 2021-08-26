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
# https://discourse.julialang.org/t/solve-ode-with-many-different-initial-conditions/66384

# %%
using DifferentialEquations
using BenchmarkTools
using StaticArrays
using Base.Threads
using PyPlot

# %%
# define system parameters
const Ip = 0.5
const N  = 2
const w  = 0.057
const T  = 2π/w
const F0 = 0.075

# define laser pulse
Fx(t) = F0 * cos(w*t/(2N))^2 * cos(w*t) * (abs(t) < N*T/2)
Fy(t) = F0 * cos(w*t/(2N))^2 * sin(w*t) * (abs(t) < N*T/2)

# define rate
W_adk(F,kd) = exp(-(2.0*(kd^2+2Ip)^1.5)/3F)

# define electron trajectory
const a = 1.0
function traj(u,p,t)
    r3i = (u[1]^2+u[2]^2+a)^(-1.5)
    du1 = u[3]
    du2 = u[4]
    du3 = -u[1]*r3i-Fx(t)
    du4 = -u[2]*r3i-Fy(t)
    @SVector [du1,du2,du3,du4]
end

# define result box
const Px_min = -2.0
const Px_max =  2.0
const Px_num =  200
const Px  = LinRange(Px_min,Px_max,Px_num)
const dPx = (Px_max-Px_min)/(Px_num-1)

const Py_min = -2.0
const Py_max =  2.0
const Py_num =  200
const Py  = LinRange(Py_min,Py_max,Py_num)
const dPy = (Py_max-Py_min)/(Py_num-1)

# define simulation parameters
const Nt = 500
const Nkd = 100
const kd_max = 2.0;

# %%
# classical trajectory simulation
Prob_nthreads = let Prob_nthreads = zeros(Float64,Px_num,Py_num,nthreads())
@time @threads for tr in LinRange(-0.9N*T/2,0.9N*T/2,Nt)
    Fxtr = Fx(tr)
    Fytr = Fy(tr)
    Ftr = hypot(Fxtr,Fytr)
    phi = atan(-Fytr,-Fxtr)
    r0 = Ip/Ftr
    x0 = r0*cos(phi)
    y0 = r0*sin(phi)
    tspan = (tr,N*T/2)
    for kd in LinRange(-kd_max,kd_max,Nkd)
        kx0 = kd*cos(phi+0.5π)
        ky0 = kd*sin(phi+0.5π)
        rate = W_adk(Ftr,kd)

        u0 = @SVector [x0,y0,kx0,ky0]
        Traj = ODEProblem(traj,u0,tspan)
        sol = solve(Traj,Tsit5(),reltol=1e-6,save_everystep=false,save_start=false)
        (x,y,px,py) = sol.u[end]

        E_inf = 0.5*(px^2+py^2)-1/(x^2+y^2+a)
        E_inf >= zero(E_inf) || continue
        pxidx = Int64(round((px-Px[1])/dPx))+1
        pyidx = Int64(round((py-Py[1])/dPy))+1
        checkbounds(Bool,Prob_nthreads,pxidx,pyidx,threadid()) || continue
        @inbounds Prob_nthreads[pxidx,pyidx,threadid()] += rate
    end
end
Prob_nthreads
end
flush(stdout)

Prob = reshape(sum(Prob_nthreads,dims=3),size(Prob_nthreads)[1:2])

# plotting
fig,ax = plt.subplots()
ax.imshow(Prob',extent=(Px[1],Px[end],Py[1],Py[end]),origin="lower")
ax.set_xlim(Px[1],Px[end])
ax.set_ylim(Py[1],Py[end])
ax.set_xlabel("\$p_x\$ (a.u.)")
ax.set_ylabel("\$p_y\$ (a.u.)")
#plt.show()

# %%
Prob_nthreads = let Prob_nthreads = zeros(Float64,Px_num,Py_num,nthreads())
@time @threads for tr in LinRange(-0.9N*T/2,0.9N*T/2,Nt)
    Fxtr = Fx(tr)
    Fytr = Fy(tr)
    Ftr = hypot(Fxtr,Fytr)
    phi = atan(-Fytr,-Fxtr)
    r0 = Ip/Ftr
    x0 = r0*cos(phi)
    y0 = r0*sin(phi)
    Traj = ODEProblem(traj,zeros(Float64,4),(tr,N*T/2))
    for kd in LinRange(-kd_max,kd_max,Nkd)
        kx0 = kd*cos(phi+0.5π)
        ky0 = kd*sin(phi+0.5π)
        rate = W_adk(Ftr,kd)

        u0 = @SVector [x0,y0,kx0,ky0]
        Traj = remake(Traj;u0=u0)
        sol = solve(Traj,Tsit5(),reltol=1e-6,save_everystep=false,save_start=false)
        (x,y,px,py) = sol.u[end]

        E_inf = 0.5*(px^2+py^2)-1/(x^2+y^2+a)
        E_inf >= zero(E_inf) || continue
        pxidx = Int64(round((px-Px[1])/dPx))+1
        pyidx = Int64(round((py-Py[1])/dPy))+1
        checkbounds(Bool,Prob_nthreads,pxidx,pyidx,threadid()) || continue
        @inbounds Prob_nthreads[pxidx,pyidx,threadid()] += rate
    end
end
Prob_nthreads
end
flush(stdout)

Prob = reshape(sum(Prob_nthreads,dims=3),size(Prob_nthreads)[1:2])

# plotting
fig,ax = plt.subplots()
ax.imshow(Prob',extent=(Px[1],Px[end],Py[1],Py[end]),origin="lower")
ax.set_xlim(Px[1],Px[end])
ax.set_ylim(Py[1],Py[end])
ax.set_xlabel("\$p_x\$ (a.u.)")
ax.set_ylabel("\$p_y\$ (a.u.)")
#plt.show()

# %%
Prob_nthreads = let Prob_nthreads = zeros(Float64,Px_num,Py_num,nthreads())
@time @threads for tr in LinRange(-0.9N*T/2,0.9N*T/2,Nt)
    Fxtr = Fx(tr)
    Fytr = Fy(tr)
    Ftr = hypot(Fxtr,Fytr)
    phi = atan(-Fytr,-Fxtr)
    r0 = Ip/Ftr
    x0 = r0*cos(phi)
    y0 = r0*sin(phi)
    Traj = ODEProblem(traj,@SVector(zeros(Float64,4)),(tr,N*T/2))
    for kd in LinRange(-kd_max,kd_max,Nkd)
        kx0 = kd*cos(phi+0.5π)
        ky0 = kd*sin(phi+0.5π)
        rate = W_adk(Ftr,kd)

        u0 = @SVector [x0,y0,kx0,ky0]
        Traj = remake(Traj;u0=u0)
        sol = solve(Traj,Tsit5(),reltol=1e-6,save_everystep=false,save_start=false)
        (x,y,px,py) = sol.u[end]

        E_inf = 0.5*(px^2+py^2)-1/(x^2+y^2+a)
        E_inf >= zero(E_inf) || continue
        pxidx = Int64(round((px-Px[1])/dPx))+1
        pyidx = Int64(round((py-Py[1])/dPy))+1
        checkbounds(Bool,Prob_nthreads,pxidx,pyidx,threadid()) || continue
        @inbounds Prob_nthreads[pxidx,pyidx,threadid()] += rate
    end
end
Prob_nthreads
end
flush(stdout)

Prob = reshape(sum(Prob_nthreads,dims=3),size(Prob_nthreads)[1:2])

# plotting
fig,ax = plt.subplots()
ax.imshow(Prob',extent=(Px[1],Px[end],Py[1],Py[end]),origin="lower")
ax.set_xlim(Px[1],Px[end])
ax.set_ylim(Py[1],Py[end])
ax.set_xlabel("\$p_x\$ (a.u.)")
ax.set_ylabel("\$p_y\$ (a.u.)")
#plt.show()

# %%
Prob_nthreads = let Prob_nthreads = zeros(Float64,Px_num,Py_num,nthreads())
Trajs = Vector{ODEProblem}(undef,nthreads())
for i=1:nthreads()
    Trajs[i] = ODEProblem(traj,zeros(Float64,4),(-N*T/2,N*T/2))
end
@show typeof(Trajs)
@time @threads for tr in LinRange(-0.9N*T/2,0.9N*T/2,Nt)
    Fxtr = Fx(tr)
    Fytr = Fy(tr)
    Ftr = hypot(Fxtr,Fytr)
    phi = atan(-Fytr,-Fxtr)
    r0 = Ip/Ftr
    x0 = r0*cos(phi)
    y0 = r0*sin(phi)
    tspan = (tr,N*T/2)
    for kd in LinRange(-kd_max,kd_max,Nkd)
        kx0 = kd*cos(phi+0.5π)
        ky0 = kd*sin(phi+0.5π)
        rate = W_adk(Ftr,kd)

        u0 = @SVector [x0,y0,kx0,ky0]
        Traj = remake(Trajs[threadid()];u0=u0,tspan=tspan)
        sol = solve(Traj,Tsit5(),reltol=1e-6,save_everystep=false,save_start=false)
        (x,y,px,py) = sol.u[end]
        
        E_inf = 0.5*(px^2+py^2)-1/(x^2+y^2+a)
        E_inf >= zero(E_inf) || continue
        pxidx = Int64(round((px-Px[1])/dPx))+1
        pyidx = Int64(round((py-Py[1])/dPy))+1
        checkbounds(Bool,Prob_nthreads,pxidx,pyidx,threadid()) || continue
        @inbounds Prob_nthreads[pxidx,pyidx,threadid()] += rate
    end
end
Prob_nthreads
end
flush(stdout)

Prob = reshape(sum(Prob_nthreads,dims=3),size(Prob_nthreads)[1:2])

# plotting
fig,ax = plt.subplots()
ax.imshow(Prob',extent=(Px[1],Px[end],Py[1],Py[end]),origin="lower")
ax.set_xlim(Px[1],Px[end])
ax.set_ylim(Py[1],Py[end])
ax.set_xlabel("\$p_x\$ (a.u.)")
ax.set_ylabel("\$p_y\$ (a.u.)")
#plt.show()

# %%
Prob_nthreads = let Prob_nthreads = zeros(Float64,Px_num,Py_num,nthreads())
prob = ODEProblem(traj,@SVector(zeros(Float64,4)),(-N*T/2,N*T/2))
Trajs = Vector{typeof(prob)}(undef,nthreads())
for i=1:nthreads()
    Trajs[i] = ODEProblem(traj,@SVector(zeros(Float64,4)),(-N*T/2,N*T/2))
end
@time @threads for tr in LinRange(-0.9N*T/2,0.9N*T/2,Nt)
    Fxtr = Fx(tr)
    Fytr = Fy(tr)
    Ftr = hypot(Fxtr,Fytr)
    phi = atan(-Fytr,-Fxtr)
    r0 = Ip/Ftr
    x0 = r0*cos(phi)
    y0 = r0*sin(phi)
    tspan = (tr,N*T/2)
    for kd in LinRange(-kd_max,kd_max,Nkd)
        kx0 = kd*cos(phi+0.5π)
        ky0 = kd*sin(phi+0.5π)
        rate = W_adk(Ftr,kd)

        u0 = @SVector [x0,y0,kx0,ky0]
        Traj = remake(Trajs[threadid()];u0=u0,tspan=tspan)
        sol = solve(Traj,Tsit5(),reltol=1e-6,save_everystep=false,save_start=false)
        (x,y,px,py) = sol.u[end]
        
        E_inf = 0.5*(px^2+py^2)-1/(x^2+y^2+a)
        E_inf >= zero(E_inf) || continue
        pxidx = Int64(round((px-Px[1])/dPx))+1
        pyidx = Int64(round((py-Py[1])/dPy))+1
        checkbounds(Bool,Prob_nthreads,pxidx,pyidx,threadid()) || continue
        @inbounds Prob_nthreads[pxidx,pyidx,threadid()] += rate
    end
end
Prob_nthreads
end
flush(stdout)

Prob = reshape(sum(Prob_nthreads,dims=3),size(Prob_nthreads)[1:2])

# plotting
fig,ax = plt.subplots()
ax.imshow(Prob',extent=(Px[1],Px[end],Py[1],Py[end]),origin="lower")
ax.set_xlim(Px[1],Px[end])
ax.set_ylim(Py[1],Py[end])
ax.set_xlabel("\$p_x\$ (a.u.)")
ax.set_ylabel("\$p_y\$ (a.u.)")
#plt.show()

# %%
Prob_nthreads = let Prob_nthreads = zeros(Float64,Px_num,Py_num,nthreads())
Trajs = [ODEProblem(traj,@SVector(zeros(Float64,4)),(-N*T/2,N*T/2)) for _ in 1:nthreads()]
@show typeof(Trajs)
@time @threads for tr in LinRange(-0.9N*T/2,0.9N*T/2,Nt)
    Fxtr = Fx(tr)
    Fytr = Fy(tr)
    Ftr = hypot(Fxtr,Fytr)
    phi = atan(-Fytr,-Fxtr)
    r0 = Ip/Ftr
    x0 = r0*cos(phi)
    y0 = r0*sin(phi)
    tspan = (tr,N*T/2)
    for kd in LinRange(-kd_max,kd_max,Nkd)
        kx0 = kd*cos(phi+0.5π)
        ky0 = kd*sin(phi+0.5π)
        rate = W_adk(Ftr,kd)

        u0 = @SVector [x0,y0,kx0,ky0]
        Traj = remake(Trajs[threadid()];u0=u0,tspan=tspan)
        sol = solve(Traj,Tsit5(),reltol=1e-6,save_everystep=false,save_start=false)
        (x,y,px,py) = sol.u[end]

        E_inf = 0.5*(px^2+py^2)-1/(x^2+y^2+a)
        E_inf >= zero(E_inf) || continue
        pxidx = Int64(round((px-Px[1])/dPx))+1
        pyidx = Int64(round((py-Py[1])/dPy))+1
        checkbounds(Bool,Prob_nthreads,pxidx,pyidx,threadid()) || continue
        @inbounds Prob_nthreads[pxidx,pyidx,threadid()] += rate
    end
end
Prob_nthreads
end
flush(stdout)

Prob = reshape(sum(Prob_nthreads,dims=3),size(Prob_nthreads)[1:2])

# plotting
fig,ax = plt.subplots()
ax.imshow(Prob',extent=(Px[1],Px[end],Py[1],Py[end]),origin="lower")
ax.set_xlim(Px[1],Px[end])
ax.set_ylim(Py[1],Py[end])
ax.set_xlabel("\$p_x\$ (a.u.)")
ax.set_ylabel("\$p_y\$ (a.u.)")
#plt.show()

# %%
