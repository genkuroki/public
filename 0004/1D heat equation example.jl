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
using Printf
using Plots
default(size = (400, 250), legend=false, titlefontsize=12)

# %%
function laplacian!(v::Vector{Float64}, u::Vector{Float64}, dx::Float64)::Nothing
    v[1] = (u[end] + u[2] - 2u[1])/dx^2
    v[2:end-1] = (u[1:end-2] + u[3:end] - 2u[2:end-1])/dx^2
    v[end] = (u[end-1] + u[1] - 2u[1end])/dx^2
    return
end

# %%
n = 200
x = range(-π, π; length=n+1)[1:end-1]
dx = step(x)
f(x) = sin(x) - cos(2x)
u = f.(x)
plot(x, u; label="", title="y = u(x)")

# %%
v = similar(u)
laplacian!(v, u, dx)
plot(x, v; label="", title="y = u''(x)")

# %%
function laplacian!(v::Vector{Float64}, u::Vector{Float64}, dx::Float64)::Nothing
    v[1] = (u[end] + u[2] - 2u[1])/dx^2
    v[2:end-1] = (u[1:end-2] + u[3:end] - 2u[2:end-1])/dx^2
    v[end] = (u[end-1] + u[1] - 2u[1end])/dx^2
    return
end

function heateq(u0::Vector{Float64}, dx::Float64, tmax::Float64, N::Int=200)::Tuple{Vector{Float64}, Matrix{Float64}}
    t = 0:dx:tmax
    dt = t[2] - t[1]
    u = Matrix{Float64}(undef, length(u0), length(t)+1)
    u[:, 1] = u0
    v = Vector{Float64}(undef, length(u0))
    for i in 2:length(t)+1
        u[:, i] = u[:, i-1]
        for _ in 1:N
            laplacian!(v, u[:, i], dx)
            u[:, i] += v*dt/N
        end
    end
    t, u
end

n = 200
x = range(-π, π; length=n+1)[1:end-1]
dx = step(x)
f(x) = sin(x) - cos(2x)
u0 = f.(x)
tmax = 1.0

@time t, u = heateq(u0, dx, tmax)
@time t, u = heateq(u0, dx, tmax)
@gif for i in 1:length(t)
    title = @sprintf("t = %4.2f", t[i])
    plot(x, u[:, i]; ylim=extrema(u0), title)
end

# %%
function laplacian2!(v::Vector{Float64}, u::Vector{Float64}, dx::Float64)::Nothing
    v[1] = (u[end] + u[2] - 2u[1])/dx^2
    @. @views v[2:end-1] = (u[1:end-2] + u[3:end] - 2u[2:end-1])/dx^2
    v[end] = (u[end-1] + u[1] - 2u[1end])/dx^2
    return
end

function heateq2(u0::Vector{Float64}, dx::Float64, tmax::Float64, N::Int=200)::Tuple{Vector{Float64}, Matrix{Float64}}
    t = 0:dx:tmax
    dt = t[2] - t[1]
    u = Matrix{Float64}(undef, length(u0), length(t)+1)
    u[:, 1] = u0
    v = Vector{Float64}(undef, length(u0))
    for i in 2:length(t)+1
        @. @views u[:, i] = u[:, i-1]
        for _ in 1:N
            @views laplacian2!(v, u[:, i], dx)
            @. @views u[:, i] += v*dt/N
        end
    end
    t, u
end

n = 200
x = range(-π, π; length=n+1)[1:end-1]
dx = step(x)
f(x) = sin(x) - cos(2x)
u0 = f.(x)
tmax = 1.0

t, u = heateq2(u0, dx, tmax)

# %%
function laplacian3!(v::Vector{Float64}, u, dx::Float64)::Nothing
    v[1] = (u[end] + u[2] - 2u[1])/dx^2
    @. @views v[2:end-1] = (u[1:end-2] + u[3:end] - 2u[2:end-1])/dx^2
    v[end] = (u[end-1] + u[1] - 2u[1end])/dx^2
    return
end

function heateq3(u0::Vector{Float64}, dx::Float64, tmax::Float64, N::Int=200)::Tuple{Vector{Float64}, Matrix{Float64}}
    t = 0:dx:tmax
    dt = t[2] - t[1]
    u = Matrix{Float64}(undef, length(u0), length(t)+1)
    u[:, 1] = u0
    v = Vector{Float64}(undef, length(u0))
    N = 100
    for i in 2:length(t)+1
        @. @views u[:, i] = u[:, i-1]
        for _ in 1:N
            @views laplacian3!(v, u[:, i], dx)
            @. @views u[:, i] += v*dt/N
        end
    end
    t, u
end

n = 200
x = range(-π, π; length=n+1)[1:end-1]
dx = step(x)
f(x) = sin(x) - cos(2x)
u0 = f.(x)
tmax = 1.0

@time t, u = heateq3(u0, dx, tmax)
@time t, u = heateq3(u0, dx, tmax)
@gif for i in 1:length(t)
    title = @sprintf("t = %4.2f", t[i])
    plot(x, u[:, i]; ylim=extrema(u0), title)
end

# %%
function laplacian4!(v, u, dx)
    v[begin] = (u[end] + u[begin+1] - 2u[begin])/dx^2
    @. @views v[begin+1:end-1] = (u[begin:end-2] + u[begin+2:end] - 2u[begin+1:end-1])/dx^2
    v[end] = (u[end-1] + u[begin] - 2u[1end])/dx^2
    return
end

function heateq4(u0, dx, tmax, N=200)
    t = 0:dx:tmax
    dt = step(t)
    u = similar(u0, length(u0), length(t)+1)
    u[:, 1] = u0
    v = similar(u0)
    N = 100
    for i in 2:length(t)+1
        @. @views u[:, i] = u[:, i-1]
        for _ in 1:N
            @views laplacian4!(v, u[:, i], dx)
            @. @views u[:, i] += v*dt/N
        end
    end
    t, u
end


n = 200
x = range(-π, π; length=n+1)[1:end-1]
dx = step(x)
f(x) = sin(x) - cos(2x)
u0 = f.(x)
tmax = 1.0

@time t, u = heateq4(u0, dx, tmax)
@time t, u = heateq4(u0, dx, tmax)
@gif for i in 1:length(t)
    title = @sprintf("t = %4.2f", t[i])
    plot(x, u[:, i]; ylim=extrema(u0), title)
end

# %%
n = 200
x = range(-Float32(π), π; length=n+1)[1:end-1]
dx = step(x)
@show typeof(x) typeof(dx)
f(x) = sin(x) - cos(2x)
u0 = f.(x)
@show typeof(u0)
tmax = Float32(1)

t, u = heateq4(u0, dx, tmax)
@show typeof(t) typeof(u)
@gif for i in 1:length(t)
    title = @sprintf("t = %4.2f", t[i])
    plot(x, u[:, i]; ylim=extrema(u0), title)
end

# %%
using StaticArrays

n = 200
x = range(-π, π; length=n+1)[1:end-1]
dx = step(x)
m = 200
noise1, noise2 = 0.3randn(m), 0.3randn(m)
f(x) = SVector{m}((1 .+ noise1)*sin(x) .- (1 .+ noise2)*cos(2x))
u0 = f.(x)

plot()
for k in 1:m
    plot!(x, (p -> p[k]).(u0); lw=0.1, color=:blue)
end
plot!()

# %%
tmax = 1.0
@time t, u = heateq4(u0, dx, tmax)
@time t, u = heateq4(u0, dx, tmax)

ylim = (minimum(minimum.(u0)), maximum(maximum.(u0)))
@gif for i in 1:length(t)
    title = @sprintf("t = %4.2f", t[i])
    plot()
    for k in 1:m
        plot!(x, (p -> p[k]).(u[:, i]); ylim, lw=0.1, color=:blue)
    end
    title!(title)
end

# %%
