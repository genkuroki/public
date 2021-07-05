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
using LinearAlgebra
dot2(v) = dot(v, v)
using Plots

n = 9
z = exp(2π*im/n)
t = range(0, 2π; length=101)
v0 = normalize([z^k for k in 0:n-1])
#v0 = normalize!(rand(ComplexF64, n))
v(t) = exp(im*t)*v0
f(t) = real(dot2(real(v(t))))
g(t) = real(dot2(imag(v(t))))
plot(t, f.(t); label="f(t)", ylim=(0, 1))
plot!(t, g.(t); label="g(t)", ls=:auto)

# %%
using LinearAlgebra
dot2(v) = dot(v, v)
using Plots

n = 9
z = exp(2π*im/n)
t = range(0, 2π; length=101)
#v0 = normalize([z^k for k in 0:n-1])
v0 = normalize!(rand(ComplexF64, n))
v(t) = exp(im*t)*v0
f(t) = real(dot2(real(v(t))))
g(t) = real(dot2(imag(v(t))))
plot(t, f.(t); ylim=(0, 1))
plot!(t, g.(t); ls=:auto)

# %%
using LinearAlgebra
using Plots

n = 10
v0 = normalize(randn(ComplexF64, n))
v(t) = exp(im*t)*v0
f(t) = real(dot2(real(v(t))))
g(t) = real(dot2(imag(v(t))))

s0 = sum(z^2 for z in v0)
n0 = norm(v0)^2
F(t) = (n0 + real(exp(im*2t) * s0))/2
G(t) = (n0 - real(exp(im*2t) * s0))/2

t = range(0, 2π; length=201)
P = plot(t, f.(t); label="f(t)", lw=2)
plot!(t, F.(t); label="F(t)", ls=:auto, lw=2, ylim=(0, 1))
plot!(t, g.(t); label="g(t)", lw=2)
plot!(t, G.(t); label="G(t)", ls=:auto, lw=2, ylim=(0, 1))
plot!(; legend=:outertopright)

# %%
using LinearAlgebra
using Plots

n = 10
v0 = normalize(randn(ComplexF64, n))
v(t) = exp(im*t)*v0
f(t) = real(dot2(real(v(t))))
g(t) = real(dot2(imag(v(t))))

s0 = sum(z^2 for z in v0)
n0 = norm(v0)^2
F(t) = (n0 + real(exp(im*2t) * s0))/2
G(t) = (n0 - real(exp(im*2t) * s0))/2

t = range(0, 2π; length=201)
P = plot(t, f.(t); label="f(t)", lw=2)
plot!(t, F.(t); label="F(t)", ls=:auto, lw=2, ylim=(0, 1))
plot!(t, g.(t); label="g(t)", lw=2)
plot!(t, G.(t); label="G(t)", ls=:auto, lw=2, ylim=(0, 1))
plot!(; legend=:outertopright)

# %%
using LinearAlgebra
using Plots

n = 10
v0 = normalize(randn(ComplexF64, n))
v(t) = exp(im*t)*v0
f(t) = real(dot2(real(v(t))))
g(t) = real(dot2(imag(v(t))))

s0 = sum(z^2 for z in v0)
n0 = norm(v0)^2
F(t) = (n0 + real(exp(im*2t) * s0))/2
G(t) = (n0 - real(exp(im*2t) * s0))/2

t = range(0, 2π; length=201)
P = plot(t, f.(t); label="f(t)", lw=2)
plot!(t, F.(t); label="F(t)", ls=:auto, lw=2, ylim=(0, 1))
plot!(t, g.(t); label="g(t)", lw=2)
plot!(t, G.(t); label="G(t)", ls=:auto, lw=2, ylim=(0, 1))
plot!(; legend=:outertopright)

# %%
using LinearAlgebra
using Plots

n = 10
v0 = normalize(randn(ComplexF64, n))
v(t) = exp(im*t)*v0
f(t) = real(dot2(real(v(t))))
g(t) = real(dot2(imag(v(t))))

s0 = sum(z^2 for z in v0)
n0 = norm(v0)^2
F(t) = (n0 + real(exp(im*2t) * s0))/2
G(t) = (n0 - real(exp(im*2t) * s0))/2

t = range(0, 2π; length=201)
P = plot(t, f.(t); label="f(t)", lw=2)
plot!(t, F.(t); label="F(t)", ls=:auto, lw=2, ylim=(0, 1))
plot!(t, g.(t); label="g(t)", lw=2)
plot!(t, G.(t); label="G(t)", ls=:auto, lw=2, ylim=(0, 1))
plot!(; legend=:outertopright)

# %%
using LinearAlgebra
using Plots

n = 10
v0 = normalize(randn(ComplexF64, n))
v(t) = exp(im*t)*v0
f(t) = real(dot2(real(v(t))))
g(t) = real(dot2(imag(v(t))))

s0 = sum(z^2 for z in v0)
n0 = norm(v0)^2
F(t) = (n0 + real(exp(im*2t) * s0))/2
G(t) = (n0 - real(exp(im*2t) * s0))/2

t = range(0, 2π; length=201)
P = plot(t, f.(t); label="f(t)", lw=2)
plot!(t, F.(t); label="F(t)", ls=:auto, lw=2, ylim=(0, 1))
plot!(t, g.(t); label="g(t)", lw=2)
plot!(t, G.(t); label="G(t)", ls=:auto, lw=2, ylim=(0, 1))
plot!(; legend=:outertopright)

# %%
