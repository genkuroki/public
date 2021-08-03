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
using Plots
using LinearAlgebra
using Fmt

# %%
f(a) = -a + √(a^2 - 2)
a = range(√2, 100, 1000)
plot(a, f; ylim=(-1.5, 0), label="")

# %%
A(a) = Float64[
     0  1
    -2 -2a
]

U(a, t) = exp(t*A(a))

# %%
@gif for a in [0:0.02:1.4; fill(1.4, 40); 1.4:0.04:5]
    t = range(0, 4π; length=1000)
    u = (t -> U(a, t)[1, 1]).(t)
    v = (t -> U(a, t)[2, 1]).(t)
    P = plot(t, [u v]; label=["u1" "v1"], ylim=(-1, 1))
    u = (t -> U(a, t)[1, 2]).(t)
    v = (t -> U(a, t)[2, 2]).(t)
    plot!(t, [u v]; label=["u2" "v2"], ylim=(-1, 1), ls=:dash)
    plot!(; title=f"a = {$a:.2f}")
end

# %%
for a in 0.1:0.1:2
    t = reverse(range(0, 10π; length=5001))
    n = accumulate(max, norm.(U.(a, t)))
    k = findfirst(t -> norm(U(a, t)) > 0.1, t)
    @show a, t[k], norm(U(a, t[k]))
end

# %%
t = range(0, 4π; length=1000)
u = (t -> U(t)[1, 1]).(t)
v = (t -> U(t)[2, 1]).(t)
plot(t, [u v])

# %%
t = range(0, 10; length=1000)
u = (t -> U(t)[1, 2]).(t)
v = (t -> U(t)[2, 2]).(t)
plot(t, [u v])

# %%
a, V = eigen(A)

# %%
V/0.40824829046386296

# %%
p = [
    -1 + 1im
    2
]

A*p - (-1 - 1im)*p

# %%
q = [
    -1 - 1im
    2
]

A*q - (-1 + 1im)*q

# %%
[minimum(real.(eigvals(Float64[0 1; -2 -2a]))) for a in 0:0.1:2]

# %%
using SymPy

# %%
@vars a

# %%
A = [
    0 1
    -2 -2
]

# %%
A^2//2

# %%
A^3//6

# %%
A^4//24

# %%
B = [
    0 1
    -2 -2a
]

# %%
B^2/factorial(2)

# %%
expand.(B^3)/factorial(3)

# %%
expand.(B^4)/factorial(4)

# %%
