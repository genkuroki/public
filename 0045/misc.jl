# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:hydrogen
#     text_representation:
#       extension: .jl
#       format_name: hydrogen
#       format_version: '1.3'
#       jupytext_version: 1.10.3
#   kernelspec:
#     display_name: Julia 1.9.3
#     language: julia
#     name: julia-1.9
# ---

# %%
f(x) = tand(x) / x
for k in  0:8
    x = parse(Float64, "1e-$k")
    @eval @show f($x)
end
@show π/180;

# %%
using LinearAlgebra
using SymPy
@vars s t u λ

# %%
A(x) = Sym[
    cos(x) -sin(x) 0
    sin(x)  cos(x) 0
    0 0 1
]

B(x) = Sym[
    1 0 0
    0 cos(x) -sin(x)
    0 sin(x)  cos(x)
]

X = A(s)*B(t)*A(u)

# %%
det(X - I).simplify()

# %%
p = det(λ*I - X).simplify()

# %%
sympy.Poly(p, λ)

# %%
1

# %%
using Plots
default(fmt=:png)

f(q, z; N=100) = sum(q^n^2 * z^n for n in -N:N)
g(q, z; N=100) = prod((1 - q^2n) * (1 + q^(2n-1)*z) * (1 + q^(2n-1)/z) for n in 1:N)

q = 0.3
z = 0.5
@eval @show f($q, $z) g($q, $z)

plot(z -> f(q, z), 0.1, 10; label="f($q, z)")
plot!(z -> g(q, z), 0.1, 10; label="g($q, z)", ls=:dash)
plot!(xguide="z")

# %%
plot!(xscale=:log10)
_tick = Any[0.1, 0.2, 0.5, 1, 2, 5, 10]
xtick = (float.(_tick), string.(_tick))
plot!(; xtick)

# %%
plot(z -> f(q, z), -10, -0.1; label="f($q, z)")
plot!(z -> g(q, z), -10, -0.1; label="g($q, z)", ls=:dash)
plot!(xguide="z")

# %%
a, b = [1 2; 3 4; 5 6]
a, b

# %%
