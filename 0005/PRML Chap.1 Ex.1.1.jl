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
VERSION

# %%
using Plots
default(size=(400, 250))

using Random
Random.seed!(4649373)

# %%
x = range(0, 1; length=10)

f(x) = sinpi(2x)
noise = 0.2randn(length(x))
t = f.(x) + noise

xs = range(extrema(x)...; length=1000)
plot(xs, f.(xs); label="", xlabel="x", ylabel="t")
scatter!(x, t; label="", color=3)

# %%
PP = []
for d in (0, 1, 3, 9)
    X = x.^(0:d)'
    w = X\t
    g(x) = evalpoly(x, w)
    P = plot(xs, f.(xs); label="")
    plot!(xs, g.(xs); label="degree $d")
    scatter!(x, t; label="")
    push!(PP, P)
end
plot(PP...; size=(800, 500))

# %%
PP = []
for d in (0, 1, 3, 9)
    X = x.^(0:d)'
    w = X'X\X't
    g(x) = evalpoly(x, w)
    P = plot(xs, f.(xs); label="")
    plot!(xs, g.(xs); label="degree $d")
    scatter!(x, t; label="")
    push!(PP, P)
end
plot(PP...; size=(800, 500))

# %%
using SymPy
n = 5
x = collect(symbols("x1:$(n+1)", real=true))

# %%
t = collect(symbols("t1:$(n+1)", real=true))

# %%
d = 3
X = x.^(0:d)'

# %%
A = X'X

# %%
T = X't

# %%
