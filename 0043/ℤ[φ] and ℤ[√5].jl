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
#     display_name: Julia 1.9.2
#     language: julia
#     name: julia-1.9
# ---

# %%
using Plots
default(fmt=:png)

# %%
ϕ = (1 + √5) / 2

# %%
z = [(a+b*ϕ, a-b/ϕ) for a in -100:100 for b in -100:100]
scatter(z; label="")
plot!(xlim=(-10, 10), ylim=(-7, 7))
x, y = range(-15, 15, 300), range(-10, 10, 200)
contour!(x, y, (x, y) -> x*y;
    levels=[-1, 1, -4, 4, -5, 5, -9, 9, -11, 11, -16, 16, -19, 19, -20, 20, -25, 25],
    c=:prism,
    #colorbar=false
)
plot!(aspectratio=1)

# %%
round.(Int, prod.(z)) .|> abs |> sort |> unique |> show

# %%
o = [0, 0]
a = [1, 1]
b = [ϕ, -1/ϕ]
scatter(Tuple.([o, a, b, a+b]))
plot!(aspectratio=1)

# %%
binom(n, k) = binomial(big(n), big(k))
binom(100, 50)

# %%
ε = 1//10
for n in (100, 200, 400, 800)
    Nall = big(2)^n
    N = sum(binom(n, k) for k in 0:n if abs(k//n - 1//2) ≥ ε)
    @show n Nall N N/Nall 1/(4.0n*ε)
    println()
end

# %%
ENV["LINES"] = 100
ENV["COLUMNS"] = 200

# %%
n = 10^4
A = [abs(a^2 - 5b^2) for a in 0:n for b in 0:n] |> sort |> unique
@show A[A .≤ 100];

# %%
z = [(a+b*√5, a-b*√5) for a in -100:100 for b in -100:100]
scatter(z; label="", ms=2, mc=:auto)
plot!(xlim=(-15, 15), ylim=(-10, 10))
x, y = range(-30, 30, 300), range(-20, 20, 200)
contour!(x, y, (x, y) -> x*y;
    levels=[A[A .≤ 100]; -A[A .≤ 100]],
    c=:prism,
    #colorbar=false
)
plot!(aspectratio=1)

# %%
@show A[A .< 1000];

# %%
B = [N for N in 0:5n^2 if !(
                mod(N, 5) == 2 || mod(N, 5) == 3
                || (mod(N, 2) == 0 && mod(N, 4) != 0)
                || (mod(N, 3) == 0 && mod(N, 9) != 0)
                || (mod(N, 7) == 0 && mod(N, 49) != 0)
                || (mod(N, 8) == 0 && mod(N, 16) != 0)
                )]
@show B[B .≤ 100]
@show A[A .≤ 100];

# %%
setdiff(A, B)

# %%
C = setdiff(B, A)
@show C[C .≤ 100];

# %%
@show C[C .≤ 1000];

# %%
p = 2
@show aa = [mod(a^2, p) for a in 0:p-1]
@show bb = [mod(-5b^2, p) for b in 0:p-1];
mod.(aa .+ bb', p) .== 0

# %%
p = 3
@show aa = [mod(a^2, p) for a in 0:p-1]
@show bb = [mod(-5b^2, p) for b in 0:p-1];
mod.(aa .+ bb', p) .== 0

# %%
p = 5
@show aa = [mod(a^2, p) for a in 0:p-1]
@show bb = [mod(-5b^2, p) for b in 0:p-1];
mod.(aa .+ bb', p) .== 0

# %%
p = 7
@show aa = [mod(a^2, p) for a in 0:p-1]
@show bb = [mod(-5b^2, p) for b in 0:p-1];
mod.(aa .+ bb', p) .== 0

# %%
p = 8
@show aa = [mod(a^2, p) for a in 0:p-1]
@show bb = [mod(-5b^2, p) for b in 0:p-1];
mod.(aa .+ bb', p) .== 0

# %%
p = 11
@show aa = [mod(a^2, p) for a in 0:p-1]
@show bb = [mod(-5b^2, p) for b in 0:p-1];
mod.(aa .+ bb', p) .== 0

# %%
p = 13
@show aa = [mod(a^2, p) for a in 0:p-1]
@show bb = [mod(-5b^2, p) for b in 0:p-1];
mod.(aa .+ bb', p) .== 0

# %%
p = 17
@show aa = [mod(a^2, p) for a in 0:p-1]
@show bb = [mod(-5b^2, p) for b in 0:p-1];
mod.(aa .+ bb', p) .== 0

# %%
p = 19
@show aa = [mod(a^2, p) for a in 0:p-1]
@show bb = [mod(-5b^2, p) for b in 0:p-1];
mod.(aa .+ bb', p) .== 0

# %%
p = 23
@show aa = [mod(a^2, p) for a in 0:p-1]
@show bb = [mod(-5b^2, p) for b in 0:p-1];
mod.(aa .+ bb', p) .== 0

# %%
p = 29
@show aa = [mod(a^2, p) for a in 0:p-1]
@show bb = [mod(-5b^2, p) for b in 0:p-1];
mod.(aa .+ bb', p) .== 0

# %%
p = 31
@show aa = [mod(a^2, p) for a in 0:p-1]
@show bb = [mod(-5b^2, p) for b in 0:p-1];
mod.(aa .+ bb', p) .== 0

# %%
p = 37
@show aa = [mod(a^2, p) for a in 0:p-1]
@show bb = [mod(-5b^2, p) for b in 0:p-1];
mod.(aa .+ bb', p) .== 0

# %%
p = 41
@show aa = [mod(a^2, p) for a in 0:p-1]
@show bb = [mod(-5b^2, p) for b in 0:p-1];
mod.(aa .+ bb', p) .== 0

# %%
p = 2
@show aa = [mod(a^2, p) for a in 0:p-1]
@show bb = [mod(-5b^2, p) for b in 0:p-1];
mod.(aa .+ bb', p) .== 0

# %%
p = 4
@show aa = [mod(a^2, p) for a in 0:p-1]
@show bb = [mod(-5b^2, p) for b in 0:p-1];
mod.(aa .+ bb', p) .== 0

# %%
p = 8
@show aa = [mod(a^2, p) for a in 0:p-1]
@show bb = [mod(-5b^2, p) for b in 0:p-1];
mod.(aa .+ bb', p) .== 0

# %%
p = 16
@show aa = [mod(a^2, p) for a in 0:p-1]
@show bb = [mod(-5b^2, p) for b in 0:p-1];
mod.(aa .+ bb', p) .== 0

# %%
