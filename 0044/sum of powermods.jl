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
ENV["LINES"] = 200
ENV["COLUMNS"] = 1000

# %%
n = 5
A = [powermod(x, k, n) for x in 1:n-1, k in 1:n-1]
display(A)
sum(A; dims=1)

# %%
n = 10
A = [powermod(x, k, n) for x in 1:n-1, k in 1:n-1]
display(A)
sum(A; dims=1)

# %%
n = 7
A = [powermod(x, k, n) for x in 1:n-1, k in 1:n-1]
display(A)
sum(A; dims=1)

# %%
n = 14
A = [powermod(x, k, n) for x in 1:n-1, k in 1:n-1]
display(A)
sum(A; dims=1)

# %%
n = 11
A = [powermod(x, k, n) for x in 1:n-1, k in 1:n-1]
display(A)
sum(A; dims=1)

# %%
n = 22
A = [powermod(x, k, n) for x in 1:n-1, k in 1:n-1]
display(A)
sum(A; dims=1)

# %%
n = 13
A = [powermod(x, k, n) for x in 1:n-1, k in 1:n-1]
display(A)
sum(A; dims=1)

# %%
n = 26
A = [powermod(x, k, n) for x in 1:n-1, k in 1:n-1]
display(A)
sum(A; dims=1)

# %%
n = 17
A = [powermod(x, k, n) for x in 1:n-1, k in 1:n-1]
display(A)
sum(A; dims=1)

# %%
n = 34
A = [powermod(x, k, n) for x in 1:n-1, k in 1:n-1]
display(A)
sum(A; dims=1)

# %%
n = 19
A = [powermod(x, k, n) for x in 1:n-1, k in 1:n-1]
display(A)
sum(A; dims=1)

# %%
n = 38
A = [powermod(x, k, n) for x in 1:n-1, k in 1:n-1]
display(A)
sum(A; dims=1)

# %%
n = 29
A = [powermod(x, k, n) for x in 1:n-1, k in 1:n-1]
display(A)
sum(A; dims=1)

# %%
n = 58
A = [powermod(x, k, n) for x in 1:n-1, k in 1:n-1]
display(A)
sum(A; dims=1)

# %%
n = 37
A = [powermod(x, k, n) for x in 1:n-1, k in 1:n-1]
display(A)
sum(A; dims=1)

# %%

# %%
n = 54
A = [powermod(x, k, n) for x in 1:n-1, k in 1:n-1]
display(A)
sum(A; dims=1)

# %%

# %%
@show n = 2^2 + 1
A = [powermod(x, k, n) for x in 1:n-1, k in 1:n-1]
#display(A)
sum(A; dims=1)

# %%
@show n = 2*(2^2 + 1)
A = [powermod(x, k, n) for x in 1:n-1, k in 1:n-1]
#display(A)
sum(A; dims=1)

# %%
@show n = 2*(2^4 + 1)
A = [powermod(x, k, n) for x in 1:n-1, k in 1:n-1]
#display(A)
sum(A; dims=1)

# %%
ENV["LINES"] = 200
ENV["COLUMNS"] = 1000
@show n = 2^8 + 1
A = [powermod(x, k, n) for x in 1:n-1, k in 1:n-1]
#display(A)
sum(A; dims=1)

# %%
ENV["LINES"] = 200
ENV["COLUMNS"] = 1000
@show n = 2*(2^8 + 1)
A = [powermod(x, k, n) for x in 1:n-1, k in 1:n-1]
#display(A)
sum(A; dims=1) |> print

# %%
ENV["LINES"] = 200
ENV["COLUMNS"] = 1000

function f(n)
    A = Vector{Int}(undef, n-1)
    y = fill(1, n-1)
    for k in 1:n-1
        y .*= 1:n-1
        y .= mod.(y, n)
        A[k] = sum(y)
    end
    A
end

# %%
[p*(p-1)÷2 for p in (2, 3, 5, 17, 257, 65537)]

# %%
@show m = 2^0
@show n = m + 1
@time A = f(n)
@show A
@show unique(A[setdiff(1:n-1, [m])])
@show A[m];

# %%
@show m = 2^1
@show n = m + 1
@time A = f(n)
@show A
@show unique(A[setdiff(1:n-1, [m])])
@show A[m];

# %%
@show m = 2^2
@show n = m + 1
@time A = f(n)
@show A
@show unique(A[setdiff(1:n-1, [m])])
@show A[m];

# %%
@show m = 2^4
@show n = m + 1
@time A = f(n)
@show A
@show unique(A[setdiff(1:n-1, [m])])
@show A[m];

# %%
@show m = 2^8
@show n = m + 1
@time A = f(n)
@show A
@show unique(A[setdiff(1:n-1, [m])])
@show A[m];

# %%
@show m = 2^16
@show n = m + 1
@time A = f(n)
@show unique(A[setdiff(1:n-1, [m])])
@show A[m];

# %%
[p + (p-1)*2p for p in (2, 3, 5, 17, 257, 65537)]

# %%
@show m = 2^1
@show n = 2(m + 1)
@time A = f(n)
@show A
@show unique(A[setdiff(1:n-1, [m,2m])])
@show unique(A[[m,2m]]);

# %%
@show m = 2^2
@show n = 2(m + 1)
@time A = f(n)
@show unique(A[setdiff(1:n-1, [m,2m])])
@show unique(A[[m,2m]]);

# %%
@show m = 2^4
@show n = 2(m + 1)
@time A = f(n)
@show unique(A[setdiff(1:n-1, [m,2m])])
@show unique(A[[m,2m]]);

# %%
@show m = 2^8
@show n = 2(m + 1)
@time A = f(n)
@show unique(A[setdiff(1:n-1, [m,2m])])
@show unique(A[[m,2m]]);

# %%
@show m = 2^16
@show n = 2(m + 1)
@time A = f(n)
@show unique(A[setdiff(1:n-1, [m,2m])])
@show unique(A[[m,2m]]);

# %%
using Primes
[(k, 2^2^k + 1, isprime(2^2^k + 1)) for k in big(0):8]

# %%
using Primes
[(k, 2^k + 1, isprime(2^k + 1)) for k in 0:50]

# %%
