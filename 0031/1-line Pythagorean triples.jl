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
#     display_name: Julia 1.7.3
#     language: julia
#     name: julia-1.7
# ---

# %%
f(n) = [(x, y, z) for x in 1:n for y in x:n for z in y:2n if x^2 + y^2 == z^2 && gcd(x, y) == 1]
@time a = f(1000)
@time a = f(1000)
@time a = f(1000)
print(a)

# %%
g(n) = [(x, y, isqrt(x^2+y^2)) for x in 1:n for y in x:n if (t = x^2 + y^2; t == isqrt(t)^2) && gcd(x, y) == 1]
@time b = g(1000)
@time b = g(1000)
@time b = g(1000)
print(b)

# %%
a == b

# %%
@time c = g(10000)
last(c, 20)

# %%
