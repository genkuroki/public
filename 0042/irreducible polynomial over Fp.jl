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
#     display_name: Julia 1.9.0
#     language: julia
#     name: julia-1.9
# ---

# %%
using Primes
using SymPy
@vars x

# %%
for n in 1:10
    expr = prod(x - big(k) for k in 1:n) + 1
    @show expr.factor()
end

# %%
function allirredpolys(d, p=2, x = symbols("x"))
    d == 1 && return Set(trunc(x+k, p) for k in 0:p-1)
    
    P = Set(trunc(evalpoly(x, (a..., 1)), p) for a in Iterators.product(fill(0:p-1, d)...))
    
    Q = [Set(trunc(evalpoly(x, (a..., 1)) * evalpoly(x, (b..., 1)), p)
            for a in Iterators.product(fill(0:p-1, k)...), b in Iterators.product(fill(0:p-1, d-k)...))
        for k in 1:(d ÷ 2)] |> splat(union)
    setdiff(P, Q)
end

# %%
[allirredpolys(d) for d in 1:6]

# %%
[allirredpolys(d, 3) for d in 1:6]

# %%
f5 = prod(x - k for k in 1:5) + 1 |> expand

# %%
c = f5(x=>0) |> N

# %%
Primes.factor(c)

# %%
sympy.factor(trunc(f5, 2), modulus=2)

# %%
f6 = prod(x - k for k in 1:6) + 1 |> expand

# %%
c = f6(x=>0) |> N

# %%
Primes.factor(c)

# %%
trunc(f6, 7)

# %%
[sympy.factor(trunc(f6, p), modulus=p) for p in primes(103)]

# %%
function isirreducible_slow(f, p=2, x=symbols("x"))
    fx = f(x)
    n = degree(fx, x) |> N
    sympy.rem(x^p^n - x, fx, modulus=p) != 0 && return false
    for t in 2:n
        mod(n, t) == 0 || continue
        Primes.isprime(t) || continue
        gcd(x^p^(n ÷ t) - x, fx, modulus=p) == 1 || return false
    end
    true
end

allpolys(d, p=2, x=symbols("x")) = (trunc(evalpoly(x, (a..., 1)), p) for a in Iterators.product(fill(0:p-1, d)...))

allirredpolys_slow(d, p=2, x=symbols("x")) = [f for f in allpolys(d, p, x) if isirreducible_slow(f, p, x)]

# %%
[allirredpolys_slow(d) |> Set for d in 1:6] == [allirredpolys(d) for d in 1:6]

# %%
[allirredpolys_slow(d, 3) |> Set for d in 1:4] == [allirredpolys(d, 3) for d in 1:4]

# %%
X = @time allirredpolys_slow(5, 3) |> Set
Y = @time allirredpolys(5, 3)
X == Y

# %%
X

# %%
XX = collect(X)

# %%
@time isirreducible_slow(XX[1], 3)

# %%
f5 = prod(x - k for k in 1:5) + 1

# %%
f5 = f5 |> expand

# %%
f5mod2 = trunc(f5, 2)

# %%
sympy.rem(x^2^5 - x, f5mod2, modulus=2)

# %%
sympy.gcd(x^2 - x, f5mod2, modulus=2)

# %% [markdown]
# Thus $f_5\operatorname{mod}2$ is irreducible over $\Bbb{F}_2$.  Therefore $f_5$ is irreducible over $\Bbb{Q}$.

# %%
f6 = prod(x - k for k in 1:6) + 1

# %%
f6 = f6 |> expand

# %%
f6mod7 = trunc(f6, 7)

# %%
Primes.factor(N(f6(0)))

# %% [markdown]
# It follows from Eisenstein's criterion of irreducibility that $f_6$ is irreducible over $\Bbb{Z}$

# %%
f6mod11 = trunc(f6, 11)

# %%
d = 6
p = 11
f6modp = trunc(f6, p)
n = p^d
k = n
expr = Sym(1)
xpow = x
for i in 0:floor(Int, log2(n))
    if isodd(k)
        expr = trunc(sympy.rem(expr * xpow, f6modp), p)
    end
    xpow = trunc(sympy.rem(xpow^2, f6modp), p)
    k = k ÷ 2
end
expr - x

# %%
d = 6
p = 11
f6modp = trunc(f6, p)
n = p^(d ÷ 2)
k = n
expr = Sym(1)
xpow = x
for i in 0:floor(Int, log2(n))
    if isodd(k)
        expr = trunc(sympy.rem(expr * xpow, f6modp), p)
    end
    xpow = trunc(sympy.rem(xpow^2, f6modp), p)
    k = k ÷ 2
end
g1 = trunc(expr - x, p)

# %%
g2 = sympy.rem(f6modp, g1, modulus=p)

# %%
d = 6
p = 11
f6modp = trunc(f6, p)
n = p^(d ÷ 3)
k = n
expr = Sym(1)
xpow = x
for i in 0:floor(Int, log2(n))
    if isodd(k)
        expr = trunc(sympy.rem(expr * xpow, f6modp), p)
    end
    xpow = trunc(sympy.rem(xpow^2, f6modp), p)
    k = k ÷ 2
end
g1 = trunc(expr - x, p)

# %%
g2 = sympy.rem(f6modp, g1, modulus=p)

# %%
g3 = sympy.rem(g1, g2, modulus=p)

# %%
g4 = sympy.rem(g2, g3, modulus=p)

# %%
0*x

# %%
trunc(Sym(-1)+x, 2)

# %%
safetrunc(x::Sym, p::Integer) = x(0) == x ? x : trunc(x, p)

modsym(x::Sym, p::Integer, f::Sym) = sympy.rem(x, f, modulus=p)

function powermodsym(x::Sym, n::Integer, p::Integer, f::Sym)
    k = n
    expr = one(x)
    xpow2 = x
    while true
        if isodd(k)
            expr = modsym(expr * xpow2, p, f)
        end
        xpow2 = modsym(xpow2^2, p, f)
        k = k ÷ 2
        k == 0 && break
    end
    expr
end

function isirreducible_rev(f, p=2, x=symbols("x"))
    fx = f(x)
    n = degree(fx, x) |> N
    sympy.rem(powermodsym(x, p^n, p, fx) - x, fx, modulus=p) != 0 && return false
    for t in 2:n
        mod(n, t) == 0 || continue
        Primes.isprime(t) || continue
        gcd(powermodsym(x, p^(n ÷ t), p, fx) - x, fx, modulus=p) == 1 || return false
    end
    true
end

allpolys(d, p=2, x=symbols("x")) = (trunc(evalpoly(x, (a..., 1)), p) for a in Iterators.product(fill(0:p-1, d)...))

allirredpolys_rev(d, p=2, x=symbols("x")) = [f for f in allpolys(d, p, x) if isirreducible_rev(f, p, x)]

# %%
[allirredpolys_rev(d) |> Set for d in 1:6] == [allirredpolys(d) for d in 1:6]

# %%
[allirredpolys_rev(d, 3) |> Set for d in 1:4] == [allirredpolys(d, 3) for d in 1:4]

# %%
X = @time allirredpolys_rev(5, 3) |> Set
Y = @time allirredpolys(5, 3)
X == Y

# %%
XX = collect(X)
@time isirreducible_slow(XX[1], 3)
@time isirreducible_slow(XX[1], 3)
@time isirreducible_slow(XX[1], 3)
@time isirreducible_rev(XX[1], 3)
@time isirreducible_rev(XX[1], 3)
@time isirreducible_rev(XX[1], 3)

# %%
@time isirreducible_slow(x^7 + x^3 + x + 1, 3)
@time isirreducible_slow(x^7 + x^3 + x + 1, 3)
@time isirreducible_slow(x^7 + x^3 + x + 1, 3)
@time isirreducible_rev(x^7 + x^3 + x + 1, 3)
@time isirreducible_rev(x^7 + x^3 + x + 1, 3)
@time isirreducible_rev(x^7 + x^3 + x + 1, 3)

# %%
f6 = prod(x - k for k in 1:6) + 1

# %%
f6mod11 = trunc(f6, 11)

# %%
# @time isirreducible_slow(f6mod11, 11) # hungs up

# %%
@time isirreducible_rev(f6mod11, 11)

# %%
f5 = prod(x - k for k in 1:5) + 1

# %%
g = trunc(f5(x+1), 2)

# %%
x2 = x^2

# %%
x4 = x2^2

# %%
x8 = sympy.rem(x4^2, g, modulus=2)

# %%
x16 = sympy.rem(x8^2, g, modulus=2)

# %%
x32 = sympy.rem(x16^2, g, modulus=2)

# %%
x32 - x

# %%
g1 = x2 - x

# %%
g2 = sympy.rem(g, g1, modulus=2)

# %%
sympy.rem(g, x^2+x+1, modulus=2)

# %%
