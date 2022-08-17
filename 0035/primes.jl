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
#     display_name: Julia 1.7.3
#     language: julia
#     name: julia-1.7
# ---

# %%
using Primes
using Primes: wheel, wheel_index, wheel_prime

function _primesmask_BitVector(limit::Int)
    limit < 7 && throw(ArgumentError("The condition limit ≥ 7 must be met."))
    n = wheel_index(limit)
    m = wheel_prime(n)
    sieve = trues(n)
    @inbounds for i = 1:wheel_index(isqrt(limit))
        if sieve[i]
            p = wheel_prime(i)
            q = p^2
            j = (i - 1) & 7 + 1
            while q ≤ m
                sieve[wheel_index(q)] = false
                q += wheel[j] * p
                j = j & 7 + 1
            end
        end
    end
    return sieve
end

function my_primes(n)
    list = [2, 3, 5]
    lo, hi = 2, n
    # http://projecteuclid.org/euclid.rmjm/1181070157
    sizehint!(list, 5 + floor(Int, hi / (log(hi) - 1.12) - lo / (log(lo) - 1.12 * (lo > 7))))
    sieve = _primesmask_BitVector(n)
    @inbounds for i = 1:length(sieve)   # don't use eachindex here
        sieve[i] && push!(list, wheel_prime(i))
    end
    list
end

@show _primesmask_BitVector(100)
@show my_primes(100);

# %%
using BenchmarkTools

# %%
P = @btime my_primes(10^5);

# %%
Q = @btime primes(10^5);

# %%
P == Q

# %%
PP = @time primes(10^10)

# %%
function naive_primes(N)
    P = typeof(N)[]
    for n in 2:N
        for i in 2:isqrt(n)
            n % i == 0 && @goto next
        end
        push!(P, n)
        @label next
    end
    P
end

# %%
naive_primes(100) == primes(100)

# %%
@btime naive_primes(10^5)

# %%
naive_primes(10^5) == primes(10^5)

# %%
function naive_first_primes(K)
    P = typeof(K)[]
    k = 0
    n = 2
    while k < K
        for i in 2:isqrt(n)
            n % i == 0 && @goto next
        end
        push!(P, n)
        k += 1
        @label next
        n += 1
    end
    P
end

# %%
using BenchmarkTools

# %%
P = @btime naive_first_primes(10^5)

# %%
using Primes
Q = @btime nextprimes(2, 10^5)
P == Q

# %%
@inline function naive_isprime(n)
    for i in 2:isqrt(n)
        n % i == 0 && return false
    end
    true
end

function naive_primes2(N)
    P = typeof(N)[]
    for n in 2:N
        naive_isprime(n) && push!(P, n)
    end
    P
end

# %%
@btime naive_primes2(10^5)

# %%
naive_primes2(10^5) == primes(10^5)

# %%
