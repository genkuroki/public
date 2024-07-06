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
#     display_name: Julia 1.10.4
#     language: julia
#     name: julia-1.10
# ---

# %% [markdown]
# 次のセルのコードは
#
# * https://github.com/JuliaMath/Primes.jl/blob/master/src/Primes.jl
#
# のコードのBitVector版である。一般にVector{Bool}版よりもBitVector版の方が計算は遅くなる。

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
typeof(_primesmask_BitVector(100))

# %%
typeof(Primes._primesmask(100))

# %%
_primesmask_BitVector(10^6) == Primes._primesmask(10^6)

# %%
my_primes(10^6) == primes(10^6)

# %%
@time S_BitVector = _primesmask_BitVector(10^9)

# %%
@time S_VectorBool = Primes._primesmask(10^9)

# %%
S_BitVector == S_VectorBool

# %%
@time P_BitVector = my_primes(10^9)

# %%
@time P_VectorBool = primes(10^9)

# %%
P_BitVector == P_VectorBool

# %% [markdown]
# 次のセルのコードは
#
# * https://github.com/JuliaMath/Primes.jl/blob/master/src/Primes.jl
#
# のコードのBitVector版である。一般にVector{Bool}版よりもBitVector版の方が計算は遅くなる。

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

@show my_primes(100)
@time My_Primes = my_primes(10^9)

# %%
# https://x.com/toku51n/status/1809063593704632457
CodeTest_c = read("CodeTest.c", String)
println("\$ cat CodeTest.c\n")
display("text/markdown", "```C\n" * CodeTest_c * "\n```\n")
println("\n\$ gcc -Wall -O3 -march=native CodeTest.c -o CodeTest.exe")
run(`gcc -Wall -O3 -march=native CodeTest.c -o CodeTest.exe`)
println("\n\$ ./CodeTest.exe")
run(`./CodeTest.exe`)

# %%
function my_primes_naive(limit)
    smax = (limit - 1) ÷ 2
    pSieve = trues(smax)
    @inbounds for i in 1:smax
        if pSieve[i]
            for j in (2i*(i+1)):(2i+1):smax
                pSieve[j] = false
            end
        end
    end
    pPrime = [2]
    sizehint!(pPrime, round(Int, 1.2*limit/log(limit)))
    @inbounds for i in 1:smax
        if pSieve[i]
            push!(pPrime, 2i + 1)
        end
    end
    pPrime
end

@show my_primes_naive(10^8) == my_primes(10^8)

@time My_Primes = my_primes_naive(10^9)
length(My_Primes)

# %%
