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

# %%
# https://x.com/toku51n/status/1809063593704632457
CodeTest_c = read("CodeTest.c", String)
println("\$ cat CodeTest.c\n")
display("text/markdown", "```C\n" * CodeTest_c * "\n```\n")
println("\n\$ gcc -Wall -O3 -march=native CodeTest.c -o CodeTest.exe")
run(`gcc -Wall -O3 -march=native CodeTest.c -o CodeTest.exe`)
println("\n\$ ./CodeTest.exe")
run(`./CodeTest.exe`)
println("\n\$ ./CodeTest.exe")
run(`./CodeTest.exe`)
println("\n\$ ./CodeTest.exe")
run(`./CodeTest.exe`)

# %%
# https://x.com/toku51n/status/1809557670296449179
ideone_usLDXm_c = read("ideone_usLDXm.c", String)
println("\$ cat ideone_usLDXm.c\n")
display("text/markdown", "```C\n" * ideone_usLDXm_c * "\n```\n")
println("\n\$ gcc -Wall -O3 -march=native ideone_usLDXm.c -o ideone_usLDXm.exe")
run(`gcc -Wall -O3 -march=native ideone_usLDXm.c -o ideone_usLDXm.exe`)

println("\n\$ echo 1000000000 | ./ideone_usLDXm.exe")
run(pipeline(`echo 1000000000`, `./ideone_usLDXm.exe`))
println("\n\$ echo 1000000000 | ./ideone_usLDXm.exe")
run(pipeline(`echo 1000000000`, `./ideone_usLDXm.exe`))
println("\n\$ echo 1000000000 | ./ideone_usLDXm.exe")
run(pipeline(`echo 1000000000`, `./ideone_usLDXm.exe`))

println("\n\$ echo 2000000000 | ./ideone_usLDXm.exe")
run(pipeline(`echo 2000000000`, `./ideone_usLDXm.exe`))
println("\n\$ echo 2000000000 | ./ideone_usLDXm.exe")
run(pipeline(`echo 2000000000`, `./ideone_usLDXm.exe`))
println("\n\$ echo 2000000000 | ./ideone_usLDXm.exe")
run(pipeline(`echo 2000000000`, `./ideone_usLDXm.exe`))

# %%
# https://x.com/dc1394/status/1809778241525485747

function my_primes_naive(limit)
    smax = (limit - 1) ÷ 2
    s2max = floor(Int, sqrt(smax)) + 1
    pSieve = trues(smax)
    @inbounds for i in 1:s2max
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

my_primes_naive(100); # compile

@time my_primes_naive(10^9);
@time my_primes_naive(10^9);
My_primes = @time my_primes_naive(10^9);
@show length(My_primes);

@time my_primes_naive(2*10^9);
@time my_primes_naive(2*10^9);
My_primes = @time my_primes_naive(2*10^9);
@show length(My_primes);

# %%
using Primes

primes(100); # compile

@time primes(10^9);
@time primes(10^9);
My_primes = @time primes(10^9);
@show length(My_primes);

@time primes(2*10^9);
@time primes(2*10^9);
My_primes = @time primes(2*10^9);
@show length(My_primes);

# %%
using BenchmarkTools
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

my_primes(100); # compile

@time my_primes(10^9);
@time my_primes(10^9);
My_primes = @time my_primes(10^9);
@show length(My_primes);

@time my_primes(2*10^9);
@time my_primes(2*10^9);
My_primes = @time my_primes(2*10^9);\
@show length(My_primes);

# %%
;"c:/Program Files/Git/usr/bin/diff" -u ideone_usLDXm.c.orig ideone_usLDXm.c

# %%
; gcc --version

# %%
versioninfo()

# %%
N = 10^10
A = @time my_primes_naive(N);
B = @time Primes.primes(N);
@show A == B;
length(A)

# %%
length(@time my_primes(N))

# %%
