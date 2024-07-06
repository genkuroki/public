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
# cf. https://x.com/toku51n/status/1809063593704632457

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

@time My_primes_naive = my_primes_naive(10^9)
length(My_primes_naive)

# %%
using Primes

@time My_primes = primes(10^9)
length(My_primes)

# %%
